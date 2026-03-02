#include "evoca.h"
#include <stdlib.h>
#include <string.h>

/*
 * Distance-ring classification for each (di,dj) offset, indexed [di+2][dj+2].
 *
 * ring_map: which of the five rings (0-4) the cell belongs to
 *   0 : dist 1   (n1, 4 cells)
 *   1 : dist √2  (n2, 4 cells)
 *   2 : dist 2   (n3, 4 cells)
 *   3 : dist √5  (n4, 8 cells)
 *   4 : dist 2√2 (n5, 4 cells)
 *  -1 : centre   (ignored for LUT; v_x read separately)
 *
 * orbit_map: D4 orbit (0-5) for fiducial-pattern matching
 *   0 : centre                    (1 cell)
 *   1 : dist-1 axis cells         (4 cells)
 *   2 : dist-2 axis cells         (4 cells)
 *   3 : dist-√2 diagonal cells    (4 cells)
 *   4 : dist-2√2 corner cells     (4 cells)
 *   5 : dist-√5 off-axis cells    (8 cells)
 */
static const int ring_map[5][5] = {
    { 4,  3,  2,  3,  4},
    { 3,  1,  0,  1,  3},
    { 2,  0, -1,  0,  2},
    { 3,  1,  0,  1,  3},
    { 4,  3,  2,  3,  4},
};

static const int orbit_map[5][5] = {
    {4, 5, 2, 5, 4},
    {5, 3, 1, 3, 5},
    {2, 1, 0, 1, 2},
    {5, 3, 1, 3, 5},
    {4, 5, 2, 5, 4},
};

/* ── Bit-packed LUT helpers ─────────────────────────────────────── */

static inline uint8_t lut_get(const uint8_t *b, int bit)
{
    return (b[bit >> 3] >> (bit & 7)) & 1;
}

/* lut_set reserved for future mutation support */
__attribute__((unused))
static inline void lut_set(uint8_t *b, int bit, uint8_t val)
{
    if (val) b[bit >> 3] |=  (uint8_t)(1u << (bit & 7));
    else     b[bit >> 3] &= ~(uint8_t)(1u << (bit & 7));
}

/* ── Grid state ─────────────────────────────────────────────────── */

static int    gN          = 0;
static float  gfood_inc   = 0.0f;
static float  gm_scale    = 1.0f;
static float  gfood_repro = 0.5f;

static uint8_t *v_curr = NULL;   /* [N*N]            */
static uint8_t *v_next = NULL;   /* [N*N]            */
static uint8_t *lut    = NULL;   /* [N*N * LUT_BYTES] */
static uint8_t *cgenom = NULL;   /* [N*N]            */
static float   *f_priv = NULL;   /* [N*N]            */
static float   *F_food = NULL;   /* [N*N]            */

/* ── Lifecycle ──────────────────────────────────────────────────── */

void evoca_init(int N, float food_inc, float m_scale, float food_repro)
{
    evoca_free();
    gN          = N;
    gfood_inc   = food_inc;
    gm_scale    = m_scale;
    gfood_repro = food_repro;

    size_t cells = (size_t)N * N;
    v_curr = calloc(cells,               sizeof(uint8_t));
    v_next = calloc(cells,               sizeof(uint8_t));
    lut    = calloc(cells * LUT_BYTES,   sizeof(uint8_t));
    cgenom = calloc(cells,               sizeof(uint8_t));
    f_priv = calloc(cells,               sizeof(float));
    F_food = calloc(cells,               sizeof(float));
}

void evoca_free(void)
{
    free(v_curr); v_curr = NULL;
    free(v_next); v_next = NULL;
    free(lut);    lut    = NULL;
    free(cgenom); cgenom = NULL;
    free(f_priv); f_priv = NULL;
    free(F_food); F_food = NULL;
    gN = 0;
}

/* ── Metaparam setters ──────────────────────────────────────────── */

void evoca_set_food_inc(float f)   { gfood_inc   = f; }
void evoca_set_m_scale(float m)    { gm_scale    = m; }
void evoca_set_food_repro(float r) { gfood_repro = r; }

/* ── Bulk setters ───────────────────────────────────────────────── */

void evoca_set_v_all(const uint8_t *v, int len)
{
    memcpy(v_curr, v, (size_t)len);
}

void evoca_set_lut_all(const uint8_t *lb)
{
    size_t cells = (size_t)gN * gN;
    for (size_t i = 0; i < cells; i++)
        memcpy(lut + i * LUT_BYTES, lb, LUT_BYTES);
}

void evoca_set_lut(int idx, const uint8_t *lb)
{
    memcpy(lut + (size_t)idx * LUT_BYTES, lb, LUT_BYTES);
}

void evoca_set_cgenom_all(uint8_t cg) { memset(cgenom, cg, (size_t)gN * gN); }

void evoca_set_f_all(float f)
{
    for (size_t i = 0; i < (size_t)gN * gN; i++) f_priv[i] = f;
}

void evoca_set_F_all(float F)
{
    for (size_t i = 0; i < (size_t)gN * gN; i++) F_food[i] = F;
}

/* ── Internal helpers ───────────────────────────────────────────── */

/*
 * Compute the LUT bit index for cell at (row,col).
 * Counts active cells in each of the 5 distance rings separately,
 * plus reads v_x (current cell state) from v_curr.
 */
static int compute_lut_bit(int row, int col)
{
    int N  = gN;
    int v_x = v_curr[row * N + col];
    int n[5] = {0, 0, 0, 0, 0};

    for (int di = -2; di <= 2; di++) {
        int r = ((row + di) % N + N) % N;
        for (int dj = -2; dj <= 2; dj++) {
            int ring = ring_map[di+2][dj+2];
            if (ring < 0) continue;          /* skip centre */
            int c = ((col + dj) % N + N) % N;
            if (v_curr[r * N + c]) n[ring]++;
        }
    }
    return LUT_IDX(v_x, n[0], n[1], n[2], n[3], n[4]);
}

/* Count matches between actual 5×5 config and fiducial pattern. */
static int fiducial_matches(int row, int col, uint8_t cg)
{
    int N = gN, matches = 0;
    for (int di = -2; di <= 2; di++) {
        int r = ((row + di) % N + N) % N;
        for (int dj = -2; dj <= 2; dj++) {
            int c      = ((col + dj) % N + N) % N;
            int orbit  = orbit_map[di+2][dj+2];
            int fid    = (cg >> orbit) & 1;
            if (v_curr[r * N + c] == fid) matches++;
        }
    }
    return matches;
}

/* ── Time step ──────────────────────────────────────────────────── */

void evoca_step(void)
{
    int    N     = gN;
    size_t cells = (size_t)N * N;

    /* Phase 1: CA state update (double-buffered) */
    for (int row = 0; row < N; row++) {
        for (int col = 0; col < N; col++) {
            int idx = row * N + col;
            int bit = compute_lut_bit(row, col);
            v_next[idx] = lut_get(lut + (size_t)idx * LUT_BYTES, bit);
        }
    }
    uint8_t *tmp = v_curr; v_curr = v_next; v_next = tmp;

    /* Phase 2: Environmental food regeneration */
    for (size_t i = 0; i < cells; i++) {
        F_food[i] += gfood_inc;
        if (F_food[i] > 1.0f) F_food[i] = 1.0f;
    }

    /* Phase 3: Eating */
    for (int row = 0; row < N; row++) {
        for (int col = 0; col < N; col++) {
            int   idx      = row * N + col;
            int   matches  = fiducial_matches(row, col, cgenom[idx]);
            float mouthful = (gm_scale / 25.0f) * matches;
            if (mouthful > F_food[idx]) mouthful = F_food[idx];
            F_food[idx] -= mouthful;
            f_priv[idx] += mouthful;
        }
    }

    /* Phase 4: Reproduction */
    for (int row = 0; row < N; row++) {
        for (int col = 0; col < N; col++) {
            int idx = row * N + col;
            if (f_priv[idx] < gfood_repro) continue;

            int   best_r = -1, best_c = -1;
            float best_f = 1e30f;
            for (int di = -1; di <= 1; di++) {
                int r = ((row + di) % N + N) % N;
                for (int dj = -1; dj <= 1; dj++) {
                    if (di == 0 && dj == 0) continue;
                    int c = ((col + dj) % N + N) % N;
                    if (f_priv[r * N + c] < best_f) {
                        best_f = f_priv[r * N + c];
                        best_r = r; best_c = c;
                    }
                }
            }
            int child = best_r * N + best_c;
            memcpy(lut + (size_t)child * LUT_BYTES,
                   lut + (size_t)idx   * LUT_BYTES, LUT_BYTES);
            cgenom[child] = cgenom[idx];
            v_curr[child] = v_curr[idx];
            float half    = f_priv[idx] * 0.5f;
            f_priv[idx]   = half;
            f_priv[child] = half;
        }
    }
}

/* ── Visualisation ──────────────────────────────────────────────── */

void evoca_colorize(int32_t *pixels, int colormode)
{
    size_t cells = (size_t)gN * gN;
    switch (colormode) {
        case 0:
            for (size_t i = 0; i < cells; i++)
                pixels[i] = v_curr[i] ? (int32_t)0xFFFFFFFF
                                      : (int32_t)0xFF000000;
            break;
        case 1:
            for (size_t i = 0; i < cells; i++) {
                uint8_t g = (uint8_t)(F_food[i] * 255.0f);
                uint8_t r = v_curr[i] ? 180 : 0;
                pixels[i] = (int32_t)(0xFF000000u
                             | ((uint32_t)r << 16) | ((uint32_t)g << 8));
            }
            break;
        case 2:
            for (size_t i = 0; i < cells; i++) {
                uint8_t b = (uint8_t)(f_priv[i] * 255.0f);
                uint8_t r = v_curr[i] ? 180 : 0;
                pixels[i] = (int32_t)(0xFF000000u
                             | ((uint32_t)r << 16) | (uint32_t)b);
            }
            break;
        default:
            for (size_t i = 0; i < cells; i++)
                pixels[i] = (int32_t)0xFF000000;
    }
}

/* ── Accessors ──────────────────────────────────────────────────── */

uint8_t *evoca_get_v(void) { return v_curr; }
float   *evoca_get_F(void) { return F_food; }
float   *evoca_get_f(void) { return f_priv; }
int      evoca_get_N(void) { return gN;     }
