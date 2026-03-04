#include "evoca.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* xorshift32 PRNG — stochastic tie-breaking in reproduction */
static uint32_t g_rng = 0x12345678u;
static inline uint32_t rng_next(void)
{
    g_rng ^= g_rng << 13;
    g_rng ^= g_rng >> 17;
    g_rng ^= g_rng << 5;
    return g_rng;
}

/*
 * Distance-ring classification for each (di,dj) offset, indexed [di+2][dj+2].
 *
 * ring_map: which of the three LUT rings (0-2) the cell belongs to
 *   0 : dist 1   (n1, 4 cells)
 *   1 : dist √2  (n2, 4 cells)
 *   2 : dist 2   (n3, 4 cells)
 *  -1 : centre or outside LUT neighbourhood (skipped)
 *
 * orbit_map: D4 orbit (0-5) for fiducial-pattern matching (full 5×5)
 *   0 : centre                    (1 cell)
 *   1 : dist-1 axis cells         (4 cells)
 *   2 : dist-2 axis cells         (4 cells)
 *   3 : dist-√2 diagonal cells    (4 cells)
 *   4 : dist-2√2 corner cells     (4 cells)
 *   5 : dist-√5 off-axis cells    (8 cells)
 */
static const int ring_map[5][5] = {
    {-1, -1,  2, -1, -1},
    {-1,  1,  0,  1, -1},
    { 2,  0, -1,  0,  2},
    {-1,  1,  0,  1, -1},
    {-1, -1,  2, -1, -1},
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

static inline void lut_flip(uint8_t *b, int bit)
{
    b[bit >> 3] ^= (uint8_t)(1u << (bit & 7));
}

/* FNV-1a hash of a LUT → 32-bit value for genome coloring.
 * Hashes as uint32_t words (LUT_BYTES=32 = 8 words, no remainder). */
static uint32_t lut_hash_fn(const uint8_t *b)
{
    const uint32_t *w = (const uint32_t *)b;
    uint32_t h = 0x811c9dc5u;
    for (int i = 0; i < LUT_BYTES / 4; i++) {
        h ^= w[i];
        h *= 0x01000193u;
    }
    return h;
}

/* Map LUT hash → ARGB color.  Wild-type hash maps to white. */
static int32_t hash_to_color(uint32_t h, uint32_t wt)
{
    if (h == wt) return (int32_t)0xFFFFFFFFu;
    return (int32_t)(0xFF000000u | (h & 0x00FFFFFFu));
}

/* Poisson sample using Knuth's algorithm with the existing PRNG. */
static int poisson_sample(float lambda)
{
    if (lambda < 1e-6f) return 0;
    float L = expf(-lambda);
    float p = 1.0f;
    int k = 0;
    do {
        k++;
        p *= (float)(rng_next() & 0xFFFFFF) / 16777216.0f;
    } while (p > L);
    return k - 1;
}

/* ── Grid state ─────────────────────────────────────────────────── */

static int    gN          = 0;
static float  gfood_inc   = 0.0f;
static float  gm_scale    = 1.0f;
static float  gfood_repro = 0.5f;
static int    ggdiff      = 0;    /* diffusion passes per step */
static float  gmu_lut     = 0.0f; /* per-bit LUT mutation rate */
static float  gmu_cgenom  = 0.0f; /* per-bit cgenom mutation rate */

static uint8_t *v_curr = NULL;   /* [N*N]            */
static uint8_t *v_next = NULL;   /* [N*N]            */
static uint8_t *lut    = NULL;   /* [N*N * LUT_BYTES] */
static uint8_t *cgenom = NULL;   /* [N*N]            */
static float   *f_priv = NULL;   /* [N*N]            */
static float   *F_food = NULL;   /* [N*N]            */
static float   *F_temp = NULL;   /* [N*N] scratch for diffusion */
static uint8_t  *births    = NULL;   /* [N*N] birth events this step */
static uint32_t *lut_color = NULL;   /* [N*N] cached ARGB per cell */
static uint32_t  wt_hash   = 0;     /* hash of wild-type LUT */

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
    F_temp = calloc(cells,               sizeof(float));
    births    = calloc(cells,               sizeof(uint8_t));
    lut_color = calloc(cells,               sizeof(uint32_t));
}

void evoca_free(void)
{
    free(v_curr); v_curr = NULL;
    free(v_next); v_next = NULL;
    free(lut);    lut    = NULL;
    free(cgenom); cgenom = NULL;
    free(f_priv); f_priv = NULL;
    free(F_food); F_food = NULL;
    free(F_temp); F_temp = NULL;
    free(births);    births    = NULL;
    free(lut_color); lut_color = NULL;
    gN = 0;
}

/* ── Metaparam setters ──────────────────────────────────────────── */

void evoca_set_food_inc(float f)   { gfood_inc   = f; }
void evoca_set_m_scale(float m)    { gm_scale    = m; }
void evoca_set_food_repro(float r) { gfood_repro = r; }
void evoca_set_gdiff(int d)        { ggdiff      = d; }
int  evoca_get_gdiff(void)         { return ggdiff;   }
void  evoca_set_mu_lut(float m)    { gmu_lut     = m; }
void  evoca_set_mu_cgenom(float m) { gmu_cgenom  = m; }
float evoca_get_mu_lut(void)       { return gmu_lut;    }
float evoca_get_mu_cgenom(void)    { return gmu_cgenom;  }

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
    wt_hash = lut_hash_fn(lb);
    int32_t wt_color = hash_to_color(wt_hash, wt_hash);  /* white */
    for (size_t i = 0; i < cells; i++)
        lut_color[i] = (uint32_t)wt_color;
}

void evoca_set_lut(int idx, const uint8_t *lb)
{
    memcpy(lut + (size_t)idx * LUT_BYTES, lb, LUT_BYTES);
    lut_color[idx] = (uint32_t)hash_to_color(lut_hash_fn(lb), wt_hash);
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
 * Counts active cells in rings n1 (dist 1), n2 (dist √2), n3 (dist 2),
 * plus reads v_x (current cell state) from v_curr.
 */
static int compute_lut_bit(int row, int col)
{
    int N  = gN;
    int v_x = v_curr[row * N + col];
    int n[3] = {0, 0, 0};

    for (int di = -2; di <= 2; di++) {
        int r = ((row + di) % N + N) % N;
        for (int dj = -2; dj <= 2; dj++) {
            int ring = ring_map[di+2][dj+2];
            if (ring < 0) continue;
            int c = ((col + dj) % N + N) % N;
            if (v_curr[r * N + c]) n[ring]++;
        }
    }
    return LUT_IDX(v_x, n[0], n[1], n[2]);
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

/* ── Food diffusion (3×3 box blur, periodic) ───────────────────── */

static void diffuse_food_once(void)
{
    int N = gN;
    for (int row = 0; row < N; row++) {
        for (int col = 0; col < N; col++) {
            float sum = 0.0f;
            for (int di = -1; di <= 1; di++) {
                int r = ((row + di) % N + N) % N;
                for (int dj = -1; dj <= 1; dj++) {
                    int c = ((col + dj) % N + N) % N;
                    sum += F_food[r * N + c];
                }
            }
            F_temp[row * N + col] = sum * (1.0f / 9.0f);
        }
    }
    /* swap buffers */
    float *tmp = F_food; F_food = F_temp; F_temp = tmp;
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

    /* Phase 2b: Food diffusion (gdiff passes of 3×3 box blur) */
    for (int d = 0; d < ggdiff; d++)
        diffuse_food_once();

    /* Phase 3: Eating */
    for (int row = 0; row < N; row++) {
        for (int col = 0; col < N; col++) {
            int   idx      = row * N + col;
            int   matches  = fiducial_matches(row, col, cgenom[idx]);
            float mouthful = (gm_scale / 25.0f) * matches;
            if (mouthful > F_food[idx]) mouthful = F_food[idx];
            float headroom = 1.0f - f_priv[idx];
            if (mouthful > headroom) mouthful = headroom;
            F_food[idx] -= mouthful;
            f_priv[idx] += mouthful;
        }
    }

    /* Phase 4: Reproduction */
    memset(births, 0, cells);
    for (int row = 0; row < N; row++) {
        for (int col = 0; col < N; col++) {
            int idx = row * N + col;
            if (f_priv[idx] < gfood_repro) continue;

            /* Reservoir-sample the neighbour with minimum f_priv.
             * Tie-breaking is uniformly random to avoid directional drift. */
            int   best_r = -1, best_c = -1;
            float best_f = 1e30f;
            int   n_best = 0;
            for (int di = -1; di <= 1; di++) {
                int r = ((row + di) % N + N) % N;
                for (int dj = -1; dj <= 1; dj++) {
                    if (di == 0 && dj == 0) continue;
                    int c = ((col + dj) % N + N) % N;
                    float f = f_priv[r * N + c];
                    if (f < best_f - 1e-7f) {
                        best_f = f; best_r = r; best_c = c; n_best = 1;
                    } else if (f < best_f + 1e-7f) {
                        if ((rng_next() % (uint32_t)++n_best) == 0) {
                            best_r = r; best_c = c;
                        }
                    }
                }
            }
            int child = best_r * N + best_c;
            memcpy(lut + (size_t)child * LUT_BYTES,
                   lut + (size_t)idx   * LUT_BYTES, LUT_BYTES);
            cgenom[child] = cgenom[idx];
            births[child] = 1;

            /* Mutate child's LUT */
            uint8_t *child_lut = lut + (size_t)child * LUT_BYTES;
            int nf = poisson_sample(gmu_lut * LUT_BITS);
            for (int f = 0; f < nf; f++)
                lut_flip(child_lut, rng_next() % LUT_BITS);

            /* Mutate child's cgenom */
            int nc = poisson_sample(gmu_cgenom * 6);
            for (int f = 0; f < nc; f++)
                cgenom[child] ^= (uint8_t)(1u << (rng_next() % 6));

            /* Update cached genome color (skip expensive hash if unchanged) */
            if (nf > 0)
                lut_color[child] = (uint32_t)hash_to_color(
                    lut_hash_fn(child_lut), wt_hash);
            else
                lut_color[child] = lut_color[idx];

            /* v_curr is dynamical state, not genome — do not copy */
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
                pixels[i] = v_curr[i] ? (int32_t)lut_color[i]
                                      : (int32_t)0xFF000000;
            break;
        case 1:
            for (size_t i = 0; i < cells; i++) {
                float fv = F_food[i]; if (fv > 1.0f) fv = 1.0f;
                uint8_t g = (uint8_t)(fv * 255.0f);
                uint8_t r = v_curr[i] ? 180 : 0;
                pixels[i] = (int32_t)(0xFF000000u
                             | ((uint32_t)r << 16) | ((uint32_t)g << 8));
            }
            break;
        case 2:
            for (size_t i = 0; i < cells; i++) {
                float fv = f_priv[i]; if (fv > 1.0f) fv = 1.0f;
                uint8_t b = (uint8_t)(fv * 255.0f);
                uint8_t r = v_curr[i] ? 180 : 0;
                pixels[i] = (int32_t)(0xFF000000u
                             | ((uint32_t)r << 16) | (uint32_t)b);
            }
            break;
        case 3:
            for (size_t i = 0; i < cells; i++) {
                if (births[i])
                    pixels[i] = v_curr[i] ? (int32_t)0xFFFFFF00u
                                          : (int32_t)0xFF808000u;
                else
                    pixels[i] = v_curr[i] ? (int32_t)0xFF444444u
                                          : (int32_t)0xFF000000u;
            }
            break;
        default:
            for (size_t i = 0; i < cells; i++)
                pixels[i] = (int32_t)0xFF000000;
    }
}

/* ── Accessors ──────────────────────────────────────────────────── */

uint8_t *evoca_get_v(void)      { return v_curr; }
float   *evoca_get_F(void)      { return F_food; }
float   *evoca_get_f(void)      { return f_priv; }
uint8_t *evoca_get_cgenom(void) { return cgenom; }
uint8_t *evoca_get_lut(void)    { return lut;    }
uint8_t *evoca_get_births(void) { return births; }
int      evoca_get_N(void)      { return gN;     }
int      evoca_get_cell_px(void) { return CELL_PX; }
