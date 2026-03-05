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
static float  gtax        = 0.0f; /* priv food decrement per step */

static uint8_t *v_curr = NULL;   /* [N*N]            */
static uint8_t *v_next = NULL;   /* [N*N]            */
static uint8_t *lut    = NULL;   /* [N*N * LUT_BYTES] */
static uint8_t *cgenom = NULL;   /* [N*N]            */
static float   *f_priv = NULL;   /* [N*N]            */
static float   *F_food = NULL;   /* [N*N]            */
static float   *F_temp = NULL;   /* [N*N] scratch for diffusion */
static uint8_t  *births    = NULL;   /* [N*N] birth events this step */
static uint32_t *lut_color = NULL;   /* [N*N] cached ARGB per cell */
static uint32_t *lut_hash_cache = NULL; /* [N*N] cached FNV-1a hash */
static uint32_t  wt_hash   = 0;     /* hash of wild-type LUT */
static uint32_t  g_step    = 0;     /* global step counter */

/* ── Reproduction age histogram ────────────────────────────────── */

#define REPRO_AGE_MAX 1024          /* bins 0..1023; overflow in last bin */

static uint32_t *last_event_step = NULL;  /* [N*N] step of birth or repro */
static uint32_t  repro_age_hist[REPRO_AGE_MAX]; /* cumulative histogram */
static uint32_t  repro_age_t0 = 0;        /* start accumulating after this step */

/* ── Activity hash table ───────────────────────────────────────── */

#define ACT_INIT_CAP 4096
#define ACT_EMPTY    0        /* key=0 is the empty sentinel */

typedef struct {
    uint64_t activity;    /* cumulative count */
    uint32_t pop_count;   /* current population */
    int32_t  color;       /* ARGB (same as lut_color) */
} act_entry_t;

static uint32_t    *act_keys = NULL;
static act_entry_t *act_vals = NULL;
static int          act_cap  = 0;
static int          act_cnt  = 0;
static int          act_ymax = 2000; /* Y-scale for activity saturation */

static void act_resize(void)
{
    int new_cap = act_cap * 2;
    uint32_t    *nk = calloc((size_t)new_cap, sizeof(uint32_t));
    act_entry_t *nv = calloc((size_t)new_cap, sizeof(act_entry_t));
    for (int i = 0; i < act_cap; i++) {
        if (act_keys[i] == ACT_EMPTY) continue;
        uint32_t slot = act_keys[i] % (uint32_t)new_cap;
        while (nk[slot] != ACT_EMPTY)
            slot = (slot + 1) % (uint32_t)new_cap;
        nk[slot] = act_keys[i];
        nv[slot] = act_vals[i];
    }
    free(act_keys); free(act_vals);
    act_keys = nk; act_vals = nv; act_cap = new_cap;
}

/* Find or insert; returns pointer to value entry. */
static act_entry_t *act_find_or_insert(uint32_t key, int32_t color)
{
    if (act_cnt * 10 >= act_cap * 7) act_resize();
    uint32_t slot = key % (uint32_t)act_cap;
    while (act_keys[slot] != ACT_EMPTY) {
        if (act_keys[slot] == key) return &act_vals[slot];
        slot = (slot + 1) % (uint32_t)act_cap;
    }
    act_keys[slot] = key;
    act_vals[slot].activity  = 0;
    act_vals[slot].pop_count = 0;
    act_vals[slot].color     = color;
    act_cnt++;
    return &act_vals[slot];
}

static void evoca_activity_init(void)
{
    act_cap = ACT_INIT_CAP;
    act_cnt = 0;
    act_keys = calloc((size_t)act_cap, sizeof(uint32_t));
    act_vals = calloc((size_t)act_cap, sizeof(act_entry_t));
}

static void evoca_activity_free(void)
{
    free(act_keys);  act_keys = NULL;
    free(act_vals);  act_vals = NULL;
    act_cap = 0;
    act_cnt = 0;
}

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
    lut_hash_cache = calloc(cells,          sizeof(uint32_t));
    last_event_step = calloc(cells,         sizeof(uint32_t));
    memset(repro_age_hist, 0, sizeof(repro_age_hist));
    g_step = 0;
    evoca_activity_init();
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
    free(lut_hash_cache); lut_hash_cache = NULL;
    free(last_event_step); last_event_step = NULL;
    evoca_activity_free();
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
void  evoca_set_tax(float t)      { gtax       = t; }
float evoca_get_tax(void)         { return gtax;      }

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
    for (size_t i = 0; i < cells; i++) {
        lut_color[i] = (uint32_t)wt_color;
        lut_hash_cache[i] = wt_hash;
    }
}

void evoca_set_lut(int idx, const uint8_t *lb)
{
    memcpy(lut + (size_t)idx * LUT_BYTES, lb, LUT_BYTES);
    uint32_t h = lut_hash_fn(lb);
    lut_color[idx] = (uint32_t)hash_to_color(h, wt_hash);
    lut_hash_cache[idx] = h;
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
    g_step++;

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

    /* Phase 2c: Tax — decrement private food; death if depleted */
    if (gtax > 0.0f) {
        for (size_t i = 0; i < cells; i++) {
            f_priv[i] -= gtax;
            if (f_priv[i] <= 0.0f) {
                f_priv[i] = 0.0f;
                /* Death: zero out LUT genome */
                memset(lut + i * LUT_BYTES, 0, LUT_BYTES);
                uint32_t dh = lut_hash_fn(lut + i * LUT_BYTES);
                lut_color[i] = (uint32_t)hash_to_color(dh, wt_hash);
                lut_hash_cache[i] = dh;
            }
        }
    }

    /* Phase 3: Eating */
    for (int row = 0; row < N; row++) {
        for (int col = 0; col < N; col++) {
            int   idx      = row * N + col;
            int   matches  = fiducial_matches(row, col, cgenom[idx]);
            float mouthful = (gm_scale / 25.0f) * matches * F_food[idx];
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

            /* Record reproduction age and reset timestamps */
            last_event_step[child] = g_step;  /* child just born */
            if (g_step >= repro_age_t0 && last_event_step[idx] >= repro_age_t0) {
                uint32_t age = g_step - last_event_step[idx];
                if (age >= REPRO_AGE_MAX) age = REPRO_AGE_MAX - 1;
                repro_age_hist[age]++;
            }
            last_event_step[idx] = g_step;    /* parent just reproduced */

            memcpy(lut + (size_t)child * LUT_BYTES,
                   lut + (size_t)idx   * LUT_BYTES, LUT_BYTES);
            cgenom[child] = cgenom[idx];
            /* Mutate child's LUT */
            uint8_t *child_lut = lut + (size_t)child * LUT_BYTES;
            int nf = poisson_sample(gmu_lut * LUT_BITS);
            for (int f = 0; f < nf; f++)
                lut_flip(child_lut, rng_next() % LUT_BITS);

            /* Mutate child's cgenom */
            int nc = poisson_sample(gmu_cgenom * 6);
            for (int f = 0; f < nc; f++)
                cgenom[child] ^= (uint8_t)(1u << (rng_next() % 6));

            /* births: 1 = normal, 2 = mutant */
            births[child] = (nf > 0 || nc > 0) ? 2 : 1;

            /* Update cached genome color + hash (skip if unchanged) */
            if (nf > 0) {
                uint32_t ch = lut_hash_fn(child_lut);
                lut_color[child] = (uint32_t)hash_to_color(ch, wt_hash);
                lut_hash_cache[child] = ch;
            } else {
                lut_color[child] = lut_color[idx];
                lut_hash_cache[child] = lut_hash_cache[idx];
            }

            /* v_curr is dynamical state, not genome — do not copy */
            float half    = f_priv[idx] * 0.5f;
            f_priv[idx]   = half;
            f_priv[child] = half;
        }
    }
}

/* ── Activity tracking ─────────────────────────────────────────── */

void evoca_activity_update(void)
{
    if (!act_keys) return;
    size_t cells = (size_t)gN * gN;

    /* Clear all pop_counts */
    for (int i = 0; i < act_cap; i++)
        if (act_keys[i] != ACT_EMPTY)
            act_vals[i].pop_count = 0;

    /* Tally alive cells */
    for (size_t i = 0; i < cells; i++) {
        if (!v_curr[i]) continue;
        act_entry_t *e = act_find_or_insert(lut_hash_cache[i],
                                            (int32_t)lut_color[i]);
        e->pop_count++;
        e->activity++;
    }
}

void evoca_activity_render_col(int32_t *col, int height)
{
    /* Fill with background */
    for (int y = 0; y < height; y++)
        col[y] = (int32_t)0xFF111111u;

    if (!act_keys || act_cnt == 0) return;

    /* Saturation formula: y = (H-1) - (H-1)*act/(act+ymax).
     * Waves rise toward the top as activity >> ymax. */
    uint64_t ymax = (uint64_t)act_ymax;

    /* Per-pixel pop tracking for priority (higher pop overwrites) */
    uint32_t ypop[height];
    memset(ypop, 0, (size_t)height * sizeof(uint32_t));

    /* Pass 1: extinct genomes (pop_count == 0) — dimmed color */
    for (int i = 0; i < act_cap; i++) {
        if (act_keys[i] == ACT_EMPTY) continue;
        if (act_vals[i].pop_count > 0) continue;
        uint64_t act = act_vals[i].activity;
        int y = (height - 1) - (int)((uint64_t)(height - 1) * act / (act + ymax));
        if (y < 0) y = 0;
        if (y >= height) y = height - 1;
        /* Dim: RGB × 0.15 */
        uint32_t c = (uint32_t)act_vals[i].color;
        uint8_t r = (uint8_t)(((c >> 16) & 0xFF) * 15 / 100);
        uint8_t g = (uint8_t)(((c >>  8) & 0xFF) * 15 / 100);
        uint8_t b = (uint8_t)(( c        & 0xFF) * 15 / 100);
        col[y] = (int32_t)(0xFF000000u | ((uint32_t)r << 16)
                           | ((uint32_t)g << 8) | b);
    }

    /* Pass 2: alive genomes — full color, higher pop wins */
    for (int i = 0; i < act_cap; i++) {
        if (act_keys[i] == ACT_EMPTY) continue;
        uint32_t pop = act_vals[i].pop_count;
        if (pop == 0) continue;
        uint64_t act = act_vals[i].activity;
        int y = (height - 1) - (int)((uint64_t)(height - 1) * act / (act + ymax));
        if (y < 0) y = 0;
        if (y >= height) y = height - 1;
        if (pop >= ypop[y]) {
            col[y] = act_vals[i].color;
            ypop[y] = pop;
        }
    }
}

int evoca_activity_get(uint32_t *keys, uint64_t *activities,
                       uint32_t *pop_counts, int32_t *colors, int max_n)
{
    int n = 0;
    if (!act_keys) return 0;
    for (int i = 0; i < act_cap && n < max_n; i++) {
        if (act_keys[i] == ACT_EMPTY) continue;
        keys[n]       = act_keys[i];
        activities[n] = act_vals[i].activity;
        pop_counts[n] = act_vals[i].pop_count;
        colors[n]     = act_vals[i].color;
        n++;
    }
    return n;
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
                if (births[i] == 2)       /* mutant birth: bright magenta */
                    pixels[i] = v_curr[i] ? (int32_t)0xFFFF00FFu
                                          : (int32_t)0xFF800080u;
                else if (births[i] == 1)  /* normal birth: yellow */
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
void     evoca_set_act_ymax(int y) { act_ymax = y > 1 ? y : 1; }
int      evoca_get_act_ymax(void)  { return act_ymax; }
uint32_t *evoca_get_repro_age_hist(void) { return repro_age_hist; }
int      evoca_get_repro_age_max(void)   { return REPRO_AGE_MAX; }
uint32_t evoca_get_step(void)            { return g_step; }
void     evoca_set_repro_age_t0(uint32_t t) { repro_age_t0 = t; }
uint32_t evoca_get_repro_age_t0(void)       { return repro_age_t0; }
void     evoca_reset_repro_age_hist(void)   { memset(repro_age_hist, 0, sizeof(repro_age_hist)); }
