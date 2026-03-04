#ifndef EVOCA_H
#define EVOCA_H

#include <stdint.h>

/*
 * EvoCA — Evolutionary Cellular Automata core
 *
 * Rule space: outer-totalistic on 3-ring neighbourhood (n1, n2, n3).
 *
 * The LUT is indexed by (v_x, n1, n2, n3) where:
 *   v_x ∈ {0,1}   — current cell state
 *   n1  ∈ {0..4}  — active cells at distance 1  (orthogonal Moore nbrs)
 *   n2  ∈ {0..4}  — active cells at distance √2 (diagonal Moore nbrs)
 *   n3  ∈ {0..4}  — active cells at distance 2
 *
 * GoL is exactly encodable: Moore count = n1+n2; GoL ignores n3.
 *
 * LUT size: 2 × 5 × 5 × 5 = 250 bits  → 32 bytes (bit-packed)
 * Flat bit index:  v_x*125 + n1*25 + n2*5 + n3
 *
 * Fiducial pattern c(x): D4-symmetric binary pattern on the 5×5 grid.
 * 6 independent bits (one per D4 orbit) in the lower 6 bits of cgenom.
 * (The fiducial still uses the full 5×5 neighbourhood for eating.)
 */

#define LUT_BITS   250   /* 2*5*5*5 */
#define LUT_BYTES   32   /* ceil(250/8) */

/* Display scale: screen pixels per simulation cell.  Change and recompile. */
#define CELL_PX  2

/* Flat bit index from per-ring counts. */
#define LUT_IDX(vx,n1,n2,n3) \
    ((vx)*125 + (n1)*25 + (n2)*5 + (n3))

/* ── Lifecycle ─────────────────────────────────────────────────────── */

void evoca_init(int N, float food_inc, float m_scale, float food_repro);
void evoca_free(void);

/* ── Metaparam setters ──────────────────────────────────────────────── */

void evoca_set_food_inc(float f);
void evoca_set_m_scale(float m);
void evoca_set_food_repro(float r);
void evoca_set_gdiff(int d);
void evoca_set_mu_lut(float m);
void evoca_set_mu_cgenom(float m);
void evoca_set_tax(float t);

/* ── Bulk setters ───────────────────────────────────────────────────── */

void evoca_set_v_all(const uint8_t *v, int len);

/* Set all cells' LUT from a LUT_BYTES-length bit-packed byte array. */
void evoca_set_lut_all(const uint8_t *lut_bytes);

/* Set one cell's LUT. */
void evoca_set_lut(int idx, const uint8_t *lut_bytes);

void evoca_set_cgenom_all(uint8_t cg);
void evoca_set_f_all(float f);
void evoca_set_F_all(float F);

/* ── Update ─────────────────────────────────────────────────────────── */

void evoca_step(void);

/* ── Visualisation ──────────────────────────────────────────────────── */

/* Fill pixels[N*N] with int32 ARGB values.
   colormode 0: cell state with LUT-hashed genome color (wild-type=white)
   colormode 1: env food F as green; alive cells tinted red
   colormode 2: private food f as blue; alive cells tinted red
   colormode 3: birth events (yellow=birth, magenta=mutant birth, dim=alive, black=dead) */
void evoca_colorize(int32_t *pixels, int colormode);

/* ── Accessors ──────────────────────────────────────────────────────── */

uint8_t *evoca_get_v(void);
float   *evoca_get_F(void);
float   *evoca_get_f(void);
uint8_t *evoca_get_cgenom(void);
uint8_t *evoca_get_lut(void);    /* [N*N * LUT_BYTES] */
uint8_t *evoca_get_births(void); /* [N*N] 0=none, 1=birth, 2=mutant birth */
int      evoca_get_N(void);
int      evoca_get_cell_px(void);
int      evoca_get_gdiff(void);
float    evoca_get_mu_lut(void);
float    evoca_get_mu_cgenom(void);
float    evoca_get_tax(void);

#endif /* EVOCA_H */
