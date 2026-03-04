# Plan: Reduce LUT from 5 rings to 3 rings (n1–n3)

## Context
The 5-ring LUT (n1–n5) is 1407 bytes/cell (~93 MB at N=256), causing slow
hashing and mutation.  Rings n4 (dist √5, 8 cells) and n5 (dist 2√2, 4 cells)
are dropped.  The reduced LUT indexes on (v_x, n1, n2, n3) only:
**2×5×5×5 = 250 bits = 32 bytes/cell** — a 44× reduction.

The FNV-1a hash is also accelerated by hashing 8 × uint32_t words instead
of iterating over individual bytes (LUT_BYTES=32 = 8×4, no remainder).

The 5-ring code is preserved on branch `bigLUT`.

## New constants
```
LUT_BITS  = 250        # 2 × 5 × 5 × 5
LUT_BYTES = 32         # ceil(250 / 8)
Flat bit index: v_x*125 + n1*25 + n2*5 + n3
```

## Files to modify

### C/evoca.h
- `LUT_BITS` → 250, `LUT_BYTES` → 32
- `LUT_IDX(vx,n1,n2,n3)` — 4 params: `(vx)*125 + (n1)*25 + (n2)*5 + (n3)`
- Update header comment: remove n4/n5 from ring table, update sizes/formula
- Update flat bit index comment

### C/evoca.c

**ring_map[5][5]** (line 35): Mark n4 (ring 3) and n5 (ring 4) positions as
-1 so they're skipped. Update comment to say "three rings (0–2)":
```c
static const int ring_map[5][5] = {
    {-1, -1,  2, -1, -1},
    {-1,  1,  0,  1, -1},
    { 2,  0, -1,  0,  2},
    {-1,  1,  0,  1, -1},
    {-1, -1,  2, -1, -1},
};
```

**compute_lut_bit()** (line 207): `int n[3] = {0,0,0}`, call
`LUT_IDX(v_x, n[0], n[1], n[2])`. The 5×5 scan loop is unchanged —
ring_map -1 entries are already skipped.

**lut_hash_fn()** (line 64): Hash as uint32_t words:
```c
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
```

**Mutation** (line 339): Uses `LUT_BITS` macro — picks up new value (250)
automatically. No code change needed.

**All other LUT_BYTES references** (alloc, memcpy, set_lut_all, set_lut,
reproduction copy): pick up new value from macro automatically.

### python/evoca_py.py
- `LUT_BITS = 250`, `LUT_BYTES = 32`
- `_S = (125, 25, 5, 1)` for `(v_x, n1, n2, n3)`
- `lut_bit_index(v_x, n1, n2, n3)` — 4 params, formula `v_x*125 + n1*25 + n2*5 + n3`
- `make_gol_lut()`: 4 nested loops (v_x, n1, n2, n3), no n4/n5. GoL rule
  unchanged: birth iff n1+n2==3 (dead), survive iff n1+n2 ∈ {2,3} (alive),
  same for all n3.
- Update module docstring: remove n4/n5 lines

### Docs/model.md
- Remove n4/n5 rows from ring table
- Update LUT size: 250 bits = 32 bytes
- Update flat bit index formula
- Update ring map diagram to show only 3 rings
- Update memory: 32 + 1 + 1 + 1 + 4 + 4 + 4 + 1 + 4 = **52 bytes/cell**
  (N=256 → ~3.4 MB, N=512 → ~13.6 MB)

### CLAUDE.md
- Remove n4/n5 from ring table
- Update LUT size and flat bit index formula

## Unchanged
- **Fiducial pattern**: cgenom, orbit_map, `fiducial_matches()` still use
  full 5×5 grid. Eating is independent of LUT neighbourhood.
- **controls.py**: no changes needed (mu_lut slider range still applies;
  expected flips = mu_lut × 250 is just smaller)

## Rebuild
```bash
gcc -O2 -Wall -fPIC -dynamiclib -o C/libevoca.dylib C/evoca.c
```

## Verification
1. Build with no warnings
2. GoL glider still advances correctly (period-4 diagonal motion)
3. Large performance improvement (32 bytes/cell vs 1407)
4. mu_lut > 0: colored mutants still appear
5. All 4 color modes + probes still work
