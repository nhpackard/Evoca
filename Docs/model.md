# EvoCA Model Reference

Evolutionary Cellular Automata: a binary 2D CA where every cell carries a
genome governing its local rule and eating behavior, with continuous food
resources and reproduction.

## Table of Contents

- [Grid State](#grid-state)
- [Global Metaparameters](#global-metaparameters)
- [Restricted Mutation](#restricted-mutation)
- [LUT Indexing (Per-Ring Counts)](#lut-indexing-per-ring-counts)
- [Fiducial Pattern c(x) and cgenom](#fiducial-pattern-cx-and-cgenom)
- [Time Step Phases](#time-step-phases)
- [Visualization (Color Modes)](#visualization-color-modes)
- [Python API](#python-api)
- [Activity Tracking](#activity-tracking)
- [Probes](#probes)
- [Reproduction Age Histogram](#reproduction-age-histogram)
- [Controls and Display](#controls-and-display)
- [LUT Helper Functions](#lut-helper-functions)
- [Build](#build)
- [Architecture Notes](#architecture-notes)

---

## Grid State

The simulation is an N x N lattice with periodic (toroidal) boundaries.
Every cell at position **x** carries:

| Field       | Type    | Description                                     |
|-------------|---------|-------------------------------------------------|
| `v(x)`      | uint8   | Binary cell state: 0 (dead) or 1 (alive)        |
| `lut(x)`    | 32 B    | Bit-packed rule LUT (250 bits)                   |
| `cgenom(x)` | uint8   | 6-bit fiducial configuration genome (D4-symmetric) |
| `f(x)`      | float   | Private food store, clamped to [0, 1]            |
| `F(x)`      | float   | Environmental food at this location, clamped to [0, 1] |

---

## Global Metaparameters

All metaparameters can be set at init or adjusted at runtime via sliders.

| Name            | Type  | Default | Range (slider) | Description                                          |
|-----------------|-------|---------|----------------|------------------------------------------------------|
| `food_inc`      | float | 0.0     | [0, 0.5]       | Environmental food added per cell per step            |
| `m_scale`       | float | 1.0     | [0, 10]        | Mouthful scale factor for eating                      |
| `food_repro`    | float | 0.5     | [0, 2]         | Private food threshold triggering reproduction        |
| `gdiff`         | int   | 0       | [0, 10]        | Food diffusion passes (3x3 box blur) per step         |
| `mu_lut`        | float | 0.0     | [0, 0.001]     | Per-bit LUT mutation probability on reproduction      |
| `mu_cgenom`     | float | 0.0     | [0, 0.05]      | Per-bit cgenom mutation probability on reproduction   |
| `tax`           | float | 0.0     | [0, 0.1]       | Private food decrement per step; death if depleted    |
| `restricted_mu` | int   | 0       | checkbox       | If 1, restrict LUT mutations to dynamically active bits |

---

## Restricted Mutation

When `restricted_mu=1`, LUT mutations are restricted to bit positions that
were actually queried during the current CA step.

During Phase 1, a 250-bit mask `lut_active` records which LUT indices were
looked up.  Before Phase 4, the set bits are extracted into an `active_bits`
array (typically 20-50 of 250 bits are active in typical dynamics).

During reproduction, instead of `Poisson(mu_lut * 250)` flips at random
positions, the mutation draws `Poisson(mu_lut * n_active)` flips from only
the active positions.  This ensures every mutation is immediately phenotypic
-- no silent mutations that change the genome hash without affecting dynamics.

When `restricted_mu=0` (default), all 250 bits are eligible (original behavior).

---

## LUT Indexing (Per-Ring Counts)

The LUT maps the local neighborhood configuration to a new cell state.
It is indexed by `(v_x, n1, n2, n3)` where each `nk` is the count of
active (v=1) cells in distance-ring k:

| Ring | Distance | Offsets                            | Cells | Max count |
|------|----------|------------------------------------|-------|-----------|
| n1   | 1        | (+-1,0), (0,+-1)                   | 4     | 4         |
| n2   | sqrt(2)  | (+-1,+-1)                          | 4     | 4         |
| n3   | 2        | (+-2,0), (0,+-2)                   | 4     | 4         |

**Flat bit index**: `v_x*125 + n1*25 + n2*5 + n3`

**Total**: 2 x 5 x 5 x 5 = **250 bits = 32 bytes** per cell.

**Ring map** (indexed by `[di+2][dj+2]`, -1 = centre or outside LUT):

```
-1  -1   2  -1  -1
-1   1   0   1  -1
 2   0  -1   0   2
-1   1   0   1  -1
-1  -1   2  -1  -1
```

GoL (B3/S23) is exactly encodable: it conditions on n1+n2 and ignores n3.

Note: the fiducial pattern for eating still uses the full 5x5
neighbourhood (all 6 D4 orbits).  Only the CA rule LUT is restricted
to 3 rings.

---

## Fiducial Pattern c(x) and cgenom

Each cell's eating behavior is governed by a **fiducial configuration
pattern** — a D4-symmetric 5x5 binary pattern encoded in 6 bits.

### D4 Orbits

The 25 positions of the 5x5 neighborhood fall into 6 orbits under the
D4 symmetry group (reflections about horizontal, vertical, and both
diagonal axes).  Each orbit is controlled by one bit of `cgenom`:

```
 4   5   2   5   4        bit 0: centre               (1 cell)
 5   3   1   3   5        bit 1: axis neighbors        (4 cells)
 2   1   0   1   2        bit 2: axis distance-2       (4 cells)
 5   3   1   3   5        bit 3: diagonal neighbors    (4 cells)
 4   5   2   5   4        bit 4: corners               (4 cells)
                           bit 5: knight-move positions (8 cells)
```

The number at each position is the **orbit index** = the cgenom bit
that controls it.  The pattern value at position (i,j) is:

    pattern[i][j] = (cgenom >> orbit_map[i][j]) & 1

### cgenom Examples

| cgenom     | Binary   | Pattern description                      |
|------------|----------|------------------------------------------|
| `0b000000` | 000000   | All zeros (matches dead cells everywhere) |
| `0b111111` | 111111   | All ones (matches alive cells everywhere) |
| `0b000010` | 000010   | Only axis neighbors = 1 (plus sign)       |
| `0b001010` | 001010   | Axis + diagonal neighbors (3x3 filled)    |
| `0b000001` | 000001   | Only centre = 1                           |
| `0b010100` | 010100   | Corners + axis-2 (ring pattern)           |

To visualize any cgenom value:

```python
from python.evoca_py import cgenom_to_pattern
pat = cgenom_to_pattern(0b001010)
print(pat)
# [[0 0 0 0 0]
#  [0 1 1 1 0]
#  [0 1 0 1 0]
#  [0 1 1 1 0]
#  [0 0 0 0 0]]
```

### How cgenom affects eating

The **fiducial match count** compares the actual cell states in the 5x5
neighborhood against the fiducial pattern:

    matches = sum over all 25 positions of (v_actual == c_fiducial)

The cell's **mouthful** is then:

    M(x) = (m_scale / 25) * matches

With cgenom=0 (all-zero fiducial), matches counts **dead** neighbors.
With cgenom=0b111111 (all-one), matches counts **alive** neighbors.

---

## Time Step Phases

Each call to `evoca_step()` executes five phases in order:

### Phase 1: CA State Update

Double-buffered.  For each cell, count active neighbors per ring,
compute the LUT bit index, look up the new state from the cell's
private LUT.  The queried bit index is OR'd into the `lut_active`
mask (for restricted mutation).  Then swap `v_curr` and `v_next`.

If `restricted_mu` is enabled, the active bit indices are extracted
into `active_bits[]` for use in Phase 4.

### Phase 2: Environmental Food Regeneration

    F(x) += food_inc,   clamped to [0, 1]

### Phase 2b: Food Diffusion

If `gdiff > 0`, perform `gdiff` passes of 3x3 box blur on F:

    F_new(x) = (1/9) * sum of F over 3x3 Moore neighborhood including centre

Each pass uses a scratch buffer and pointer swap (double-buffered).
Periodic boundary conditions.

### Phase 2c: Tax and Death

If `tax > 0`:

    f(x) -= tax,   clamped to 0

If `f(x)` reaches 0, the cell's LUT is zeroed out (all-zero genome:
every neighborhood maps to state 0, so the cell always dies on the
next CA step).  This creates survival pressure: cells must eat enough
to offset the tax or lose their evolved rule.

### Phase 3: Eating

For each cell:
1. Compute fiducial matches (0-25) between 5x5 neighborhood and cgenom
2. `mouthful = (m_scale / 25) * matches * F(x)` (proportional to available food)
3. Clamp to headroom: `mouthful = min(mouthful, 1.0 - f(x))`
4. Transfer: `F(x) -= mouthful`,  `f(x) += mouthful`

Private food f(x) is hard-capped at 1.0.

### Phase 4: Reproduction and Mutation

For each cell where `f(x) >= food_repro`:
1. Find the Moore neighbor (8 cells) with the lowest f(x').
   Ties broken by uniform random (xorshift32 PRNG).
2. Copy parent genome to child: LUT, cgenom.
   (v_curr is dynamical state, NOT copied.)
3. **Mutate child's LUT**: if `restricted_mu`, draw n_flips ~
   Poisson(mu_lut * n_active) and flip only active bits; otherwise
   draw n_flips ~ Poisson(mu_lut * 250) and flip random bits.
4. **Mutate child's cgenom**: draw n_flips ~ Poisson(mu_cgenom * 6),
   flip that many random bits in the child's cgenom.
5. Update child's cached genome color (FNV-1a hash of LUT).
6. Split food: `f(parent) = f(child) = f(parent) / 2`

The `births` array is set for each child cell (used by colormode 3).

---

## Visualization (Color Modes)

Five modes, selectable via dropdown or the `colormode` parameter:

| Mode | Name       | Channel mapping (ARGB)                                        |
|------|------------|---------------------------------------------------------------|
| 0    | `state`    | alive = LUT-hashed genome color (wild-type = white), dead = black |
| 1    | `env-food` | Green = F(x)*255, Red = 180 if alive else 0                   |
| 2    | `priv-food`| Blue = f(x)*255, Red = 180 if alive else 0                    |
| 3    | `births`   | birth+alive = bright yellow, birth+dead = dim yellow, alive = dim grey, dead = black |

All food values are clamped to [0, 1] before the float-to-uint8 cast.

### Genome coloring (mode 0)

Each cell caches an ARGB color derived from an FNV-1a hash of its 32-byte
LUT.  The hash computed by `set_lut_all()` is stored as the **wild-type hash**;
any cell whose LUT hashes to this value displays as white.  Mutant genomes
get pseudo-random colors from the lower 24 bits of their hash.

---

## Activity Tracking

### LUT Activity

Tracks cumulative presence of each distinct LUT genome over time.
Each genome is identified by its FNV-1a hash (cached in `lut_hash_cache[N*N]`).
An open-addressing hash table maps hash -> `{activity, pop_count, color}`.

- `evoca_activity_update()`: clears pop_counts, scans alive cells, increments activity
- `evoca_activity_render_col(col, height)`: renders one column of the scrolling
  strip chart. Alive genomes in full color; extinct genomes dimmed (RGB x 0.15).

**Y-axis saturation formula** (from genelife):

    y = (H-1) - (H-1) * act / (act + ymax)

This is a hyperbolic saturation curve: `act=0` maps to the bottom, `act=ymax`
maps to mid-chart, and `act->inf` approaches the top.  Tunable via `act_ymax`
slider (default 2000).

Enable with `probes={'activity': True}`.

### Cgenom Activity

Mirrors LUT activity for the 6-bit cgenom (fiducial eating pattern).  Since
there are only 2^6 = 64 possible cgenoms, uses fixed-size arrays instead of
a hash table.  Wild-type cgenom is colored white; mutants get FNV-1a hash
colors.  Separate `cg_act_ymax` slider.

Enable with `probes={'cg_activity': True}`.

---

## Probes

Probes are optional strip-chart windows enabled via the `probes` dict
parameter to `run_with_controls()`.  Each probe renders in its own SDL
window, stacked to the left of the main window.

| Probe name       | Window size   | What it shows                                          |
|------------------|---------------|--------------------------------------------------------|
| `env_food`       | 512 x 128     | Mean +/- std of environmental food F(x) over time      |
| `priv_food`      | 512 x 128     | Mean +/- std of private food f(x) over time            |
| `births`         | 512 x 128     | Mean +/- std of births array over time                 |
| `activity`       | 512 x 256     | LUT genome activity (scrolling hash-colored strip)     |
| `cg_activity`    | 512 x 256     | Cgenom activity (scrolling hash-colored strip)         |
| `lut_complexity` | 512 x 128     | Stacked area: green=n1 only, yellow=n1+n2, red=n1+n2+n3 |

Example enabling multiple probes:

```python
run_with_controls(sim, probes={
    'activity': True,
    'cg_activity': True,
    'lut_complexity': True,
    'env_food': True,
})
```

### LUT Complexity Probe

Classifies each alive cell's LUT by the minimum ring set it depends on:

- **Level 1** (green): rule depends only on (v_x, n1) — constant across n2 and n3
- **Level 2** (yellow): rule depends on (v_x, n1, n2) — constant across n3
- **Level 3** (red): rule depends on all three rings (v_x, n1, n2, n3)

Dead cells are excluded from counts.  The stacked area chart shows population
fractions at each complexity level over time.  Useful for studying whether
evolution discovers higher-ring dependencies starting from simple n1-only rules.

---

## Reproduction Age Histogram

Tracks the distribution of time between successive reproduction events
(or birth-to-first-reproduction).  A per-cell timestamp `last_event_step[N*N]`
records each cell's most recent birth or reproduction.  At each reproduction,
`age = step - last_event_step[parent]` is binned into `repro_age_hist[1024]`.

A configurable `repro_age_t0` (default 0) skips the transient: only events
where both the current step and the parent's last event are >= t0 are counted.

```python
sim.set_repro_age_t0(1000)     # start accumulating after step 1000
sim.reset_repro_age_hist()     # clear histogram
hist = sim.get_repro_age_hist() # numpy array, 1024 bins
```

---

## Python API

### EvoCA class  (`python/evoca_py.py`)

#### Lifecycle

```python
sim = EvoCA(lib_path=None)
    # Load the shared library (auto-finds C/libevoca.dylib or .so)

sim.init(N, food_inc=0.0, m_scale=1.0, food_repro=0.5, gdiff=0,
         mu_lut=0.0, mu_cgenom=0.0, tax=0.0, restricted_mu=0)
    # Allocate N x N lattice, set metaparameters.
    # All grids initialized to zero.

sim.free()
    # Deallocate C-side memory.
```

#### Metaparameter Setters

```python
sim.update_food_inc(f)       # float
sim.update_m_scale(m)        # float
sim.update_food_repro(r)     # float
sim.update_gdiff(d)          # int
sim.update_mu_lut(m)         # float, per-bit LUT mutation rate
sim.update_mu_cgenom(m)      # float, per-bit cgenom mutation rate
sim.update_tax(t)            # float, priv food decrement per step
sim.update_restricted_mu(r)  # int (0 or 1)
sim.update_act_ymax(y)       # int, Y-scale for LUT activity chart
sim.update_cg_act_ymax(y)    # int, Y-scale for cgenom activity chart
```

Each setter updates both the Python attribute (`sim.food_inc`, etc.)
and the C global.

#### Grid Setters

```python
sim.set_v(v_array)
    # Set cell states from (N, N) or flat uint8 array.

sim.set_lut_all(lut_bytes)
    # Set ALL cells' LUT from a single LUT_BYTES-length uint8 array.

sim.set_lut(idx, lut_bytes)
    # Set one cell's LUT (idx = flat cell index).

sim.set_cgenom_all(cg)
    # Set all cells' fiducial genome (6-bit value, masked to 0x3F).

sim.set_cgenom_random()
    # Set each cell's cgenom to a random value in [0, 63].
    # No wild-type: all 64 cgenoms get distinct hash-based colors.

sim.set_lut_random(n_init=3)
    # Set each cell's LUT to an independent random rule.
    # n_init: number of rings the rule conditions on (1, 2, or 3).
    #   1: depends on (v_x, n1) only — 10 independent bits per LUT
    #   2: depends on (v_x, n1, n2) — 50 independent bits
    #   3: depends on (v_x, n1, n2, n3) — all 250 bits independent

sim.set_f_all(f)
    # Set all cells' private food to float f.

sim.set_F_all(F)
    # Set all cells' environmental food to float F.

sim.set_F_random(lo=0.0, hi=1.0)
    # Set env food to uniform random floats in [lo, hi].
```

#### Step and Colorize

```python
sim.step()
    # Execute one time step (all 5 phases).

sim.colorize(pixels, colormode=0)
    # Fill a (N*N,) int32 numpy array with ARGB pixel values.
    # colormode: 0=state, 1=env-food, 2=priv-food, 3=births
```

#### Getters

```python
sim.get_v()       -> np.ndarray   # (N, N) uint8, cell states (copy)
sim.get_F()       -> np.ndarray   # (N, N) float32, env food (copy)
sim.get_f()       -> np.ndarray   # (N, N) float32, private food (copy)
sim.get_cgenom()  -> np.ndarray   # (N, N) uint8, fiducial genomes (copy)
sim.get_births()  -> np.ndarray   # (N, N) uint8, birth events last step (copy)
sim.get_lut(idx)  -> np.ndarray   # (LUT_BYTES,) uint8, one cell's LUT (copy)
sim.get_step()    -> int           # current global step counter
sim.N             -> int           # lattice size
sim.cell_px       -> int           # CELL_PX compile-time constant

# Activity and diagnostics
sim.get_activity(max_n=4096) -> dict   # {'hash', 'activity', 'pop_count', 'color'}
sim.get_cg_activity()        -> dict   # {'activity', 'pop_count', 'color'} (64 entries)
sim.get_lut_complexity()     -> dict   # {'n1': count, 'n2': count, 'n3': count}
sim.get_repro_age_hist()     -> np.ndarray  # (1024,) uint32 histogram
sim.set_repro_age_t0(t)                # set step threshold for histogram accumulation
sim.reset_repro_age_hist()             # clear histogram

# Params export
sim.params()      -> dict         # current metaparameters as dict
sim.params_str()  -> str          # copy-pasteable sim.init(...) call
```

#### Stored Attributes

After `init()` or the corresponding setter, these Python attributes
reflect current values:

    sim.food_inc, sim.m_scale, sim.food_repro, sim.gdiff,
    sim.mu_lut, sim.mu_cgenom, sim.tax, sim.restricted_mu, sim.cgenom

---

## Controls and Display

### run_with_controls  (`python/controls.py`)

```python
run_with_controls(sim, cell_px=None, colormode=0, paused=False, probes=None)
    -> threading.Thread
```

Opens an SDL2 window and displays ipywidgets controls below the
notebook cell.  Returns immediately (non-blocking).

**Parameters**:

| Parameter   | Type  | Default        | Description                          |
|-------------|-------|----------------|--------------------------------------|
| `sim`       | EvoCA | (required)     | Initialized EvoCA instance           |
| `cell_px`   | int   | `sim.cell_px`  | Screen pixels per simulation cell    |
| `colormode` | int   | 0              | Initial color mode (0/1/2/3)         |
| `paused`    | bool  | False          | Start in paused state                |
| `probes`    | dict  | None           | Probe names to enable (see [Probes](#probes) for full list) |

**Widgets**:
- **Pause/Run** toggle button
- **Step** button (advances one step when paused)
- **Quit** button (closes SDL window and stops simulation)
- **food_inc** slider: [0, 0.5], step 0.001
- **m_scale** slider: [0, 10], step 0.1
- **food_repro** slider: [0, 2], step 0.05
- **gdiff** slider: [0, 10], step 1
- **mu_lut** slider: [0, 0.001], step 0.00001
- **mu_cgenom** slider: [0, 0.05], step 0.001
- **tax** slider: [0, 0.1], step 0.001
- **act_ymax** slider: [100, 100000] (when activity probe enabled)
- **cg_act_ymax** slider: [100, 100000] (when cg_activity probe enabled)
- **restricted_mu** checkbox: toggle restricted mutation
- **Color** dropdown: state / env-food / priv-food / births
- **Save Plots** button: saves probe strip charts to PNG (when paused)
- **Export Params** button: prints copy-pasteable `sim.init(...)` call
- **Fiducial pattern**: 5x5 matplotlib grid displayed above widgets

**SDL2 window title** shows: time step, FPS (when running), color mode,
PAUSED indicator.

**Keyboard** (SDL2 window): Q or Esc to quit.

**Slider behavior**: dragging any slider auto-pauses the simulation;
200 ms after the last touch, it auto-resumes (unless already paused).

### Architecture

SDL2 requires the main thread on macOS.  Since Jupyter's kernel runs
cells in a worker thread, SDL2 is launched in a **subprocess**
(`python/sdl_worker.py`).  Pixel data and control signals flow via
POSIX shared memory (zero-copy):

```
Main process (Jupyter kernel)
+-- Jupyter event loop    -> ipywidgets callbacks
+-- Simulation thread     -> sim.step() + sim.colorize() -> pixel_shm
+-- Reader thread         -> relays SDL subprocess stdout
+-- subprocess (sdl_worker.py)
    +-- SDL2 main thread  -> reads pixel_shm, renders window
```

**ctrl_shm layout** (5 x int32):

| Index | Name      | Values                        |
|-------|-----------|-------------------------------|
| 0     | quit      | 0 = running, 1 = exit         |
| 1     | colormode | 0 / 1 / 2                     |
| 2     | step_cnt  | current time step             |
| 3     | fps * 10  | FPS scaled for int storage     |
| 4     | paused    | 0 = running, 1 = paused        |

---

## LUT Helper Functions

### make_gol_lut  (`python/evoca_py.py`)

```python
make_gol_lut() -> np.ndarray   # LUT_BYTES-length uint8
```

Builds Conway's Game of Life (B3/S23) as a bit-packed LUT.
- Dead cell (v_x=0): birth iff Moore count n1+n2 == 3
- Alive cell (v_x=1): survive iff Moore count n1+n2 in {2, 3}
- n3 is a don't-care (all values set identically).

With mutation=0 (all cells share this LUT), the simulation runs exact GoL.

### pack_lut / unpack_lut

```python
pack_lut(bits) -> np.ndarray
    # Pack LUT_BITS-length uint8 0/1 array into LUT_BYTES bytes.

unpack_lut(packed) -> np.ndarray
    # Unpack LUT_BYTES bytes into LUT_BITS-length uint8 0/1 array.
```

### cgenom_to_pattern

```python
cgenom_to_pattern(cg) -> np.ndarray   # (5, 5) uint8
    # Expand a 6-bit cgenom value to a 5x5 binary array.
```

### lut_bit_index

```python
lut_bit_index(v_x, n1, n2, n3) -> int
    # Compute the flat bit index into the LUT.
```

### Constants

```python
LUT_BITS  = 250      # bits per LUT
LUT_BYTES = 32       # bytes per LUT (ceil(250/8))
ORBIT_MAP            # (5, 5) uint8 array of orbit indices
```

---

## Build

```bash
# macOS
gcc -O2 -Wall -fPIC -dynamiclib -o C/libevoca.dylib C/evoca.c

# Linux
gcc -O2 -Wall -fPIC -shared -o C/libevoca.so C/evoca.c

# or just:
make
```

The compile-time constant `CELL_PX` (default 2) in `C/evoca.h` sets
display scaling.  Change it and recompile.

---

## Architecture Notes

### Performance (CA step only, no food/repro, GoL LUT)

- N=256: ~150 fps
- N=512: ~30 fps

### Memory

Per cell: 32 (LUT) + 1 (cgenom) + 1 (v_curr) + 1 (v_next) + 4 (f_priv)
+ 4 (F_food) + 4 (F_temp) + 1 (births) + 4 (lut_color) + 4 (lut_hash_cache)
+ 4 (last_event_step) = **60 bytes**.

- N=256: ~3.4 MB
- N=512: ~13.6 MB

### SDL2 on macOS

SDL2 video init (Cocoa/NSApplication) requires the actual main thread
(thread 0).  Jupyter's ipykernel runs sync cells in an executor thread,
so threading.Thread for SDL2 crashes the kernel.  The fix: run SDL2 in a
subprocess, which has its own main thread.
