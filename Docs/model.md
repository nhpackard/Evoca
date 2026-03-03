# EvoCA Model Reference

Evolutionary Cellular Automata: a binary 2D CA where every cell carries a
genome governing its local rule and eating behavior, with continuous food
resources and reproduction.

## Table of Contents

- [Grid State](#grid-state)
- [Global Metaparameters](#global-metaparameters)
- [LUT Indexing (Per-Ring Counts)](#lut-indexing-per-ring-counts)
- [Fiducial Pattern c(x) and cgenom](#fiducial-pattern-cx-and-cgenom)
- [Time Step Phases](#time-step-phases)
- [Visualization (Color Modes)](#visualization-color-modes)
- [Python API](#python-api)
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
| `lut(x)`    | 1407 B  | Bit-packed rule LUT (11250 bits)                 |
| `cgenom(x)` | uint8   | 6-bit fiducial configuration genome (D4-symmetric) |
| `f(x)`      | float   | Private food store, clamped to [0, 1]            |
| `F(x)`      | float   | Environmental food at this location, clamped to [0, 1] |

---

## Global Metaparameters

All metaparameters can be set at init or adjusted at runtime via sliders.

| Name         | Type  | Default | Range (slider) | Description                                       |
|--------------|-------|---------|----------------|---------------------------------------------------|
| `food_inc`   | float | 0.0     | [0, 0.5]       | Environmental food added per cell per step         |
| `m_scale`    | float | 1.0     | [0, 10]        | Mouthful scale factor for eating                   |
| `food_repro` | float | 0.5     | [0, 2]         | Private food threshold triggering reproduction     |
| `gdiff`      | int   | 0       | [0, 10]        | Food diffusion passes (3x3 box blur) per step      |

---

## LUT Indexing (Per-Ring Counts)

The LUT maps the local neighborhood configuration to a new cell state.
It is indexed by `(v_x, n1, n2, n3, n4, n5)` where each `nk` is the
count of active (v=1) cells in distance-ring k:

| Ring | Distance | Offsets                            | Cells | Max count |
|------|----------|------------------------------------|-------|-----------|
| n1   | 1        | (+-1,0), (0,+-1)                   | 4     | 4         |
| n2   | sqrt(2)  | (+-1,+-1)                          | 4     | 4         |
| n3   | 2        | (+-2,0), (0,+-2)                   | 4     | 4         |
| n4   | sqrt(5)  | (+-2,+-1), (+-1,+-2)               | 8     | 8         |
| n5   | 2*sqrt(2)| (+-2,+-2)                          | 4     | 4         |

**Flat bit index**: `v_x*5625 + n1*1125 + n2*225 + n3*45 + n4*5 + n5`

**Total**: 2 x 5 x 5 x 5 x 9 x 5 = **11,250 bits = 1,407 bytes** per cell.

**Ring map** (indexed by `[di+2][dj+2]`, -1 = centre):

```
 4   3   2   3   4
 3   1   0   1   3
 2   0  -1   0   2
 3   1   0   1   3
 4   3   2   3   4
```

**Design rationale**: The original spec proposed a weighted sum
S = A + B*sqrt(2) + C*sqrt(5) with A = n1 + 2*n3.  This conflates n1
and n3: the same A can arise from different Moore counts, making GoL
unrepresentable.  Separate per-ring counts remove this ambiguity.
GoL (B3/S23) is exactly encodable: it conditions on n1+n2 and ignores
n3, n4, n5.

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
private LUT.  Then swap `v_curr` and `v_next`.

### Phase 2: Environmental Food Regeneration

    F(x) += food_inc,   clamped to [0, 1]

### Phase 2b: Food Diffusion

If `gdiff > 0`, perform `gdiff` passes of 3x3 box blur on F:

    F_new(x) = (1/9) * sum of F over 3x3 Moore neighborhood including centre

Each pass uses a scratch buffer and pointer swap (double-buffered).
Periodic boundary conditions.

### Phase 3: Eating

For each cell:
1. Compute fiducial matches (0-25) between 5x5 neighborhood and cgenom
2. `mouthful = (m_scale / 25) * matches`
3. Clamp to available food: `mouthful = min(mouthful, F(x))`
4. Clamp to headroom: `mouthful = min(mouthful, 1.0 - f(x))`
5. Transfer: `F(x) -= mouthful`,  `f(x) += mouthful`

Private food f(x) is hard-capped at 1.0.

### Phase 4: Reproduction

For each cell where `f(x) >= food_repro`:
1. Find the Moore neighbor (8 cells) with the lowest f(x').
   Ties broken by uniform random (xorshift32 PRNG).
2. Copy parent genome to child: LUT, cgenom, v_curr.
3. Split food: `f(parent) = f(child) = f(parent) / 2`

---

## Visualization (Color Modes)

Three modes, selectable via dropdown or the `colormode` parameter:

| Mode | Name       | Channel mapping (ARGB)                          |
|------|------------|-------------------------------------------------|
| 0    | `state`    | alive = white (0xFFFFFFFF), dead = black         |
| 1    | `env-food` | Green = F(x)*255, Red = 180 if alive else 0      |
| 2    | `priv-food`| Blue = f(x)*255, Red = 180 if alive else 0       |

All food values are clamped to [0, 1] before the float-to-uint8 cast.

---

## Python API

### EvoCA class  (`python/evoca_py.py`)

#### Lifecycle

```python
sim = EvoCA(lib_path=None)
    # Load the shared library (auto-finds C/libevoca.dylib or .so)

sim.init(N, food_inc=0.0, m_scale=1.0, food_repro=0.5, gdiff=0)
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
    # colormode: 0=state, 1=env-food, 2=priv-food
```

#### Getters

```python
sim.get_v()   -> np.ndarray   # (N, N) uint8, cell states (copy)
sim.get_F()   -> np.ndarray   # (N, N) float32, env food (copy)
sim.get_f()   -> np.ndarray   # (N, N) float32, private food (copy)
sim.N         -> int           # lattice size
sim.cell_px   -> int           # CELL_PX compile-time constant
```

#### Stored Attributes

After `init()` or the corresponding setter, these Python attributes
reflect current values:

    sim.food_inc, sim.m_scale, sim.food_repro, sim.gdiff, sim.cgenom

---

## Controls and Display

### run_with_controls  (`python/controls.py`)

```python
run_with_controls(sim, cell_px=None, colormode=0, paused=False)
    -> threading.Thread
```

Opens an SDL2 window and displays ipywidgets controls below the
notebook cell.  Returns immediately (non-blocking).

**Parameters**:

| Parameter   | Type  | Default        | Description                          |
|-------------|-------|----------------|--------------------------------------|
| `sim`       | EvoCA | (required)     | Initialized EvoCA instance           |
| `cell_px`   | int   | `sim.cell_px`  | Screen pixels per simulation cell    |
| `colormode` | int   | 0              | Initial color mode (0/1/2)           |
| `paused`    | bool  | False          | Start in paused state                |

**Widgets**:
- **Pause/Run** toggle button
- **Step** button (advances one step when paused)
- **Quit** button (closes SDL window and stops simulation)
- **food_inc** slider: [0, 0.5], step 0.001
- **m_scale** slider: [0, 10], step 0.1
- **food_repro** slider: [0, 2], step 0.05
- **gdiff** slider: [0, 10], step 1
- **Color** dropdown: state / env-food / priv-food
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
- Outer rings (n3, n4, n5) are don't-cares (all combinations set).

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
lut_bit_index(v_x, n1, n2, n3, n4, n5) -> int
    # Compute the flat bit index into the LUT.
```

### Constants

```python
LUT_BITS  = 11250    # bits per LUT
LUT_BYTES = 1407     # bytes per LUT (ceil(11250/8))
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

Per cell: 1407 (LUT) + 1 (cgenom) + 1 (v_curr) + 1 (v_next) + 4 (f_priv)
+ 4 (F_food) + 4 (F_temp) = **1422 bytes**.

- N=256: ~93 MB
- N=512: ~373 MB

### SDL2 on macOS

SDL2 video init (Cocoa/NSApplication) requires the actual main thread
(thread 0).  Jupyter's ipykernel runs sync cells in an executor thread,
so threading.Thread for SDL2 crashes the kernel.  The fix: run SDL2 in a
subprocess, which has its own main thread.
