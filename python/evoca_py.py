"""
evoca_py.py — Python ctypes wrapper for the EvoCA C library.

LUT indexing
------------
The LUT conditions on (v_x, n1, n2, n3):
  v_x ∈ {0,1}  — current cell state
  n1  ∈ {0..4} — active cells at distance 1   (4 orthogonal Moore nbrs)
  n2  ∈ {0..4} — active cells at distance √2  (4 diagonal Moore nbrs)
  n3  ∈ {0..4} — active cells at distance 2   (4 axis-aligned)

Flat bit index:  v_x*125 + n1*25 + n2*5 + n3
Total bits:  2*5*5*5 = 250  →  32 bytes bit-packed per cell.

GoL is exactly encodable: it conditions on n1+n2 and ignores n3.
"""

import ctypes
import os

import numpy as np

# ── LUT constants ──────────────────────────────────────────────────────
LUT_BITS  = 250   # 2*5*5*5
LUT_BYTES = 32    # ceil(250/8)

# Ring strides for flat index
_S = (125, 25, 5, 1)   # (v_x, n1, n2, n3)


def lut_bit_index(v_x, n1, n2, n3):
    return v_x*125 + n1*25 + n2*5 + n3


# D4-symmetric orbit map for the 5×5 fiducial pattern.
# 6 independent orbits (bits 0–5 of cgenom).
ORBIT_MAP = np.array([
    [4, 5, 2, 5, 4],
    [5, 3, 1, 3, 5],
    [2, 1, 0, 1, 2],
    [5, 3, 1, 3, 5],
    [4, 5, 2, 5, 4],
], dtype=np.uint8)


def cgenom_to_pattern(cg):
    """Expand a 6-bit cgenom value to a 5×5 binary numpy array."""
    return ((cg >> ORBIT_MAP) & 1).astype(np.uint8)


def _find_lib():
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.join(here, "..")
    for name in ("libevoca.dylib", "libevoca.so"):
        p = os.path.realpath(os.path.join(root, "C", name))
        if os.path.exists(p):
            return p
    raise FileNotFoundError(
        "libevoca not found. Run `make` (or the gcc command in CLAUDE.md) first.")


class EvoCA:
    """Thin ctypes wrapper around the EvoCA shared library."""

    def __init__(self, lib_path=None):
        self._lib        = ctypes.CDLL(lib_path or _find_lib())
        self._N          = 0
        self.food_inc    = 0.0
        self.m_scale     = 0.0
        self.food_repro  = 0.5
        self.gdiff       = 0
        self.mu_lut      = 0.0
        self.mu_cgenom   = 0.0
        self.tax         = 0.0
        self.restricted_mu = 0
        self.cgenom      = 0
        self._setup_signatures()

    def _setup_signatures(self):
        L = self._lib
        L.evoca_init.argtypes           = [ctypes.c_int, ctypes.c_float,
                                            ctypes.c_float, ctypes.c_float]
        L.evoca_init.restype            = None
        L.evoca_free.argtypes           = []
        L.evoca_free.restype            = None
        L.evoca_set_food_inc.argtypes   = [ctypes.c_float]
        L.evoca_set_food_inc.restype    = None
        L.evoca_set_m_scale.argtypes    = [ctypes.c_float]
        L.evoca_set_m_scale.restype     = None
        L.evoca_set_food_repro.argtypes = [ctypes.c_float]
        L.evoca_set_food_repro.restype  = None
        L.evoca_set_gdiff.argtypes      = [ctypes.c_int]
        L.evoca_set_gdiff.restype       = None
        L.evoca_get_gdiff.argtypes      = []
        L.evoca_get_gdiff.restype       = ctypes.c_int
        L.evoca_set_mu_lut.argtypes     = [ctypes.c_float]
        L.evoca_set_mu_lut.restype      = None
        L.evoca_set_mu_cgenom.argtypes  = [ctypes.c_float]
        L.evoca_set_mu_cgenom.restype   = None
        L.evoca_get_mu_lut.argtypes     = []
        L.evoca_get_mu_lut.restype      = ctypes.c_float
        L.evoca_get_mu_cgenom.argtypes  = []
        L.evoca_get_mu_cgenom.restype   = ctypes.c_float
        L.evoca_set_tax.argtypes        = [ctypes.c_float]
        L.evoca_set_tax.restype         = None
        L.evoca_get_tax.argtypes        = []
        L.evoca_get_tax.restype         = ctypes.c_float
        L.evoca_set_restricted_mu.argtypes = [ctypes.c_int]
        L.evoca_set_restricted_mu.restype  = None
        L.evoca_get_restricted_mu.argtypes = []
        L.evoca_get_restricted_mu.restype  = ctypes.c_int
        L.evoca_set_v_all.argtypes      = [ctypes.POINTER(ctypes.c_uint8),
                                            ctypes.c_int]
        L.evoca_set_v_all.restype       = None
        L.evoca_set_lut_all.argtypes    = [ctypes.POINTER(ctypes.c_uint8)]
        L.evoca_set_lut_all.restype     = None
        L.evoca_set_lut.argtypes        = [ctypes.c_int,
                                            ctypes.POINTER(ctypes.c_uint8)]
        L.evoca_set_lut.restype         = None
        L.evoca_set_cgenom_all.argtypes = [ctypes.c_uint8]
        L.evoca_set_cgenom_all.restype  = None
        L.evoca_set_f_all.argtypes      = [ctypes.c_float]
        L.evoca_set_f_all.restype       = None
        L.evoca_set_F_all.argtypes      = [ctypes.c_float]
        L.evoca_set_F_all.restype       = None
        L.evoca_step.argtypes           = []
        L.evoca_step.restype            = None
        L.evoca_colorize.argtypes       = [ctypes.POINTER(ctypes.c_int32),
                                            ctypes.c_int]
        L.evoca_colorize.restype        = None
        L.evoca_get_v.argtypes          = []
        L.evoca_get_v.restype           = ctypes.POINTER(ctypes.c_uint8)
        L.evoca_get_F.argtypes          = []
        L.evoca_get_F.restype           = ctypes.POINTER(ctypes.c_float)
        L.evoca_get_f.argtypes          = []
        L.evoca_get_f.restype           = ctypes.POINTER(ctypes.c_float)
        L.evoca_get_cgenom.argtypes     = []
        L.evoca_get_cgenom.restype      = ctypes.POINTER(ctypes.c_uint8)
        L.evoca_get_lut.argtypes        = []
        L.evoca_get_lut.restype         = ctypes.POINTER(ctypes.c_uint8)
        L.evoca_get_births.argtypes     = []
        L.evoca_get_births.restype      = ctypes.POINTER(ctypes.c_uint8)
        L.evoca_get_N.argtypes          = []
        L.evoca_get_N.restype           = ctypes.c_int
        L.evoca_get_cell_px.argtypes    = []
        L.evoca_get_cell_px.restype     = ctypes.c_int
        L.evoca_set_act_ymax.argtypes   = [ctypes.c_int]
        L.evoca_set_act_ymax.restype    = None
        L.evoca_get_act_ymax.argtypes   = []
        L.evoca_get_act_ymax.restype    = ctypes.c_int
        L.evoca_get_repro_age_hist.argtypes = []
        L.evoca_get_repro_age_hist.restype  = ctypes.POINTER(ctypes.c_uint32)
        L.evoca_get_repro_age_max.argtypes  = []
        L.evoca_get_repro_age_max.restype   = ctypes.c_int
        L.evoca_get_step.argtypes           = []
        L.evoca_get_step.restype            = ctypes.c_uint32
        L.evoca_set_repro_age_t0.argtypes   = [ctypes.c_uint32]
        L.evoca_set_repro_age_t0.restype    = None
        L.evoca_get_repro_age_t0.argtypes   = []
        L.evoca_get_repro_age_t0.restype    = ctypes.c_uint32
        L.evoca_reset_repro_age_hist.argtypes = []
        L.evoca_reset_repro_age_hist.restype  = None
        L.evoca_activity_update.argtypes  = []
        L.evoca_activity_update.restype   = None
        L.evoca_activity_render_col.argtypes = [
            ctypes.POINTER(ctypes.c_int32), ctypes.c_int]
        L.evoca_activity_render_col.restype  = None
        L.evoca_activity_get.argtypes   = [
            ctypes.POINTER(ctypes.c_uint32),
            ctypes.POINTER(ctypes.c_uint64),
            ctypes.POINTER(ctypes.c_uint32),
            ctypes.POINTER(ctypes.c_int32),
            ctypes.c_int]
        L.evoca_activity_get.restype    = ctypes.c_int
        # Cgenom activity
        L.evoca_cg_activity_update.argtypes  = []
        L.evoca_cg_activity_update.restype   = None
        L.evoca_cg_activity_render_col.argtypes = [
            ctypes.POINTER(ctypes.c_int32), ctypes.c_int]
        L.evoca_cg_activity_render_col.restype  = None
        L.evoca_cg_activity_get.argtypes = [
            ctypes.POINTER(ctypes.c_uint64),
            ctypes.POINTER(ctypes.c_uint32),
            ctypes.POINTER(ctypes.c_int32)]
        L.evoca_cg_activity_get.restype  = ctypes.c_int
        L.evoca_set_cg_act_ymax.argtypes = [ctypes.c_int]
        L.evoca_set_cg_act_ymax.restype  = None
        L.evoca_get_cg_act_ymax.argtypes = []
        L.evoca_get_cg_act_ymax.restype  = ctypes.c_int
        # LUT complexity
        L.evoca_lut_complexity_counts.argtypes = [
            ctypes.POINTER(ctypes.c_uint32)]
        L.evoca_lut_complexity_counts.restype  = None
        L.evoca_lut_complexity_render_col.argtypes = [
            ctypes.POINTER(ctypes.c_int32), ctypes.c_int]
        L.evoca_lut_complexity_render_col.restype  = None

    # ── Lifecycle ──────────────────────────────────────────────────────

    def init(self, N, food_inc=0.0, m_scale=1.0, food_repro=0.5, gdiff=0,
             mu_lut=0.0, mu_cgenom=0.0, tax=0.0, restricted_mu=0):
        stop = getattr(self, '_stop_display', None)
        if stop is not None:
            stop()
            self._stop_display = None
        self._N         = N
        self.food_inc   = float(food_inc)
        self.m_scale    = float(m_scale)
        self.food_repro = float(food_repro)
        self.gdiff      = int(gdiff)
        self.mu_lut     = float(mu_lut)
        self.mu_cgenom  = float(mu_cgenom)
        self.tax        = float(tax)
        self.restricted_mu = int(restricted_mu)
        self._lib.evoca_init(N, self.food_inc, self.m_scale, self.food_repro)
        self._lib.evoca_set_gdiff(self.gdiff)
        self._lib.evoca_set_mu_lut(self.mu_lut)
        self._lib.evoca_set_mu_cgenom(self.mu_cgenom)
        self._lib.evoca_set_tax(self.tax)
        self._lib.evoca_set_restricted_mu(self.restricted_mu)

    def free(self):
        stop = getattr(self, '_stop_display', None)
        if stop is not None:
            stop()
            self._stop_display = None
        self._lib.evoca_free()
        self._N = 0

    def __del__(self):
        try:
            if self._N:
                self.free()
        except Exception:
            pass

    # ── Metaparam setters ──────────────────────────────────────────────

    def update_food_inc(self, f):
        self.food_inc = float(f)
        self._lib.evoca_set_food_inc(self.food_inc)

    def update_m_scale(self, m):
        self.m_scale = float(m)
        self._lib.evoca_set_m_scale(self.m_scale)

    def update_food_repro(self, r):
        self.food_repro = float(r)
        self._lib.evoca_set_food_repro(self.food_repro)

    def update_gdiff(self, d):
        self.gdiff = int(d)
        self._lib.evoca_set_gdiff(self.gdiff)

    def update_mu_lut(self, m):
        self.mu_lut = float(m)
        self._lib.evoca_set_mu_lut(self.mu_lut)

    def update_mu_cgenom(self, m):
        self.mu_cgenom = float(m)
        self._lib.evoca_set_mu_cgenom(self.mu_cgenom)

    def update_tax(self, t):
        self.tax = float(t)
        self._lib.evoca_set_tax(self.tax)

    def update_restricted_mu(self, r):
        self.restricted_mu = int(r)
        self._lib.evoca_set_restricted_mu(self.restricted_mu)

    def update_act_ymax(self, y):
        self._lib.evoca_set_act_ymax(int(y))

    def update_cg_act_ymax(self, y):
        self._lib.evoca_set_cg_act_ymax(int(y))

    # ── Params export ─────────────────────────────────────────────────

    _DEFAULTS = dict(food_inc=0.0, m_scale=1.0, food_repro=0.5,
                      gdiff=0, mu_lut=0.0, mu_cgenom=0.0, tax=0.0,
                      restricted_mu=0)

    def params(self):
        """Return current metaparameters as a dict suitable for init(**d)."""
        return dict(N=self._N, food_inc=self.food_inc, m_scale=self.m_scale,
                    food_repro=self.food_repro, gdiff=self.gdiff,
                    mu_lut=self.mu_lut, mu_cgenom=self.mu_cgenom,
                    tax=self.tax, restricted_mu=self.restricted_mu)

    def params_str(self):
        """Return a copy-pasteable sim.init(...) call with defaults annotated."""
        p = self.params()
        lines = ["sim.init("]
        for k, v in p.items():
            val = repr(v)
            if k in self._DEFAULTS and v != self._DEFAULTS[k]:
                lines.append(f"    {k}={val},   # default: {self._DEFAULTS[k]!r}")
            else:
                lines.append(f"    {k}={val},")
        lines.append(")")
        return "\n".join(lines)

    # ── Grid setters ──────────────────────────────────────────────────

    def set_v(self, v_array):
        arr = np.ascontiguousarray(v_array.ravel(), dtype=np.uint8)
        self._lib.evoca_set_v_all(
            arr.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8)), len(arr))

    def set_lut_all(self, lut_bytes):
        arr = np.ascontiguousarray(lut_bytes, dtype=np.uint8)
        assert len(arr) == LUT_BYTES, f"need {LUT_BYTES} bytes, got {len(arr)}"
        self._lib.evoca_set_lut_all(
            arr.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8)))

    def set_lut(self, idx, lut_bytes):
        arr = np.ascontiguousarray(lut_bytes, dtype=np.uint8)
        assert len(arr) == LUT_BYTES
        self._lib.evoca_set_lut(
            int(idx), arr.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8)))

    def set_cgenom_all(self, cg):
        self.cgenom = int(cg) & 0x3F
        self._lib.evoca_set_cgenom_all(self.cgenom)

    def set_cgenom_random(self):
        """Set each cell's cgenom to a random value in [0, 63].
        No wild-type: all 64 cgenoms get distinct hash-based colors."""
        self._lib.evoca_set_cgenom_all(0xFF)   # wt=0xFF → no match → all colored
        ptr = self._lib.evoca_get_cgenom()
        arr = np.ctypeslib.as_array(ptr, shape=(self._N * self._N,))
        arr[:] = np.random.randint(0, 64, self._N * self._N, dtype=np.uint8)

    def set_lut_random(self, n_init=3):
        """Set each cell's LUT to an independent random rule.

        n_init: number of rings the rule conditions on (1, 2, or 3).
          1: depends on (v_x, n1) — 10 independent bits per LUT
          2: depends on (v_x, n1, n2) — 50 independent bits
          3: depends on (v_x, n1, n2, n3) — all 250 bits independent
        """
        N2 = self._N * self._N

        if n_init >= 3:
            bits = np.random.randint(0, 2, (N2, 2, 5, 5, 5), dtype=np.uint8)
        elif n_init == 2:
            core = np.random.randint(0, 2, (N2, 2, 5, 5, 1), dtype=np.uint8)
            bits = np.tile(core, (1, 1, 1, 1, 5))
        else:  # n_init == 1
            core = np.random.randint(0, 2, (N2, 2, 5, 1, 1), dtype=np.uint8)
            bits = np.tile(core, (1, 1, 1, 5, 5))

        flat = bits.reshape(N2, LUT_BITS)

        # Bit-pack to (N2, 32) bytes
        padded = np.zeros((N2, 256), dtype=np.uint8)
        padded[:, :LUT_BITS] = flat
        grouped = padded.reshape(N2, 32, 8)
        packed = np.zeros((N2, 32), dtype=np.uint8)
        for b in range(8):
            packed |= grouped[:, :, b] << b

        packed = np.ascontiguousarray(packed)
        for i in range(N2):
            self.set_lut(i, packed[i])

    def set_f_all(self, f):
        self._lib.evoca_set_f_all(float(f))

    def set_F_all(self, F):
        self._lib.evoca_set_F_all(float(F))

    def set_F_random(self, lo=0.0, hi=1.0):
        """Set env food F(x) to uniform random values in [lo, hi]."""
        ptr = self._lib.evoca_get_F()
        arr = np.ctypeslib.as_array(ptr, shape=(self._N * self._N,))
        arr[:] = np.random.uniform(lo, hi, self._N * self._N).astype(np.float32)

    # ── Step and colorize ─────────────────────────────────────────────

    def step(self):
        self._lib.evoca_step()

    def colorize(self, pixels, colormode=0):
        """Fill a (N*N,) int32 numpy array with ARGB values in-place."""
        arr = np.ascontiguousarray(pixels, dtype=np.int32)
        self._lib.evoca_colorize(
            arr.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
            int(colormode))
        return arr

    # ── Getters ───────────────────────────────────────────────────────

    def get_v(self):
        ptr = self._lib.evoca_get_v()
        return np.ctypeslib.as_array(ptr, shape=(self._N * self._N,)) \
                 .copy().reshape(self._N, self._N)

    def get_F(self):
        ptr = self._lib.evoca_get_F()
        return np.ctypeslib.as_array(ptr, shape=(self._N * self._N,)) \
                 .copy().reshape(self._N, self._N)

    def get_f(self):
        ptr = self._lib.evoca_get_f()
        return np.ctypeslib.as_array(ptr, shape=(self._N * self._N,)) \
                 .copy().reshape(self._N, self._N)

    def get_cgenom(self):
        """Return (N, N) uint8 array of fiducial genomes (6-bit values)."""
        ptr = self._lib.evoca_get_cgenom()
        return np.ctypeslib.as_array(ptr, shape=(self._N * self._N,)) \
                 .copy().reshape(self._N, self._N)

    def get_births(self):
        """Return (N, N) uint8 array of birth events from last step."""
        ptr = self._lib.evoca_get_births()
        return np.ctypeslib.as_array(ptr, shape=(self._N * self._N,)) \
                 .copy().reshape(self._N, self._N)

    def get_lut(self, idx):
        """Return LUT_BYTES-length uint8 array for cell at flat index idx."""
        ptr = self._lib.evoca_get_lut()
        arr = np.ctypeslib.as_array(ptr, shape=(self._N * self._N * LUT_BYTES,))
        off = int(idx) * LUT_BYTES
        return arr[off:off + LUT_BYTES].copy()

    def get_lut_complexity(self):
        """Return LUT complexity counts: {n1: count, n2: count, n3: count}."""
        counts = np.zeros(3, dtype=np.uint32)
        self._lib.evoca_lut_complexity_counts(
            counts.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)))
        return {'n1': int(counts[0]), 'n2': int(counts[1]), 'n3': int(counts[2])}

    def get_repro_age_hist(self):
        """Return reproduction age histogram as a numpy array (copy)."""
        n = self._lib.evoca_get_repro_age_max()
        ptr = self._lib.evoca_get_repro_age_hist()
        return np.ctypeslib.as_array(ptr, shape=(n,)).copy()

    def set_repro_age_t0(self, t):
        """Set the step at which reproduction age accumulation begins."""
        self._lib.evoca_set_repro_age_t0(int(t))

    def reset_repro_age_hist(self):
        """Clear the reproduction age histogram."""
        self._lib.evoca_reset_repro_age_hist()

    def get_step(self):
        """Return the current global step counter."""
        return int(self._lib.evoca_get_step())

    def get_cg_activity(self):
        """Return cgenom activity table as dict of arrays (64 entries)."""
        acts = np.zeros(64, dtype=np.uint64)
        pops = np.zeros(64, dtype=np.uint32)
        cols = np.zeros(64, dtype=np.int32)
        self._lib.evoca_cg_activity_get(
            acts.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),
            pops.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
            cols.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)))
        return {'activity': acts, 'pop_count': pops, 'color': cols}

    def get_activity(self, max_n=4096):
        """Return activity table as dict of arrays.

        Returns dict with keys: 'hash', 'activity', 'pop_count', 'color',
        each a numpy array of length n (actual number of entries).
        """
        keys = np.zeros(max_n, dtype=np.uint32)
        acts = np.zeros(max_n, dtype=np.uint64)
        pops = np.zeros(max_n, dtype=np.uint32)
        cols = np.zeros(max_n, dtype=np.int32)
        n = self._lib.evoca_activity_get(
            keys.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
            acts.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),
            pops.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
            cols.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
            max_n)
        return {'hash': keys[:n], 'activity': acts[:n],
                'pop_count': pops[:n], 'color': cols[:n]}

    @property
    def N(self):
        return self._N

    @property
    def cell_px(self):
        return self._lib.evoca_get_cell_px()


# ── LUT construction helpers ───────────────────────────────────────────

def pack_lut(bits):
    """Pack a LUT_BITS-length uint8 array of 0/1 values into LUT_BYTES bytes."""
    packed = np.zeros(LUT_BYTES, dtype=np.uint8)
    for i, b in enumerate(bits[:LUT_BITS]):
        if b:
            packed[i >> 3] |= np.uint8(1 << (i & 7))
    return packed


def unpack_lut(packed):
    """Unpack LUT_BYTES bytes into a LUT_BITS-length uint8 0/1 array."""
    bits = np.zeros(LUT_BITS, dtype=np.uint8)
    for i in range(LUT_BITS):
        bits[i] = (packed[i >> 3] >> (i & 7)) & 1
    return bits


def make_gol_lut():
    """
    Build a bit-packed LUT implementing exact Conway's Game of Life.

    GoL rule (B3/S23) on the Moore neighbourhood:
      dead  cell (v_x=0): birth  iff n1+n2 == 3
      alive cell (v_x=1): survive iff n1+n2 in {2, 3}
      n3 is ignored — GoL does not use it.

    Returns a LUT_BYTES-length uint8 array (bit-packed).
    """
    bits = np.zeros(LUT_BITS, dtype=np.uint8)
    for v_x in range(2):
        for n1 in range(5):
            for n2 in range(5):
                moore = n1 + n2
                if v_x == 0:
                    new_state = 1 if moore == 3 else 0
                else:
                    new_state = 1 if moore in (2, 3) else 0
                for n3 in range(5):
                    bits[lut_bit_index(v_x, n1, n2, n3)] = new_state
    return pack_lut(bits)
