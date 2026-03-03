"""
controls.py — ipywidgets control panel + SDL2 display for EvoCA.

Architecture
------------
SDL2 requires macOS's actual main thread (thread 0), which belongs to the
Jupyter kernel's event loop.  Running SDL2 in any other thread crashes the
kernel.  The fix: run SDL2 in a *subprocess* that has its own main thread.

  Main process (Jupyter kernel)
  ├── Jupyter event loop  →  ipywidgets callbacks fire here
  ├── Sim thread          →  sim.step() + sim.colorize() → shared memory
  ├── Reader thread       →  relays SDL worker stdout to terminal
  └── subprocess (sdl_worker.py)
      └── SDL2 main thread  →  reads shared memory, renders window

Pixel data flows via POSIX shared memory with no copying:
  sim.colorize(pixels_shm, …) writes directly into shared memory;
  the subprocess reads the same buffer.

Usage (Jupyter cell)
--------------------
    from python.controls import run_with_controls
    run_with_controls(sim_sdl)
    # Cell returns immediately; widgets appear below; SDL2 window opens.
    # Click Quit or press Q/Esc in the SDL2 window to stop.
    # Then call sim_sdl.free() in the next cell.
"""

import atexit
import os
import sys
import subprocess
import threading
import time

import numpy as np
import ipywidgets as widgets
import matplotlib.pyplot as plt
from IPython.display import display as ipy_display
from .evoca_py import cgenom_to_pattern
from multiprocessing.shared_memory import SharedMemory

import ctypes

COLOR_MODES      = ["state", "env-food", "priv-food", "births"]
_SLIDER_RESUME_S = 0.20   # seconds after last slider touch before auto-resume
_WORKER          = os.path.join(os.path.dirname(__file__), "sdl_worker.py")

# ctrl_shm layout (5 × int32)
_QUIT, _CMODE, _STEP, _FPS10, _PAUSED = 0, 1, 2, 3, 4

# Probe strip-chart constants
PROBE_W = 512    # pixels = time steps visible
PROBE_H = 128    # pixel height of each strip chart

# Map probe name → (C getter function name, ctype element type)
_PROBE_GETTER = {
    'env_food':  ('evoca_get_F',      ctypes.c_float),
    'priv_food': ('evoca_get_f',      ctypes.c_float),
    'births':    ('evoca_get_births', ctypes.c_uint8),
}

# Module-level handle: stop any previous session before starting a new one.
_active_stop = None


def run_with_controls(sim, cell_px=None, colormode=0, paused=False, probes=None):
    """
    Display ipywidgets controls and open an SDL2 simulation window.

    Returns immediately (non-blocking).  The simulation runs in a background
    thread; the SDL2 display runs in a subprocess.

    Parameters
    ----------
    sim       : initialised EvoCA instance
    cell_px   : screen pixels per cell (default: sim.cell_px from CELL_PX #define)
    colormode : initial colour mode (0=state, 1=env-food, 2=priv-food)
    paused    : if True, start in paused state
    probes    : dict of probe names to enable, e.g. {'env_food': True, 'priv_food': True}

    Returns
    -------
    threading.Thread — the simulation thread (can be .join()-ed if desired)
    """
    # ── Tear down any previous session ────────────────────────────
    global _active_stop
    if _active_stop is not None:
        _active_stop()
        _active_stop = None

    N  = sim.N
    px = cell_px if cell_px is not None else sim.cell_px

    # ── Shared memory ─────────────────────────────────────────────
    pixel_shm = SharedMemory(create=True, size=N * N * 4)
    ctrl_shm  = SharedMemory(create=True, size=5 * 4)   # 5 int32

    # numpy views into shared memory
    pixels = np.ndarray((N * N,), dtype=np.int32, buffer=pixel_shm.buf)
    ctrl   = np.ndarray((5,),     dtype=np.int32, buffer=ctrl_shm.buf)
    ctrl[:] = [0, colormode, 0, 0, int(paused)]

    # ── Probe setup ─────────────────────────────────────────────────
    probe_names = [k for k, v in (probes or {}).items() if v and k in _PROBE_GETTER]
    n_probes    = len(probe_names)
    probe_shm   = None
    probe_cursor = None       # int32 view into probe_shm
    probe_means  = []         # list of float32[PROBE_W] views
    probe_stds   = []         # list of float32[PROBE_W] views

    # Build list of (C getter func, numpy dtype) for fast access in sim thread
    probe_getters = []   # list of (getter_fn, np.dtype)
    for pname in probe_names:
        fn_name, ctype = _PROBE_GETTER[pname]
        getter = getattr(sim._lib, fn_name)
        getter.argtypes = []
        getter.restype  = ctypes.POINTER(ctype)
        probe_getters.append((getter, np.dtype(ctype)))

    if n_probes > 0:
        probe_shm_size = 4 + n_probes * 2 * PROBE_W * 4
        probe_shm = SharedMemory(create=True, size=probe_shm_size)
        _pbuf = np.ndarray((probe_shm_size,), dtype=np.uint8, buffer=probe_shm.buf)
        _pbuf[:] = 0
        probe_cursor = np.ndarray((1,), dtype=np.int32, buffer=probe_shm.buf)
        off = 4
        for _ in range(n_probes):
            probe_means.append(np.ndarray((PROBE_W,), dtype=np.float32,
                                          buffer=probe_shm.buf, offset=off))
            off += PROBE_W * 4
            probe_stds.append(np.ndarray((PROBE_W,), dtype=np.float32,
                                         buffer=probe_shm.buf, offset=off))
            off += PROBE_W * 4

    # ── SDL2 subprocess ───────────────────────────────────────────
    cmd = [sys.executable, _WORKER,
           pixel_shm.name, ctrl_shm.name, str(N), str(px)]
    if n_probes > 0:
        cmd += [probe_shm.name, ",".join(probe_names)]
    sdl_proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    )

    # Relay subprocess output to the terminal (for diagnostics)
    def _reader():
        for line in sdl_proc.stdout:
            print(line.decode(errors='replace'), end='', flush=True)

    threading.Thread(target=_reader, name="evoca-sdl-reader", daemon=True).start()

    # ── Shared sim state ──────────────────────────────────────────
    st    = dict(paused=bool(paused), running=True, colormode=colormode, step_cnt=0)
    # _alive guards widget callbacks against use-after-free:
    # set to False just before shared memory is released.
    _alive        = [True]
    _cleanup_done = [False]

    # ── Guaranteed cleanup (sim thread or atexit, whichever fires first) ──
    def _do_cleanup():
        if _cleanup_done[0]:
            return
        _cleanup_done[0] = True
        _alive[0] = False
        if sdl_proc.poll() is None:
            sdl_proc.terminate()
            try:
                sdl_proc.wait(timeout=2)
            except Exception:
                pass
        all_shm = [pixel_shm, ctrl_shm]
        if probe_shm is not None:
            all_shm.append(probe_shm)
        for shm in all_shm:
            try:
                shm.unlink()
            except Exception:
                pass
            try:
                shm.close()
            except Exception:
                pass

    atexit.register(_do_cleanup)

    def _stop():
        """Signal sim thread to quit and wait for cleanup."""
        global _active_stop
        st['running'] = False
        if _alive[0]:
            ctrl[_QUIT] = 1
        # Give sim thread time to finish and call _do_cleanup
        time.sleep(0.1)
        _do_cleanup()
        _active_stop = None

    _active_stop = _stop
    sim._stop_display = _stop

    # ── ipywidgets ─────────────────────────────────────────────────
    btn_pause = widgets.ToggleButton(
        value=bool(paused), description="Run" if paused else "Pause",
        button_style="", layout=widgets.Layout(width="90px"))
    btn_step  = widgets.Button(
        description="Step",  layout=widgets.Layout(width="70px"))
    btn_quit  = widgets.Button(
        description="Quit",  button_style="danger",
        layout=widgets.Layout(width="70px"))
    btn_save  = widgets.Button(
        description="Save Plots", layout=widgets.Layout(width="100px"))

    sl_kw = dict(continuous_update=True,
                 style={"description_width": "90px"},
                 layout=widgets.Layout(width="440px"))
    sl_food_inc = widgets.FloatSlider(
        value=sim.food_inc,   min=0.0, max=0.5,  step=0.001,
        description="food_inc:",   readout_format=".3f", **sl_kw)
    sl_m_scale = widgets.FloatSlider(
        value=sim.m_scale,    min=0.0, max=10.0, step=0.1,
        description="m_scale:",    readout_format=".2f", **sl_kw)
    sl_food_repro = widgets.FloatSlider(
        value=sim.food_repro, min=0.0, max=2.0,  step=0.05,
        description="food_repro:", readout_format=".2f", **sl_kw)
    sl_gdiff = widgets.IntSlider(
        value=sim.gdiff, min=0, max=10, step=1,
        description="gdiff:",
        style={"description_width": "90px"},
        layout=widgets.Layout(width="440px"))

    color_dd   = widgets.Dropdown(
        options=COLOR_MODES, value=COLOR_MODES[colormode],
        description="Color:",
        layout=widgets.Layout(width="200px"))
    status_lbl = widgets.Label(value="Starting…")

    # ── Fiducial pattern display ─────────────────────────────────────
    pat = cgenom_to_pattern(sim.cgenom)
    fig, ax = plt.subplots(figsize=(2, 2))
    ax.imshow(pat, cmap='Greys', vmin=0, vmax=1, interpolation='nearest')
    for edge in range(6):
        ax.axhline(edge - 0.5, color='gray', linewidth=0.5)
        ax.axvline(edge - 0.5, color='gray', linewidth=0.5)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f'fiducial c(x)   cg=0b{sim.cgenom:06b}', fontsize=9)
    plt.tight_layout()
    plt.show()

    ipy_display(widgets.VBox([
        widgets.HBox([btn_pause, btn_step, btn_quit, btn_save]),
        sl_food_inc, sl_m_scale, sl_food_repro, sl_gdiff,
        widgets.HBox([color_dd, status_lbl]),
    ]))

    # ── Widget callbacks ──────────────────────────────────────────
    _guard = [False]

    def _set_paused(p):
        if _guard[0] or not _alive[0]:
            return
        _guard[0] = True
        st['paused']          = p
        ctrl[_PAUSED]         = int(p)
        btn_pause.value       = p
        btn_pause.description = "Run" if p else "Pause"
        _guard[0] = False

    def on_pause_toggle(change):
        _set_paused(change['new'])

    def on_step(_):
        if not _alive[0]:
            return
        if st['paused']:
            sim.step()
            st['step_cnt'] += 1
            ctrl[_STEP] = st['step_cnt']
            sim.colorize(pixels, st['colormode'])
            if n_probes > 0:
                cur = int(probe_cursor[0])
                for pi, (getter, dt) in enumerate(probe_getters):
                    ptr = getter()
                    arr = np.ctypeslib.as_array(ptr, shape=(N * N,))
                    farr = arr.astype(np.float64) if dt != np.float32 else arr
                    probe_means[pi][cur] = farr.mean()
                    probe_stds[pi][cur]  = farr.std()
                probe_cursor[0] = (cur + 1) % PROBE_W
            status_lbl.value = f"t={st['step_cnt']}  (paused)"

    def on_quit(_):
        st['running'] = False
        if _alive[0]:
            ctrl[_QUIT] = 1
        status_lbl.value = "Stopped — call sim.free() when ready"

    def on_color(change):
        if not _alive[0]:
            return
        cm = COLOR_MODES.index(change['new'])
        st['colormode'] = cm
        ctrl[_CMODE]    = cm

    def on_save(_):
        if not _alive[0] or n_probes == 0:
            return
        cur = int(probe_cursor[0])
        for i, pname in enumerate(probe_names):
            # Reconstruct time-ordered arrays from ring buffer
            m = np.roll(probe_means[i], -cur)
            s = np.roll(probe_stds[i], -cur)
            t = np.arange(PROBE_W)
            fig_s, ax_s = plt.subplots(figsize=(8, 2.5))
            ax_s.fill_between(t, m - s, m + s, alpha=0.3)
            ax_s.plot(t, m, linewidth=0.8)
            ax_s.set_title(pname)
            ax_s.set_xlabel("t (relative)")
            fig_s.tight_layout()
            fname = f"probe_{pname}.png"
            fig_s.savefig(fname, dpi=150)
            plt.close(fig_s)
        status_lbl.value = f"Saved {n_probes} probe plot(s)"

    btn_pause.observe(on_pause_toggle, names='value')
    btn_step.on_click(on_step)
    btn_quit.on_click(on_quit)
    btn_save.on_click(on_save)
    color_dd.observe(on_color, names='value')

    # Slider drag: pause on first touch, auto-resume after 200 ms idle
    def _make_slider_cb(attr, slider):
        _timer      = [None]
        _was_paused = [False]

        def on_value(change):
            if not _alive[0]:
                return
            if not st['paused']:
                _was_paused[0] = False
                _set_paused(True)
            else:
                _was_paused[0] = True
            getattr(sim, f"update_{attr}")(change['new'])
            if _timer[0] is not None:
                _timer[0].cancel()

            def _resume():
                if not _alive[0]:
                    return
                if not _was_paused[0]:
                    st['paused']  = False
                    ctrl[_PAUSED] = 0
                    # btn_pause visual won't sync from this thread — acceptable

            _timer[0] = threading.Timer(_SLIDER_RESUME_S, _resume)
            _timer[0].start()

        slider.observe(on_value, names='value')

    _make_slider_cb("food_inc",   sl_food_inc)
    _make_slider_cb("m_scale",    sl_m_scale)
    _make_slider_cb("food_repro", sl_food_repro)
    _make_slider_cb("gdiff",      sl_gdiff)

    # ── Simulation thread ─────────────────────────────────────────
    def _sim_thread():
        FPS_ALPHA = 0.05
        fps       = 0.0
        t_last    = time.perf_counter()

        while st['running'] and ctrl[_QUIT] == 0:

            # Check if subprocess died unexpectedly (once, after a few steps)
            if st['step_cnt'] == 5:
                rc = sdl_proc.poll()
                if rc is not None:
                    print(f"EvoCA: SDL worker exited early (rc={rc})", flush=True)
                    status_lbl.value = f"SDL worker crashed (rc={rc}) — check terminal"

            if st['paused']:
                # Colorize so Step button results are shown; then rest.
                sim.colorize(pixels, st['colormode'])
                t_last = time.perf_counter()   # reset timer; avoid fps spike on resume
                time.sleep(0.01)
                continue

            sim.step()
            st['step_cnt'] += 1

            # colorize writes directly into shared pixel memory
            sim.colorize(pixels, st['colormode'])

            # Record probe data
            if n_probes > 0:
                cur = int(probe_cursor[0])
                for pi, (getter, dt) in enumerate(probe_getters):
                    ptr = getter()
                    arr = np.ctypeslib.as_array(ptr, shape=(N * N,))
                    farr = arr.astype(np.float64) if dt != np.float32 else arr
                    probe_means[pi][cur] = farr.mean()
                    probe_stds[pi][cur]  = farr.std()
                probe_cursor[0] = (cur + 1) % PROBE_W

            # FPS (only meaningful when running)
            t_now = time.perf_counter()
            dt    = t_now - t_last
            t_last = t_now
            if dt > 0:
                fps = FPS_ALPHA * (1.0 / dt) + (1.0 - FPS_ALPHA) * fps

            # Update ctrl_shm so subprocess can update its window title
            sc           = st['step_cnt']
            ctrl[_STEP]  = sc
            ctrl[_FPS10] = int(fps * 10)

            # Update status label (infrequent — tolerate ZMQ non-safety)
            if sc % 20 == 0:
                status_lbl.value = f"t={sc}  fps={fps:.1f}"

        # Cleanup — also called by atexit if this thread doesn't finish in time
        st['running'] = False
        if _alive[0]:
            ctrl[_QUIT] = 1
        sdl_proc.wait(timeout=3)
        _do_cleanup()

    t = threading.Thread(target=_sim_thread, name="evoca-sim", daemon=True)
    t.start()
    return t
