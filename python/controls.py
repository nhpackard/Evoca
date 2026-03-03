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

COLOR_MODES      = ["state", "env-food", "priv-food"]
_SLIDER_RESUME_S = 0.20   # seconds after last slider touch before auto-resume
_WORKER          = os.path.join(os.path.dirname(__file__), "sdl_worker.py")

# ctrl_shm layout (5 × int32)
_QUIT, _CMODE, _STEP, _FPS10, _PAUSED = 0, 1, 2, 3, 4


def run_with_controls(sim, cell_px=None, colormode=0, paused=False):
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

    Returns
    -------
    threading.Thread — the simulation thread (can be .join()-ed if desired)
    """
    N  = sim.N
    px = cell_px if cell_px is not None else sim.cell_px

    # ── Shared memory ─────────────────────────────────────────────
    pixel_shm = SharedMemory(create=True, size=N * N * 4)
    ctrl_shm  = SharedMemory(create=True, size=5 * 4)   # 5 int32

    # numpy views into shared memory
    pixels = np.ndarray((N * N,), dtype=np.int32, buffer=pixel_shm.buf)
    ctrl   = np.ndarray((5,),     dtype=np.int32, buffer=ctrl_shm.buf)
    ctrl[:] = [0, colormode, 0, 0, int(paused)]

    # ── SDL2 subprocess ───────────────────────────────────────────
    sdl_proc = subprocess.Popen(
        [sys.executable, _WORKER,
         pixel_shm.name, ctrl_shm.name, str(N), str(px)],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
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
        for shm in (pixel_shm, ctrl_shm):
            try:
                shm.unlink()
            except Exception:
                pass
            try:
                shm.close()
            except Exception:
                pass

    atexit.register(_do_cleanup)

    # ── ipywidgets ─────────────────────────────────────────────────
    btn_pause = widgets.ToggleButton(
        value=bool(paused), description="Run" if paused else "Pause",
        button_style="", layout=widgets.Layout(width="90px"))
    btn_step  = widgets.Button(
        description="Step",  layout=widgets.Layout(width="70px"))
    btn_quit  = widgets.Button(
        description="Quit",  button_style="danger",
        layout=widgets.Layout(width="70px"))

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
        widgets.HBox([btn_pause, btn_step, btn_quit]),
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

    btn_pause.observe(on_pause_toggle, names='value')
    btn_step.on_click(on_step)
    btn_quit.on_click(on_quit)
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
