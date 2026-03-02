"""
display.py — SDL2 display and metaparam widget for EvoCA.

Architecture (mirrors genelife):
  - C fills an int32 ARGB pixel buffer via evoca_colorize()
  - Python reshapes and copies it to the SDL2 surface via pixels2d()
  - Widget row below the grid: colored blocks, click-to-adjust params
  - Window title shows current metaparam values and fps

Usage
-----
    from python.display import run
    from python.evoca_py import EvoCA, make_gol_lut
    sim = EvoCA()
    sim.init(512, food_inc=0.01, m_scale=1.0, food_repro=0.5)
    sim.set_lut_all(make_gol_lut())
    ...
    run(sim)           # blocks until window is closed

Keyboard shortcuts
------------------
    SPACE  — pause / unpause
    C      — cycle colour mode  (state | env-food | private-food)
    S      — single step (when paused)
    Q / Esc — quit
"""

import ctypes
import time

import numpy as np
import sdl2
import sdl2.ext

# ── Widget geometry ─────────────────────────────────────────────────────
WIDGET_H  = 72    # total height of widget strip below grid (pixels)
ROW_H     = 18    # height of each param row
BTN_W     = 20    # width of the -/+ click zones at each end of a param row

# ── Metaparam definitions ────────────────────────────────────────────────
# (attr_name, step, min, max, colour_RGB)
PARAMS = [
    ("food_inc",   0.002, 0.0,  0.5,  (50,  200, 80)),   # green
    ("m_scale",    0.1,   0.0,  10.0, (80,  130, 220)),   # blue
    ("food_repro", 0.05,  0.0,  2.0,  (220, 100, 50)),    # orange
]
N_PARAMS = len(PARAMS)

# colour-mode labels (displayed in window title)
COLOR_MODES = ["state", "env-food", "priv-food"]


# ── Drawing helpers ─────────────────────────────────────────────────────

def fill_rect(surface, rgb, x, y, w, h):
    """Fill a rectangle [x,y,w,h] on surface with colour (r,g,b)."""
    sdl2.ext.fill(surface, sdl2.ext.Color(rgb[0], rgb[1], rgb[2]),
                  (x, y, w, h))


def draw_widgets(surface, sim, N, colormode, paused):
    """Render the widget strip below the N×N grid."""
    W = N
    y0 = N   # top of widget area

    # Background
    fill_rect(surface, (30, 30, 30), 0, y0, W, WIDGET_H)

    for i, (name, step, lo, hi, col) in enumerate(PARAMS):
        y = y0 + i * ROW_H
        val = getattr(sim, name)
        frac = max(0.0, min(1.0, (val - lo) / (hi - lo) if hi > lo else 0.0))

        bar_x  = BTN_W
        bar_w  = W - 2 * BTN_W
        filled = int(frac * bar_w)

        # "-" button (left)
        fill_rect(surface, (col[0]//3, col[1]//3, col[2]//3),
                  0, y, BTN_W, ROW_H - 1)
        # filled bar (value indicator)
        fill_rect(surface, col, bar_x, y, filled, ROW_H - 1)
        # empty bar
        fill_rect(surface, (col[0]//6, col[1]//6, col[2]//6),
                  bar_x + filled, y, bar_w - filled, ROW_H - 1)
        # "+" button (right)
        fill_rect(surface, (min(255,col[0]*2), min(255,col[1]*2), min(255,col[2]*2)),
                  W - BTN_W, y, BTN_W, ROW_H - 1)

    # Bottom action row: colour-mode blocks + PAUSE + indicator
    y = y0 + N_PARAMS * ROW_H
    seg = W // (len(COLOR_MODES) + 2)

    # colour-mode buttons
    cm_cols = [(200, 200, 200), (50, 220, 80), (80, 80, 220)]
    for ci, (label, cc) in enumerate(zip(COLOR_MODES, cm_cols)):
        hi_col = cc if ci == colormode else (cc[0]//4, cc[1]//4, cc[2]//4)
        fill_rect(surface, hi_col, ci * seg, y, seg - 1, ROW_H - 1)

    # PAUSE button
    pause_col = (220, 200, 50) if paused else (80, 80, 40)
    fill_rect(surface, pause_col,
              len(COLOR_MODES) * seg, y, seg - 1, ROW_H - 1)

    # QUIT button
    fill_rect(surface, (180, 40, 40),
              (len(COLOR_MODES) + 1) * seg, y, seg - 1, ROW_H - 1)


def handle_widget_click(sim, N, mx, my):
    """
    Process a mouse click in the widget area.
    Returns (colormode_delta, pause_toggle, quit_requested).
    """
    W = N
    y0 = N
    colormode_delta = 0
    pause_toggle    = False
    quit_requested  = False

    rel_y = my - y0
    if rel_y < 0:
        return colormode_delta, pause_toggle, quit_requested

    row = rel_y // ROW_H

    if row < N_PARAMS:
        name, step, lo, hi, _ = PARAMS[row]
        if mx < BTN_W:
            new_val = max(lo, getattr(sim, name) - step)
        elif mx >= W - BTN_W:
            new_val = min(hi, getattr(sim, name) + step)
        else:
            return colormode_delta, pause_toggle, quit_requested
        # update via the sim's setter method
        getattr(sim, f"update_{name}")(new_val)

    elif row == N_PARAMS:   # action row
        seg = W // (len(COLOR_MODES) + 2)
        btn = mx // seg
        if btn < len(COLOR_MODES):
            colormode_delta = btn   # select specific colormode
        elif btn == len(COLOR_MODES):
            pause_toggle = True
        else:
            quit_requested = True

    return colormode_delta, pause_toggle, quit_requested


# ── Main run loop ────────────────────────────────────────────────────────

def run(sim, colormode=0, steps=-1):
    """
    Open an SDL2 window and run the EvoCA simulation interactively.

    Parameters
    ----------
    sim        : initialised EvoCA instance
    colormode  : initial colour mode (0=state, 1=env-food, 2=priv-food)
    steps      : max steps to run (-1 = unlimited)
    """
    N = sim.N
    W = N
    H = N + WIDGET_H

    # Pixel buffer filled by C colorize
    cgolg = np.zeros(N * N, dtype=np.int32)

    sdl2.ext.init()
    window = sdl2.ext.Window(
        "EvoCA", (W, H), (500, 60),
        sdl2.SDL_WINDOW_SHOWN | sdl2.SDL_WINDOW_INPUT_FOCUS |
        sdl2.SDL_WINDOW_MOUSE_FOCUS)
    surface = sdl2.ext.Window.get_surface(window)
    sdl2.SDL_SetSurfaceBlendMode(surface, sdl2.SDL_BLENDMODE_NONE)

    # cgrid is indexed [x][y] (width × height) — matches SDL surface layout
    cgrid = sdl2.ext.pixels2d(surface)

    paused   = False
    step_cnt = 0
    t_last   = time.perf_counter()
    fps      = 0.0
    fps_alpha = 0.1   # exponential smoothing factor

    running = True
    while running:

        # ── Event handling ─────────────────────────────────────────
        for event in sdl2.ext.get_events():
            if event.type == sdl2.SDL_QUIT:
                running = False

            elif event.type == sdl2.SDL_KEYDOWN:
                k = event.key.keysym.sym
                if k in (sdl2.SDLK_q, sdl2.SDLK_ESCAPE):
                    running = False
                elif k == sdl2.SDLK_SPACE:
                    paused = not paused
                elif k == sdl2.SDLK_c:
                    colormode = (colormode + 1) % len(COLOR_MODES)
                elif k == sdl2.SDLK_s and paused:
                    # single step when paused
                    sim.step()
                    step_cnt += 1

            elif event.type == sdl2.SDL_MOUSEBUTTONDOWN:
                mx = event.button.x
                my = event.button.y
                if my >= N:   # click in widget area
                    cm_btn, pt, quit_req = handle_widget_click(
                        sim, N, mx, my)
                    if quit_req:
                        running = False
                    elif pt:
                        paused = not paused
                    elif cm_btn >= 0:
                        colormode = cm_btn % len(COLOR_MODES)

        if not running:
            break

        # ── Simulation step ────────────────────────────────────────
        if not paused:
            sim.step()
            step_cnt += 1
            if steps > 0 and step_cnt >= steps:
                running = False

        # ── Render grid ────────────────────────────────────────────
        sim.colorize(cgolg, colormode)
        # Transpose: cgolg is row-major [row*N+col]; SDL surface is [x][y]
        cgrid[:, :N] = cgolg.reshape(N, N).T

        # ── Render widgets ─────────────────────────────────────────
        draw_widgets(surface, sim, N, colormode, paused)

        # ── FPS measurement and window title ───────────────────────
        t_now = time.perf_counter()
        dt = t_now - t_last
        t_last = t_now
        if dt > 0:
            fps = fps_alpha * (1.0 / dt) + (1 - fps_alpha) * fps

        window.title = (
            f"EvoCA  t={step_cnt}"
            f"  food_inc={sim.food_inc:.4f}"
            f"  m={sim.m_scale:.3f}"
            f"  repro={sim.food_repro:.3f}"
            f"  color={COLOR_MODES[colormode]}"
            f"  {'PAUSED' if paused else f'fps={fps:.1f}'}"
        )

        sdl2.ext.Window.refresh(window)

    sdl2.ext.quit()
