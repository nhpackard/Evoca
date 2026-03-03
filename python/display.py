"""
display.py — Minimal SDL2 display for EvoCA (no widgets).

Scales grid by sim.cell_px (from CELL_PX in evoca.h).
Drift fix: pixels2d() is acquired fresh every frame and released
(del) before refresh(), avoiding stale-view / SDL-lock conflicts.

Usage
-----
    from python.display import run
    run(sim)           # blocks until window is closed

Keyboard
--------
    SPACE     pause / unpause
    C         cycle colour mode
    S         single step (when paused)
    Q / Esc   quit
"""

import time

import numpy as np
import sdl2
import sdl2.ext

COLOR_MODES = ["state", "env-food", "priv-food", "births"]


def run(sim, cell_px=None, colormode=0, steps=-1):
    """
    Open a plain SDL2 window and run the simulation (blocking).

    Parameters
    ----------
    sim       : initialised EvoCA instance
    cell_px   : screen pixels per cell (default: sim.cell_px from C #define)
    colormode : initial colour mode (0=state, 1=env-food, 2=priv-food)
    steps     : max steps to run (-1 = unlimited)
    """
    N  = sim.N
    px = cell_px if cell_px is not None else sim.cell_px
    W  = N * px
    H  = N * px

    cgolg = np.zeros(N * N, dtype=np.int32)

    sdl2.ext.init()
    window  = sdl2.ext.Window("EvoCA", (W, H), (500, 60))
    surface = sdl2.ext.Window.get_surface(window)
    sdl2.SDL_SetSurfaceBlendMode(surface, sdl2.SDL_BLENDMODE_NONE)

    paused    = False
    step_cnt  = 0
    t_last    = time.perf_counter()
    fps       = 0.0
    fps_alpha = 0.1

    running = True
    while running:

        # ── Events ─────────────────────────────────────────────────
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
                    sim.step()
                    step_cnt += 1

        if not running:
            break

        # ── Simulation step ────────────────────────────────────────
        if not paused:
            sim.step()
            step_cnt += 1
            if steps > 0 and step_cnt >= steps:
                running = False

        # ── Render — fresh pixels2d each frame to prevent drift ────
        sim.colorize(cgolg, colormode)
        cgrid = sdl2.ext.pixels2d(surface)
        if px == 1:
            cgrid[:, :] = cgolg.reshape(N, N).T
        else:
            cgrid[:, :] = np.repeat(
                np.repeat(cgolg.reshape(N, N), px, axis=0),
                px, axis=1).T
        del cgrid   # release SDL surface lock before refresh

        # ── FPS and title ──────────────────────────────────────────
        t_now = time.perf_counter()
        dt = t_now - t_last
        t_last = t_now
        if dt > 0:
            fps = fps_alpha * (1.0 / dt) + (1.0 - fps_alpha) * fps

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
