"""
sdl_worker.py — SDL2 display subprocess for EvoCA.

Launched by controls.run_with_controls() via subprocess.Popen.
All output goes to the terminal where Jupyter was started.

Usage (internal):
    python sdl_worker.py <pixel_shm> <ctrl_shm> <N> <px> [<probe_shm> <probe_names_csv>]

ctrl_shm layout (5 × int32)
    [0] quit      1 = exit
    [1] colormode 0/1/2
    [2] step_cnt
    [3] fps × 10
    [4] paused    0/1

probe_shm layout (if present)
    [0]      int32   cursor (write position, 0..PROBE_W-1)
    [4..]    per probe: float32[PROBE_W] mean, float32[PROBE_W] std
"""
import sys
import ctypes
import traceback

PROBE_W = 512
PROBE_H = 128
TITLE_BAR_H = 28   # macOS title bar estimate

def _c(argb):
    """Convert ARGB uint32 to np.int32 (avoids numpy deprecation warning)."""
    import numpy as np
    return np.int32(argb if argb < 0x80000000 else argb - 0x100000000)

# Probe color palette (ARGB) — converted to int32 once at import
BG_COLOR      = _c(0xFF1A1A1A)
CURSOR_COLOR  = _c(0xFF444444)
_PROBE_COLORS = {
    'env_food':  (_c(0xFF00CC00), _c(0xFF003300)),   # green
    'priv_food': (_c(0xFF4488FF), _c(0xFF112244)),   # blue
    'births':    (_c(0xFFFFFF00), _c(0xFF444400)),   # yellow
}
_DEFAULT_COLORS = (_c(0xFFCCCCCC), _c(0xFF333333))   # grey fallback


def _render_probe(dst, pitch_i32, mean_buf, std_buf, cursor, color_mean, color_band,
                   y_fixed=None):
    """Render one strip chart into an SDL surface pixel array.

    y_fixed : optional (lo, hi) tuple to fix the Y axis range.
              If None, auto-scale from data.
    """
    import numpy as np

    # Clear to background
    dst[:PROBE_H, :PROBE_W] = BG_COLOR

    # Build time-ordered view: oldest on left, newest on right
    m = np.roll(mean_buf, -cursor)
    s = np.roll(std_buf, -cursor)

    lo_arr = m - s
    hi_arr = m + s
    # Only scale from non-zero data (buffer may be partially filled)
    filled = (m != 0) | (s != 0)
    if not filled.any():
        return

    if y_fixed is not None:
        y_min, y_max = y_fixed
    else:
        y_min = float(lo_arr[filled].min())
        y_max = float(hi_arr[filled].max())
    if y_max - y_min < 1e-8:
        y_min -= 0.01
        y_max += 0.01
    margin = (y_max - y_min) * 0.05
    y_min -= margin
    y_max += margin
    scale = (PROBE_H - 1) / (y_max - y_min)

    for x in range(PROBE_W):
        mv = m[x]
        sv = s[x]
        if mv == 0.0 and sv == 0.0 and not filled[x]:
            continue
        # Y pixel (0 = top = y_max, PROBE_H-1 = bottom = y_min)
        y_mean = int((y_max - mv) * scale)
        y_lo   = int((y_max - (mv + sv)) * scale)
        y_hi   = int((y_max - (mv - sv)) * scale)
        y_lo   = max(0, min(PROBE_H - 1, y_lo))
        y_hi   = max(0, min(PROBE_H - 1, y_hi))
        y_mean = max(0, min(PROBE_H - 1, y_mean))
        # Draw std band
        for y in range(y_lo, y_hi + 1):
            dst[y, x] = color_band
        # Draw mean line (overwrite band)
        dst[y_mean, x] = color_mean

    # Draw cursor position (right edge of newest data) as dim vertical line
    cx = PROBE_W - 1
    for y in range(PROBE_H):
        dst[y, cx] = CURSOR_COLOR


def main():
    if len(sys.argv) < 5:
        print("EvoCA SDL: bad args", flush=True)
        sys.exit(1)

    pixel_shm_name = sys.argv[1]
    ctrl_shm_name  = sys.argv[2]
    N  = int(sys.argv[3])
    px = int(sys.argv[4])
    W, H = N * px, N * px

    # Optional probe args (positional: argv[5]=probe_shm, argv[6]=names_csv)
    probe_shm_name = sys.argv[5] if len(sys.argv) > 6 else None
    probe_names    = sys.argv[6].split(",") if len(sys.argv) > 6 else []
    n_probes       = len(probe_names)

    # Optional activity probe (--activity=<shm_name>)
    activity_shm_name = None
    cg_activity_shm_name = None
    lut_complexity_shm_name = None
    for arg in sys.argv:
        if arg.startswith("--activity="):
            activity_shm_name = arg[len("--activity="):]
        elif arg.startswith("--cg-activity="):
            cg_activity_shm_name = arg[len("--cg-activity="):]
        elif arg.startswith("--lut-complexity="):
            lut_complexity_shm_name = arg[len("--lut-complexity="):]

    ACT_H = 2 * PROBE_H  # 256

    print(f"EvoCA SDL: starting  N={N} px={px}  probes={probe_names}", flush=True)

    import numpy as np
    from multiprocessing.shared_memory import SharedMemory
    import sdl2

    # Open shared memory
    try:
        pixel_shm = SharedMemory(name=pixel_shm_name)
        ctrl_shm  = SharedMemory(name=ctrl_shm_name)
    except Exception as e:
        print(f"EvoCA SDL: SharedMemory open failed: {e}", flush=True)
        sys.exit(1)

    pixels = np.ndarray((N * N,), dtype=np.int32, buffer=pixel_shm.buf)
    ctrl   = np.ndarray((5,),     dtype=np.int32, buffer=ctrl_shm.buf)

    # Open probe shared memory
    probe_shm    = None
    probe_cursor = None
    probe_means  = []
    probe_stds   = []
    if n_probes > 0 and probe_shm_name:
        try:
            probe_shm = SharedMemory(name=probe_shm_name)
        except Exception as e:
            print(f"EvoCA SDL: probe SharedMemory open failed: {e}", flush=True)
            n_probes = 0
            probe_names = []
        if probe_shm is not None:
            probe_cursor = np.ndarray((1,), dtype=np.int32, buffer=probe_shm.buf)
            off = 4
            for _ in range(n_probes):
                probe_means.append(np.ndarray((PROBE_W,), dtype=np.float32,
                                              buffer=probe_shm.buf, offset=off))
                off += PROBE_W * 4
                probe_stds.append(np.ndarray((PROBE_W,), dtype=np.float32,
                                             buffer=probe_shm.buf, offset=off))
                off += PROBE_W * 4

    # Open activity shared memory
    activity_shm     = None
    activity_cursor  = None
    activity_pixels  = None
    if activity_shm_name:
        try:
            activity_shm = SharedMemory(name=activity_shm_name)
            activity_cursor = np.ndarray((1,), dtype=np.int32,
                                         buffer=activity_shm.buf)
            activity_pixels = np.ndarray((ACT_H, PROBE_W), dtype=np.int32,
                                         buffer=activity_shm.buf, offset=4)
            print(f"EvoCA SDL: activity shm opened ({ACT_H}x{PROBE_W})",
                  flush=True)
        except Exception as e:
            print(f"EvoCA SDL: activity SharedMemory open failed: {e}",
                  flush=True)
            activity_shm_name = None

    # Open cgenom activity shared memory
    cg_activity_shm     = None
    cg_activity_cursor  = None
    cg_activity_pixels  = None
    if cg_activity_shm_name:
        try:
            cg_activity_shm = SharedMemory(name=cg_activity_shm_name)
            cg_activity_cursor = np.ndarray((1,), dtype=np.int32,
                                            buffer=cg_activity_shm.buf)
            cg_activity_pixels = np.ndarray((ACT_H, PROBE_W), dtype=np.int32,
                                            buffer=cg_activity_shm.buf, offset=4)
            print(f"EvoCA SDL: cg_activity shm opened ({ACT_H}x{PROBE_W})",
                  flush=True)
        except Exception as e:
            print(f"EvoCA SDL: cg_activity SharedMemory open failed: {e}",
                  flush=True)
            cg_activity_shm_name = None

    # Open LUT complexity shared memory
    lut_complexity_shm     = None
    lut_complexity_cursor  = None
    lut_complexity_pixels  = None
    if lut_complexity_shm_name:
        try:
            lut_complexity_shm = SharedMemory(name=lut_complexity_shm_name)
            lut_complexity_cursor = np.ndarray((1,), dtype=np.int32,
                                               buffer=lut_complexity_shm.buf)
            lut_complexity_pixels = np.ndarray((PROBE_H, PROBE_W), dtype=np.int32,
                                               buffer=lut_complexity_shm.buf, offset=4)
            print(f"EvoCA SDL: lut_complexity shm opened ({PROBE_H}x{PROBE_W})",
                  flush=True)
        except Exception as e:
            print(f"EvoCA SDL: lut_complexity SharedMemory open failed: {e}",
                  flush=True)
            lut_complexity_shm_name = None

    COLOR_MODES = ["state", "env-food", "priv-food", "births"]

    # ── SDL2 init ─────────────────────────────────────────────────
    if sdl2.SDL_Init(sdl2.SDL_INIT_VIDEO) != 0:
        print(f"EvoCA SDL: SDL_Init failed: {sdl2.SDL_GetError()}", flush=True)
        ctrl[0] = 1
        sys.exit(1)

    print("EvoCA SDL: SDL_Init OK", flush=True)

    # ── Screen geometry for window placement ──────────────────────
    dm = sdl2.SDL_DisplayMode()
    sdl2.SDL_GetCurrentDisplayMode(0, ctypes.byref(dm))
    scr_w, scr_h = dm.w, dm.h

    # Main window: upper-right corner
    main_x = scr_w - W
    main_y = 0

    window_p = sdl2.SDL_CreateWindow(
        b"EvoCA",
        main_x, main_y,
        W, H,
        sdl2.SDL_WINDOW_SHOWN,
    )
    if not window_p:
        print(f"EvoCA SDL: SDL_CreateWindow failed: {sdl2.SDL_GetError()}", flush=True)
        ctrl[0] = 1
        sdl2.SDL_Quit()
        sys.exit(1)

    sdl2.SDL_RaiseWindow(window_p)
    print("EvoCA SDL: window created and raised", flush=True)

    surface_p = sdl2.SDL_GetWindowSurface(window_p)
    if not surface_p:
        print(f"EvoCA SDL: SDL_GetWindowSurface failed: {sdl2.SDL_GetError()}", flush=True)
        ctrl[0] = 1
        sdl2.SDL_DestroyWindow(window_p)
        sdl2.SDL_Quit()
        sys.exit(1)

    sdl2.SDL_SetSurfaceBlendMode(surface_p, sdl2.SDL_BLENDMODE_NONE)

    surf       = surface_p.contents
    pitch_i32  = surf.pitch // 4
    pixels_ptr = ctypes.cast(surf.pixels, ctypes.POINTER(ctypes.c_int32))
    dst_flat   = np.ctypeslib.as_array(pixels_ptr, shape=(H * pitch_i32,))
    dst        = dst_flat.reshape(H, pitch_i32)

    # ── Probe windows ────────────────────────────────────────────
    # Fixed Y ranges for known probes (both food fields are in [0, 1])
    _Y_FIXED = {
        'env_food':  (0.0, 1.0),
        'priv_food': (0.0, 1.0),
        'births':    (0.0, 1.0),
    }

    probe_windows  = []   # SDL_Window pointers
    probe_surfaces = []   # SDL_Surface pointers
    probe_dsts     = []   # numpy pixel arrays
    probe_colors   = []   # (color_mean, color_band) tuples
    probe_yfix     = []   # fixed Y ranges (or None for auto-scale)

    # Stack probe windows top-down.  After creating each window, read back
    # its actual Y position so the next window is placed correctly below it
    # (macOS may adjust positions for the menu bar and window decorations).
    next_probe_y = 0   # requested Y for the next probe window
    real_title_h = TITLE_BAR_H  # updated from first window's actual position

    for i, pname in enumerate(probe_names):
        pw_x = main_x - PROBE_W
        pw = sdl2.SDL_CreateWindow(
            pname.encode(),
            pw_x, next_probe_y,
            PROBE_W, PROBE_H,
            sdl2.SDL_WINDOW_SHOWN,
        )
        if not pw:
            print(f"EvoCA SDL: probe window '{pname}' failed", flush=True)
            continue

        # Read back actual content-area position
        actual_y = ctypes.c_int(0)
        sdl2.SDL_GetWindowPosition(pw, None, ctypes.byref(actual_y))
        # Measure the window's title-bar decoration height (top border)
        if i == 0:
            top_border = ctypes.c_int(0)
            ret = sdl2.SDL_GetWindowBordersSize(
                pw, ctypes.byref(top_border), None, None, None)
            if ret == 0 and top_border.value > 0:
                real_title_h = top_border.value
            print(f"EvoCA SDL: title bar = {real_title_h}px", flush=True)
        print(f"EvoCA SDL: probe '{pname}' y={actual_y.value}", flush=True)
        # Next window: below this content + room for next window's title bar
        next_probe_y = actual_y.value + PROBE_H + real_title_h

        ps = sdl2.SDL_GetWindowSurface(pw)
        if not ps:
            sdl2.SDL_DestroyWindow(pw)
            continue
        sdl2.SDL_SetSurfaceBlendMode(ps, sdl2.SDL_BLENDMODE_NONE)
        psurf     = ps.contents
        pp_i32    = psurf.pitch // 4
        pp_ptr    = ctypes.cast(psurf.pixels, ctypes.POINTER(ctypes.c_int32))
        pd_flat   = np.ctypeslib.as_array(pp_ptr, shape=(PROBE_H * pp_i32,))
        pd        = pd_flat.reshape(PROBE_H, pp_i32)

        probe_windows.append(pw)
        probe_surfaces.append(ps)
        probe_dsts.append(pd)
        probe_colors.append(_PROBE_COLORS.get(pname, _DEFAULT_COLORS))
        probe_yfix.append(_Y_FIXED.get(pname))

    print(f"EvoCA SDL: {len(probe_windows)} probe window(s) created", flush=True)

    # ── Activity window ──────────────────────────────────────────
    act_window_p  = None
    act_surface_p = None
    act_dst       = None
    if activity_shm is not None:
        aw_x = main_x - PROBE_W
        aw = sdl2.SDL_CreateWindow(
            b"activity",
            aw_x, next_probe_y,
            PROBE_W, ACT_H,
            sdl2.SDL_WINDOW_SHOWN,
        )
        if aw:
            actual_y = ctypes.c_int(0)
            sdl2.SDL_GetWindowPosition(aw, None, ctypes.byref(actual_y))
            next_probe_y = actual_y.value + ACT_H + real_title_h
            aps = sdl2.SDL_GetWindowSurface(aw)
            if aps:
                sdl2.SDL_SetSurfaceBlendMode(aps, sdl2.SDL_BLENDMODE_NONE)
                asurf   = aps.contents
                ap_i32  = asurf.pitch // 4
                ap_ptr  = ctypes.cast(asurf.pixels,
                                      ctypes.POINTER(ctypes.c_int32))
                ad_flat = np.ctypeslib.as_array(ap_ptr,
                                                shape=(ACT_H * ap_i32,))
                act_dst       = ad_flat.reshape(ACT_H, ap_i32)
                act_window_p  = aw
                act_surface_p = aps
                print("EvoCA SDL: activity window created", flush=True)
            else:
                sdl2.SDL_DestroyWindow(aw)
        else:
            print("EvoCA SDL: activity window creation failed", flush=True)

    # ── Cgenom activity window ──────────────────────────────────
    cg_act_window_p  = None
    cg_act_surface_p = None
    cg_act_dst       = None
    if cg_activity_shm is not None:
        caw_x = main_x - PROBE_W
        caw = sdl2.SDL_CreateWindow(
            b"cg_activity",
            caw_x, next_probe_y,
            PROBE_W, ACT_H,
            sdl2.SDL_WINDOW_SHOWN,
        )
        if caw:
            actual_y = ctypes.c_int(0)
            sdl2.SDL_GetWindowPosition(caw, None, ctypes.byref(actual_y))
            next_probe_y = actual_y.value + ACT_H + real_title_h
            caps = sdl2.SDL_GetWindowSurface(caw)
            if caps:
                sdl2.SDL_SetSurfaceBlendMode(caps, sdl2.SDL_BLENDMODE_NONE)
                casurf   = caps.contents
                cap_i32  = casurf.pitch // 4
                cap_ptr  = ctypes.cast(casurf.pixels,
                                       ctypes.POINTER(ctypes.c_int32))
                cad_flat = np.ctypeslib.as_array(cap_ptr,
                                                 shape=(ACT_H * cap_i32,))
                cg_act_dst       = cad_flat.reshape(ACT_H, cap_i32)
                cg_act_window_p  = caw
                cg_act_surface_p = caps
                print("EvoCA SDL: cg_activity window created", flush=True)
            else:
                sdl2.SDL_DestroyWindow(caw)
        else:
            print("EvoCA SDL: cg_activity window creation failed", flush=True)

    # ── LUT complexity window ──────────────────────────────────
    lc_window_p  = None
    lc_surface_p = None
    lc_dst       = None
    if lut_complexity_shm is not None:
        lcw_x = main_x - PROBE_W
        lcw = sdl2.SDL_CreateWindow(
            b"lut_complexity",
            lcw_x, next_probe_y,
            PROBE_W, PROBE_H,
            sdl2.SDL_WINDOW_SHOWN,
        )
        if lcw:
            actual_y = ctypes.c_int(0)
            sdl2.SDL_GetWindowPosition(lcw, None, ctypes.byref(actual_y))
            next_probe_y = actual_y.value + PROBE_H + real_title_h
            lps = sdl2.SDL_GetWindowSurface(lcw)
            if lps:
                sdl2.SDL_SetSurfaceBlendMode(lps, sdl2.SDL_BLENDMODE_NONE)
                lsurf   = lps.contents
                lp_i32  = lsurf.pitch // 4
                lp_ptr  = ctypes.cast(lsurf.pixels,
                                      ctypes.POINTER(ctypes.c_int32))
                ld_flat = np.ctypeslib.as_array(lp_ptr,
                                                shape=(PROBE_H * lp_i32,))
                lc_dst       = ld_flat.reshape(PROBE_H, lp_i32)
                lc_window_p  = lcw
                lc_surface_p = lps
                print("EvoCA SDL: lut_complexity window created", flush=True)
            else:
                sdl2.SDL_DestroyWindow(lcw)
        else:
            print("EvoCA SDL: lut_complexity window creation failed", flush=True)

    print("EvoCA SDL: entering main loop", flush=True)

    event = sdl2.SDL_Event()

    while ctrl[0] == 0:

        # Process SDL events
        while sdl2.SDL_PollEvent(ctypes.byref(event)):
            if event.type == sdl2.SDL_QUIT:
                ctrl[0] = 1
            elif event.type == sdl2.SDL_KEYDOWN:
                k = event.key.keysym.sym
                if k in (sdl2.SDLK_q, sdl2.SDLK_ESCAPE):
                    ctrl[0] = 1

        if ctrl[0]:
            break

        # Render main window
        sdl2.SDL_LockSurface(surface_p)
        src = pixels.reshape(N, N)
        if px == 1:
            dst[:N, :N] = src
        else:
            dst[:H, :W] = np.repeat(np.repeat(src, px, axis=0), px, axis=1)
        sdl2.SDL_UnlockSurface(surface_p)
        sdl2.SDL_UpdateWindowSurface(window_p)

        # Render probe windows
        if probe_cursor is not None:
            cur = int(probe_cursor[0])
            for i in range(len(probe_windows)):
                sdl2.SDL_LockSurface(probe_surfaces[i])
                cm, cb = probe_colors[i]
                _render_probe(probe_dsts[i], probe_dsts[i].shape[1],
                              probe_means[i], probe_stds[i], cur, cm, cb,
                              y_fixed=probe_yfix[i])
                sdl2.SDL_UnlockSurface(probe_surfaces[i])
                sdl2.SDL_UpdateWindowSurface(probe_windows[i])

        # Render activity window (scroll so newest column is on the right)
        if act_window_p is not None and activity_pixels is not None:
            sdl2.SDL_LockSurface(act_surface_p)
            cur_act = int(activity_cursor[0])
            act_dst[:ACT_H, :PROBE_W] = np.roll(activity_pixels, -cur_act,
                                                  axis=1)
            sdl2.SDL_UnlockSurface(act_surface_p)
            sdl2.SDL_UpdateWindowSurface(act_window_p)

        # Render cgenom activity window
        if cg_act_window_p is not None and cg_activity_pixels is not None:
            sdl2.SDL_LockSurface(cg_act_surface_p)
            cur_cga = int(cg_activity_cursor[0])
            cg_act_dst[:ACT_H, :PROBE_W] = np.roll(cg_activity_pixels,
                                                     -cur_cga, axis=1)
            sdl2.SDL_UnlockSurface(cg_act_surface_p)
            sdl2.SDL_UpdateWindowSurface(cg_act_window_p)

        # Render LUT complexity window
        if lc_window_p is not None and lut_complexity_pixels is not None:
            sdl2.SDL_LockSurface(lc_surface_p)
            cur_lc = int(lut_complexity_cursor[0])
            lc_dst[:PROBE_H, :PROBE_W] = np.roll(lut_complexity_pixels,
                                                   -cur_lc, axis=1)
            sdl2.SDL_UnlockSurface(lc_surface_p)
            sdl2.SDL_UpdateWindowSurface(lc_window_p)

        # Window title
        step   = int(ctrl[2])
        mode   = COLOR_MODES[min(int(ctrl[1]), 3)]
        paused = bool(ctrl[4])
        if paused:
            title = f"EvoCA  PAUSED  t={step}  color={mode}"
            sdl2.SDL_Delay(16)
        else:
            fps   = ctrl[3] / 10.0
            title = f"EvoCA  t={step}  fps={fps:.1f}  color={mode}"
            sdl2.SDL_Delay(4)   # cap render ~250 fps; prevents 100% CPU spin
        sdl2.SDL_SetWindowTitle(window_p, title.encode())

    print("EvoCA SDL: exiting cleanly", flush=True)
    if lc_window_p is not None:
        sdl2.SDL_DestroyWindow(lc_window_p)
    if cg_act_window_p is not None:
        sdl2.SDL_DestroyWindow(cg_act_window_p)
    if act_window_p is not None:
        sdl2.SDL_DestroyWindow(act_window_p)
    for pw in probe_windows:
        sdl2.SDL_DestroyWindow(pw)
    sdl2.SDL_DestroyWindow(window_p)
    sdl2.SDL_Quit()
    pixel_shm.close()
    ctrl_shm.close()
    if probe_shm is not None:
        probe_shm.close()
    if activity_shm is not None:
        activity_shm.close()
    if cg_activity_shm is not None:
        cg_activity_shm.close()
    if lut_complexity_shm is not None:
        lut_complexity_shm.close()


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
