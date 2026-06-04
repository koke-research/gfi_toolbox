"""
gfi_engine.py
=============
GFI Engine — core hydrological and topological computations.

Public API
----------
raster_to_array            QgsRasterLayer -> (array, cellsize, bounds)
build_continuous_channel   channel mask from CI threshold + downstream propagation
compute_strahler_order     Strahler order internally (no GRASS/TopoToolbox needed)
get_flow_moves             encoding-agnostic D8 lookup table
trace_flow_step            single downstream step along FlowDir
hillslope_to_river         maps every hillslope pixel to its nearest channel pixel
river_to_confluence        maps every channel pixel to its next downstream confluence
compute_gfi_v1             GFI version 1.0  =  ln(hr / H)
calibrate_gfi              ROC calibration -> optimal threshold, a, MargArea
compute_gfi_v2             GFI version 2.0  (hierarchical backwater update)
compute_performance        confusion-matrix metrics (TPR, TNR, CSI, Kappa, AUC)
validate_gfi_inputs        pre-flight sanity checks
"""

import numpy as np
from collections import deque


# =============================================================================
# RASTER I/O
# =============================================================================

def raster_to_array(qgs_layer, feedback=None):
    """
    Convert a QgsRasterLayer to a float32 NumPy array.

    FIX: the previous version called provider.block(1, extent) without
    dimensions and cast the raw QByteArray bytes to float32, producing
    garbage values (3e+348, 0, etc.).  The correct call is
    provider.block(band, extent, cols, rows) followed by block.value(r, c).

    Returns
    -------
    data     : np.float32 2-D array  (NaN where NoData)
    cellsize : float
    bounds   : dict  {xmin, ymax, cellsize}
    """
    provider = qgs_layer.dataProvider()
    extent   = provider.extent()
    rows     = qgs_layer.height()
    cols     = qgs_layer.width()

    if feedback:
        feedback.pushInfo(f"  Reading {rows}×{cols} raster…")

    block = provider.block(1, extent, cols, rows)

    data = np.empty((rows, cols), dtype=np.float64)
    for i in range(rows):
        for j in range(cols):
            data[i, j] = block.value(i, j)
    data = data.astype(np.float32)

    if feedback:
        feedback.pushInfo(f"  Raw range: [{np.nanmin(data):.4g}, {np.nanmax(data):.4g}]")

    nodata_val = provider.sourceNoDataValue(1)
    if nodata_val is not None and not np.isnan(float(nodata_val)):
        mask = np.isclose(data, float(nodata_val), rtol=1e-5, atol=1e-5)
        data[mask] = np.nan
        if feedback:
            feedback.pushInfo(f"  NoData masked: {np.count_nonzero(mask)} px")

    data[data >  1e20] = np.nan
    data[data < -1e10] = np.nan

    cellsize = (extent.xMaximum() - extent.xMinimum()) / cols

    if feedback:
        feedback.pushInfo(f"  Cell size: {cellsize:.4g}  valid px: {np.count_nonzero(~np.isnan(data))}")

    bounds = {'xmin': extent.xMinimum(), 'ymax': extent.yMaximum(), 'cellsize': cellsize}
    return data, cellsize, bounds


# =============================================================================
# FLOW DIRECTION UTILITIES
# =============================================================================

def get_flow_moves(FlowDir_clean):
    """
    Auto-detect GRASS (1-8) or ESRI D8 (powers-of-2) encoding.

    Returns
    -------
    moves      : dict  fd_value -> (dr, dc)
    diag_keys  : set of fd values that are diagonal (distance = cellsize*√2)
    FlowDir_rt : int16 array ready for indexing
    encoding   : str label
    """
    valid = FlowDir_clean[FlowDir_clean != 0]
    if len(valid) == 0:
        return None, None, FlowDir_clean.astype(np.int16), "UNKNOWN"

    if np.max(np.abs(valid)) <= 8:
        encoding   = "GRASS r.watershed (1-8)"
        FlowDir_rt = np.abs(FlowDir_clean).astype(np.int16)
        moves = {
            1: (-1,  1),  # NE
            2: (-1,  0),  # N
            3: (-1, -1),  # NW
            4: ( 0, -1),  # W
            5: ( 1, -1),  # SW
            6: ( 1,  0),  # S
            7: ( 1,  1),  # SE
            8: ( 0,  1),  # E
        }
        diag_keys = {1, 3, 5, 7}
    else:
        encoding   = "ESRI D8 (1-128)"
        FlowDir_rt = FlowDir_clean.astype(np.int16)
        moves = {
            1:   ( 0,  1),  # E
            2:   ( 1,  1),  # SE
            4:   ( 1,  0),  # S
            8:   ( 1, -1),  # SW
            16:  ( 0, -1),  # W
            32:  (-1, -1),  # NW
            64:  (-1,  0),  # N
            128: (-1,  1),  # NE
        }
        diag_keys = {2, 8, 32, 128}

    return moves, diag_keys, FlowDir_rt, encoding


def trace_flow_step(FlowDir_rt, moves, diag_keys, r, c, Ld, cellsize):
    """
    One downstream step from pixel (r, c).

    Returns (nr, nc, Ld_new) on success, or (None, None, Ld) at a sink.
    """
    fd = int(FlowDir_rt[r, c])
    if fd not in moves:
        return None, None, Ld
    dr, dc = moves[fd]
    dist   = cellsize * (1.4142135623730951 if fd in diag_keys else 1.0)
    return r + dr, c + dc, Ld + dist


# =============================================================================
# CHANNEL NETWORK
# =============================================================================

def build_continuous_channel(FlowAcc, FlowDir, Slope, cellsize, th, feedback=None):
    """
    Build a continuous binary channel mask.

    1. CI = A [m²] × tan(S)^1.7  ;  seed pixels where CI > th
    2. Trace each seed downstream, filling any gaps.

    Returns
    -------
    channel : int8 2-D array  (1 = channel, 0 = hillslope)
    """
    rows, cols = FlowDir.shape
    if feedback:
        feedback.pushInfo("  Building channel network…")

    FD  = np.nan_to_num(FlowDir, nan=0.0)
    FA  = np.nan_to_num(FlowAcc, nan=0.0)
    SL  = np.nan_to_num(Slope,   nan=0.001)

    Sl_tan       = np.tan(np.radians(np.clip(SL, 0.0, 89.0)))
    Sl_tan       = np.where(Sl_tan < 1e-6, 1e-6, Sl_tan)
    ci           = FA * (cellsize ** 2) * (Sl_tan ** 1.7)
    initial_mask = ci > th

    if feedback:
        feedback.pushInfo(f"    CI range: [{ci.min():.2e}, {ci.max():.2e}]  th={th:.2e}")
        feedback.pushInfo(f"    Seed pixels: {np.count_nonzero(initial_mask)}")

    moves, diag_keys, FD_rt, enc = get_flow_moves(FD)
    if feedback:
        feedback.pushInfo(f"    Encoding: {enc}")
    if moves is None:
        return np.zeros((rows, cols), dtype=np.int8)

    channel   = np.zeros((rows, cols), dtype=np.int8)
    max_steps = rows + cols

    for r0, c0 in zip(*np.where(initial_mask)):
        r, c, steps = int(r0), int(c0), 0
        while 0 <= r < rows and 0 <= c < cols and steps < max_steps:
            if channel[r, c] == 1 and (r != r0 or c != c0):
                break
            channel[r, c] = 1
            nr, nc, _ = trace_flow_step(FD_rt, moves, diag_keys, r, c, 0.0, cellsize)
            if nr is None:
                break
            r, c = nr, nc
            steps += 1

    if feedback:
        feedback.pushInfo(f"    Channel pixels: {np.count_nonzero(channel)}")
    return channel


# =============================================================================
# STRAHLER ORDER
# =============================================================================

def compute_strahler_order(channel, FlowDir, feedback=None):
    """
    Compute Strahler stream order for every channel pixel without external libs.

    Algorithm: topological sort (BFS from headwaters toward outlet).
      - headwater     → order 1
      - all tributaries equal max  → max + 1
      - otherwise                  → max

    Returns
    -------
    strahler : int16 2-D array  (0 outside channel)
    """
    rows, cols = channel.shape
    if feedback:
        feedback.pushInfo("  Computing Strahler order…")

    FD             = np.nan_to_num(FlowDir, nan=0.0)
    moves, diag_keys, FD_rt, _ = get_flow_moves(FD)
    if moves is None:
        return np.zeros((rows, cols), dtype=np.int16)

    ch_rows, ch_cols = np.where(channel == 1)

    downstream_r = np.full((rows, cols), -1, dtype=np.int32)
    downstream_c = np.full((rows, cols), -1, dtype=np.int32)
    in_count     = np.zeros((rows, cols), dtype=np.int16)

    for r, c in zip(ch_rows, ch_cols):
        nr, nc, _ = trace_flow_step(FD_rt, moves, diag_keys, r, c, 0.0, 1.0)
        if nr is not None and 0 <= nr < rows and 0 <= nc < cols and channel[nr, nc] == 1:
            downstream_r[r, c] = nr
            downstream_c[r, c] = nc
            in_count[nr, nc] += 1

    strahler      = np.zeros((rows, cols), dtype=np.int16)
    max_order_in  = np.zeros((rows, cols), dtype=np.int16)
    max2_order_in = np.zeros((rows, cols), dtype=np.int16)
    pending       = in_count.copy()

    queue = deque()
    for r, c in zip(ch_rows, ch_cols):
        if in_count[r, c] == 0:
            strahler[r, c] = 1
            queue.append((r, c))

    while queue:
        r, c = queue.popleft()
        nr, nc = int(downstream_r[r, c]), int(downstream_c[r, c])
        if nr == -1:
            continue
        s = int(strahler[r, c])
        if s > max_order_in[nr, nc]:
            max2_order_in[nr, nc] = max_order_in[nr, nc]
            max_order_in[nr, nc]  = s
        elif s > max2_order_in[nr, nc]:
            max2_order_in[nr, nc] = s
        pending[nr, nc] -= 1
        if pending[nr, nc] == 0:
            m1 = int(max_order_in[nr, nc])
            m2 = int(max2_order_in[nr, nc])
            strahler[nr, nc] = m1 + 1 if m1 == m2 else m1
            queue.append((nr, nc))

    if feedback:
        feedback.pushInfo(f"    Max Strahler order: {int(np.max(strahler))}")
    return strahler


# =============================================================================
# HILLSLOPE-TO-RIVER MAPPING
# =============================================================================

def hillslope_to_river(DEM, channel, FlowDir, cellsize, feedback=None):
    """
    For every non-channel pixel, trace downstream to the first channel pixel
    and record its (row, col).

    Returns
    -------
    ROW_channel : float32 2-D array  (NaN where unreachable)
    COL_channel : float32 2-D array
    """
    rows, cols = channel.shape
    if feedback:
        feedback.pushInfo("  Hillslope-to-river mapping…")

    FD             = np.nan_to_num(FlowDir, nan=0.0)
    moves, diag_keys, FD_rt, _ = get_flow_moves(FD)
    if moves is None:
        nan_g = np.full((rows, cols), np.nan, dtype=np.float32)
        return nan_g, nan_g.copy()

    ROW_ch    = np.full((rows, cols), np.nan, dtype=np.float32)
    COL_ch    = np.full((rows, cols), np.nan, dtype=np.float32)
    max_steps = rows + cols

    for i in range(1, rows - 1):
        for j in range(1, cols - 1):
            if channel[i, j] != 0 or np.isnan(DEM[i, j]):
                continue
            r, c, steps = i, j, 0
            while 0 < r < rows - 1 and 0 < c < cols - 1 and steps < max_steps:
                if channel[r, c] == 1:
                    ROW_ch[i, j] = r
                    COL_ch[i, j] = c
                    break
                nr, nc, _ = trace_flow_step(FD_rt, moves, diag_keys, r, c, 0.0, cellsize)
                if nr is None:
                    break
                r, c = nr, nc
                steps += 1

    if feedback:
        feedback.pushInfo(f"    Hillslope pixels mapped: {np.count_nonzero(~np.isnan(ROW_ch))}")
    return ROW_ch, COL_ch


# =============================================================================
# RIVER-TO-CONFLUENCE MAPPING
# =============================================================================

def river_to_confluence(channel, strahler, FlowDir, cellsize, feedback=None):
    """
    For every channel pixel, trace downstream until Strahler order increases
    (= confluence) and record that pixel's (row, col).

    Returns
    -------
    ROW_confluence : float32 2-D array  (NaN for outlet-order pixels)
    COL_confluence : float32 2-D array
    """
    rows, cols = channel.shape
    if feedback:
        feedback.pushInfo("  River-to-confluence mapping…")

    FD             = np.nan_to_num(FlowDir, nan=0.0)
    moves, diag_keys, FD_rt, _ = get_flow_moves(FD)
    if moves is None:
        nan_g = np.full((rows, cols), np.nan, dtype=np.float32)
        return nan_g, nan_g.copy()

    max_order      = int(np.max(strahler))
    ROW_conf       = np.full((rows, cols), np.nan, dtype=np.float32)
    COL_conf       = np.full((rows, cols), np.nan, dtype=np.float32)
    max_steps      = rows + cols

    for i, j in zip(*np.where(channel == 1)):
        s0 = int(strahler[i, j])
        if s0 == 0 or s0 >= max_order:
            continue
        r, c, steps = i, j, 0
        while 0 < r < rows - 1 and 0 < c < cols - 1 and steps < max_steps:
            nr, nc, _ = trace_flow_step(FD_rt, moves, diag_keys, r, c, 0.0, cellsize)
            if nr is None:
                break
            if strahler[nr, nc] > s0:
                ROW_conf[i, j] = nr
                COL_conf[i, j] = nc
                break
            if strahler[nr, nc] == s0:
                r, c = nr, nc
                steps += 1
            else:
                ROW_conf[i, j] = r
                COL_conf[i, j] = c
                break

    if feedback:
        feedback.pushInfo(f"    Channel pixels with confluence: {np.count_nonzero(~np.isnan(ROW_conf))}")
    return ROW_conf, COL_conf


# =============================================================================
# GFI 1.0
# =============================================================================

def compute_gfi_v1(DEM, FlowAcc, channel, ROW_channel, COL_channel,
                   cellsize, n, feedback=None):
    """
    GFI version 1.0  =  ln(hr / H)

    Formulation (Samela et al. 2017)
    ---------------------------------
    hr = ( (A_river + 1) * cellsize² / 10⁶ ) ^ n
    H  = DEM(pixel) - DEM(nearest channel pixel)      [≥ 0.001 m floor]
    GFI = ln(hr / H)

    For channel pixels themselves, H is set to 0.001 m (effectively in-channel).

    Parameters
    ----------
    DEM         : float32 2-D array
    FlowAcc     : float32 2-D array  (cells or area — will be converted via cellsize)
    channel     : int8 2-D array
    ROW_channel : float32 2-D array from hillslope_to_river
    COL_channel : float32 2-D array from hillslope_to_river
    cellsize    : float
    n           : float  (calibration exponent, default 0.3544)

    Returns
    -------
    GFI_v1    : float32 2-D array
    Hr_matrix : float32 2-D array  (reference water depth at every pixel)
    H_matrix  : float32 2-D array  (hillslope height above channel)
    Ariver    : float32 2-D array  (flow accumulation of mapped channel pixel)
    """
    rows, cols = DEM.shape
    if feedback:
        feedback.pushInfo("  Computing GFI 1.0…")

    valid_idx = ~np.isnan(ROW_channel)  # hillslope pixels that reached a channel
    R_idx = ROW_channel[valid_idx].astype(int)
    C_idx = COL_channel[valid_idx].astype(int)

    # H — height above nearest channel pixel
    H = np.full((rows, cols), np.nan, dtype=np.float32)
    H[valid_idx] = DEM[valid_idx] - DEM[R_idx, C_idx]
    H[channel > 0] = 0.001          # channel pixels: H → ~0
    H[H <= 0]      = 0.001          # clamp depressions

    # Ariver — flow accumulation of mapped channel pixel
    Ariver = np.full((rows, cols), np.nan, dtype=np.float32)
    Ariver[valid_idx]   = FlowAcc[R_idx, C_idx]
    Ariver[channel > 0] = FlowAcc[channel > 0]

    # hr — reference water depth
    Hr = (((Ariver + 1.0) * (cellsize ** 2)) / 1_000_000.0) ** n

    # GFI 1.0
    with np.errstate(divide='ignore', invalid='ignore'):
        GFI_v1 = np.where(
            (Hr > 0) & (H > 0),
            np.log(Hr / H).astype(np.float32),
            np.nan
        ).astype(np.float32)

    if feedback:
        finite = GFI_v1[np.isfinite(GFI_v1)]
        feedback.pushInfo(f"    GFI 1.0 range: [{finite.min():.3f}, {finite.max():.3f}]")

    return GFI_v1, Hr, H, Ariver


# =============================================================================
# CALIBRATION  — method: ROCcurve_maggiore (MATLAB translation)
# =============================================================================

def _auc_maggiore(fpr_arr, tpr_arr):
    """
    AUC via the rectangular averaging method from areaundercurve.m.

    Sorts by FPR, appends the (1,1) closing point, then averages the
    upper (right-edge) and lower (left-edge) rectangular sums.
    """
    sort_idx = np.argsort(fpr_arr)
    x = np.concatenate([fpr_arr[sort_idx], [1.0]])
    y = np.concatenate([tpr_arr[sort_idx], [1.0]])

    x_diff        = np.diff(x)
    x_diff_full   = np.concatenate([[x[0]], x_diff])   # prepend first x value

    auc_upper = float(np.sum(y * x_diff_full))                          # upper rectangles
    auc_lower = float(np.sum(np.concatenate([[0.0], y[:-1]]) * x_diff_full))  # lower rectangles
    return (auc_upper + auc_lower) / 2.0


def calibrate_gfi(GFI, flood_map, ROW_channel, COL_channel, channel,
                  version_label="1.0", roc_steps=0.01, feedback=None):
    """
    ROC calibration following ROCcurve_maggiore.m (Samela / Manfreda group).

    Algorithm
    ---------
    1. Build MargArea: project flood_map from each hillslope pixel's nearest
       channel pixel.  MargArea acts as SPATIAL MASK — pixels where MargArea
       is not NaN define the evaluation domain (marginal hazard area).
    2. Normalise GFI to [-1, 1]  (identical to MATLAB matrix_norm).
    3. Sweep thresholds t from -1 to +1 in steps of roc_steps.
    4. At each t: classify pixel as flooded if GFI_norm >= t.
    5. Compute FP, FN, VP, VN inside the MargArea mask.
    6. Compute FPR and TPR; optimise by minimising  FPR + FNR  (= FPR + (1-TPR)),
       i.e. the distance to the perfect point (0, 1) on the ROC plane.
    7. Denormalise the optimal normalised threshold back to original GFI units.
    8. Derive  a = 1 / exp(t_original)   (= exp(-t_original)).
    9. AUC via the rectangular averaging method (areaundercurve.m).

    Parameters
    ----------
    GFI           : float32 2-D array
    flood_map     : float32 2-D array, binary (0 / 1)
    ROW_channel   : float32 2-D array from hillslope_to_river
    COL_channel   : float32 2-D array from hillslope_to_river
    channel       : int8   2-D array
    version_label : str   for log messages
    roc_steps     : float threshold sweep step (default 0.01; smaller = finer)

    Returns
    -------
    a             : float  calibration coefficient  exp(-t_original)
    t_original    : float  optimal threshold in original GFI units
    MargArea      : float32 2-D array  (spatial evaluation mask)
    roc_data      : dict   {fpr, tpr, thresholds_norm, thresholds_orig,
                             auc, opt_t_norm, opt_fpr, opt_tpr, opt_dist}
    """
    rows, cols = GFI.shape
    if feedback:
        feedback.pushInfo(f"  Calibrating GFI {version_label} "
                          f"(step={roc_steps})…")

    # ── 1. Build MargArea ─────────────────────────────────────────────────────
    valid_idx = ~np.isnan(ROW_channel)
    R_idx = ROW_channel[valid_idx].astype(int)
    C_idx = COL_channel[valid_idx].astype(int)

    MargArea = np.full((rows, cols), np.nan, dtype=np.float32)
    MargArea[valid_idx]   = flood_map[R_idx, C_idx]
    MargArea[channel > 0] = flood_map[channel > 0]

    # Evaluation mask: finite GFI + valid flood truth + inside MargArea
    mask_2d = MargArea > 0
    n_valid = int(np.count_nonzero(mask_2d))

    if feedback:
        n_pos = int(np.count_nonzero(flood_map[mask_2d] == 1))
        n_neg = n_valid - n_pos
        feedback.pushInfo(f"    Evaluation pixels: {n_valid}  "
                          f"(flooded={n_pos}, dry={n_neg})")

    _nan_roc = {'fpr': np.array([0.0, 1.0]), 'tpr': np.array([0.0, 1.0]),
                'thresholds_norm': np.array([np.nan]),
                'thresholds_orig': np.array([np.nan]),
                'auc': np.nan, 'opt_t_norm': np.nan,
                'opt_fpr': np.nan, 'opt_tpr': np.nan, 'opt_dist': np.nan}

    if n_valid < 10:
        if feedback:
            feedback.pushInfo("    WARNING: too few valid pixels — skipping calibration.")
        return 1.0, 0.0, MargArea, _nan_roc

    gfi_vals   = GFI[mask_2d].astype(np.float64)
    truth_vals = flood_map[mask_2d].astype(np.int8)

    P = int(np.sum(truth_vals == 1))
    N = int(np.sum(truth_vals == 0))
    if P == 0 or N == 0:
        if feedback:
            feedback.pushInfo("    WARNING: flood map has only one class — ROC undefined.")
        return 1.0, 0.0, MargArea, _nan_roc

    # ── 2. Normalise GFI to [-1, 1] ──────────────────────────────────────────
    gfi_min = float(np.nanmin(GFI[mask_2d]))
    gfi_max = float(np.nanmax(GFI[mask_2d]))
    if gfi_max == gfi_min:
        if feedback:
            feedback.pushInfo("    WARNING: GFI has no variation — cannot normalise.")
        return 1.0, 0.0, MargArea, _nan_roc

    gfi_norm = 2.0 * ((gfi_vals - gfi_min) / (gfi_max - gfi_min) - 0.5)

    # ── 3-6. Threshold sweep ──────────────────────────────────────────────────
    thresholds_norm = np.arange(-1.0, 1.0 + roc_steps, roc_steps)

    fpr_arr  = np.empty(len(thresholds_norm), dtype=np.float64)
    tpr_arr  = np.empty(len(thresholds_norm), dtype=np.float64)
    fnr_arr  = np.empty(len(thresholds_norm), dtype=np.float64)

    F_optim  = np.inf
    opt_t_norm  = thresholds_norm[0]
    opt_t_orig  = 0.0
    opt_fpr     = 1.0
    opt_tpr     = 0.0

    for i, t in enumerate(thresholds_norm):
        predicted = (gfi_norm >= t).astype(np.int8)

        vp = int(np.sum((predicted == 1) & (truth_vals == 1)))  # TP
        vn = int(np.sum((predicted == 0) & (truth_vals == 0)))  # TN
        fp = int(np.sum((predicted == 1) & (truth_vals == 0)))  # FP
        fn = int(np.sum((predicted == 0) & (truth_vals == 1)))  # FN

        fpr = fp / (fp + vn) if (fp + vn) > 0 else 0.0
        fnr = fn / (fn + vp) if (fn + vp) > 0 else 0.0
        tpr = 1.0 - fnr

        fpr_arr[i] = fpr
        tpr_arr[i] = tpr
        fnr_arr[i] = fnr

        # Optimise: minimise FPR + FNR  (distance to perfect point (0,1))
        dist = fpr + fnr
        if dist < F_optim:
            F_optim    = dist
            opt_t_norm = t
            opt_fpr    = fpr
            opt_tpr    = tpr
            # Denormalise: t_orig = ((t_norm + 1) / 2) * (max - min) + min
            opt_t_orig = ((t + 1.0) / 2.0) * (gfi_max - gfi_min) + gfi_min

    # ── 7-8. Calibration coefficient ─────────────────────────────────────────
    a = float(np.exp(-opt_t_orig))

    # ── 9. AUC (Maggiore rectangular method) ─────────────────────────────────
    auc_val = _auc_maggiore(fpr_arr, tpr_arr)

    # Denormalise all thresholds for reference
    thresholds_orig = ((thresholds_norm + 1.0) / 2.0) * (gfi_max - gfi_min) + gfi_min

    roc_data = {
        'fpr':              fpr_arr,
        'tpr':              tpr_arr,
        'thresholds_norm':  thresholds_norm,
        'thresholds_orig':  thresholds_orig,
        'auc':              auc_val,
        'opt_t_norm':       opt_t_norm,
        'opt_fpr':          opt_fpr,
        'opt_tpr':          opt_tpr,
        'opt_dist':         F_optim,
    }

    if feedback:
        feedback.pushInfo(
            f"    GFI {version_label} — AUC={auc_val:.4f}  "
            f"opt_t_norm={opt_t_norm:.4f}  opt_t_orig={opt_t_orig:.4f}  "
            f"a={a:.6f}  FPR={opt_fpr:.4f}  TPR={opt_tpr:.4f}")

    return a, opt_t_orig, MargArea, roc_data


# =============================================================================
# PERFORMANCE METRICS
# =============================================================================

def compute_performance(flood_pred, flood_map, MargArea, version_label="1.0",
                        roc_data=None, feedback=None):
    """
    Confusion-matrix metrics for a binary flood prediction raster.

    Evaluation domain: pixels inside MargArea mask with valid flood truth.
    Truth labels come directly from flood_map (not MargArea).

    Metrics
    -------
    TPR  / Sensitivity   = TP / (TP + FN)
    TNR  / Specificity   = TN / (TN + FP)
    FPR  / Fall-out      = FP / (FP + TN)  = 1 - TNR
    FNR  / Miss rate     = FN / (FN + TP)  = 1 - TPR
    CSI  / Threat Score  = TP / (TP + FP + FN)
    F1                   = 2·TP / (2·TP + FP + FN)
    Kappa (Cohen's)      = (p_o - p_e) / (1 - p_e)
    Bias                 = (TP + FP) / (TP + FN)
    AUC                  = from roc_data

    Parameters
    ----------
    flood_pred    : float32 2-D array  (binary prediction: 1 = flooded, 0 = dry)
    flood_map     : float32 2-D array  binary truth (0 / 1)
    MargArea      : float32 2-D array  spatial mask (not used as truth)
    version_label : str
    roc_data      : dict from calibrate_gfi

    Returns
    -------
    metrics : dict
    """
    valid = ~np.isnan(MargArea) & ~np.isnan(flood_map) & ~np.isnan(flood_pred)

    pred  = flood_pred[valid].astype(int)
    truth = flood_map[valid].astype(int)

    TP    = int(np.sum((pred == 1) & (truth == 1)))
    TN    = int(np.sum((pred == 0) & (truth == 0)))
    FP    = int(np.sum((pred == 1) & (truth == 0)))
    FN    = int(np.sum((pred == 0) & (truth == 1)))
    total = TP + TN + FP + FN

    def _r(num, den): return round(num / den, 4) if den > 0 else np.nan

    TPR   = _r(TP, TP + FN)
    TNR   = _r(TN, TN + FP)
    FPR   = _r(FP, FP + TN)
    FNR   = _r(FN, FN + TP)
    CSI   = _r(TP, TP + FP + FN)
    F1    = _r(2 * TP, 2 * TP + FP + FN)
    Bias  = _r(TP + FP, TP + FN)

    p_o   = (TP + TN) / total if total > 0 else np.nan
    p_e   = (((TP + FP) * (TP + FN)) + ((TN + FN) * (TN + FP))) / (total ** 2) \
            if total > 0 else np.nan
    Kappa = round((p_o - p_e) / (1.0 - p_e), 4) \
            if (p_e is not None and not np.isnan(p_e) and p_e < 1.0) else np.nan

    AUC   = round(float(roc_data['auc']), 4) \
            if (roc_data is not None and not np.isnan(roc_data['auc'])) else np.nan

    metrics = {
        'version':           version_label,
        'TP': TP, 'TN': TN, 'FP': FP, 'FN': FN,
        'TPR_sensitivity':   TPR,
        'TNR_specificity':   TNR,
        'FPR':               FPR,
        'FNR':               FNR,
        'CSI_threat_score':  CSI,
        'F1_score':          F1,
        'Kappa':             Kappa,
        'Bias':              Bias,
        'AUC':               AUC,
    }

    if feedback:
        feedback.pushInfo(
            f"  GFI {version_label} — AUC={AUC}  CSI={CSI}  "
            f"Kappa={Kappa}  F1={F1}  "
            f"TPR={TPR}  TNR={TNR}  FPR={FPR}  FNR={FNR}")

    return metrics


# =============================================================================
# GFI 2.0  — hierarchical backwater update
# =============================================================================

def compute_gfi_v2(DEM, FlowAcc, channel, strahler,
                   ROW_channel, COL_channel,
                   ROW_confluence, COL_confluence,
                   cellsize, n, a_GFIv1,
                   max_iter=6, feedback=None):
    """
    GFI version 2.0 — hierarchical confluence backwater propagation.

    Algorithm (Albertini et al. 2022)
    -----------------------------------
    For each channel pixel, check whether backwater from the next confluence
    (higher Strahler order) produces a deeper flood depth than the local hr.
    If so, adopt the confluence's hydraulic properties and look further
    downstream.  Iterate up to max_iter times.

    After updating channel hydraulics, re-map to hillslope pixels as in v1.

    Parameters
    ----------
    a_GFIv1  : float  calibration coefficient from GFI v1 ROC optimisation
                       (a = exp(-optimal_threshold))
    max_iter : int     number of upstream propagation iterations (default 6)

    Returns
    -------
    GFI_v2     : float32 2-D array
    Hr_v2      : float32 2-D array
    H_v2       : float32 2-D array
    Ariver_v2  : float32 2-D array
    """
    rows, cols = DEM.shape
    if feedback:
        feedback.pushInfo(f"  Computing GFI 2.0 (max_iter={max_iter})…")

    # Initialise river-network arrays from v1
    Ariver_net    = FlowAcc.copy().astype(np.float32)
    DEM_river_net = DEM.copy().astype(np.float32)

    hr_net = (((Ariver_net + 1.0) * (cellsize ** 2)) / 1_000_000.0) ** n
    # Initial water depth estimate (used only for comparison inside loop)
    WD_net = np.maximum(0.0, hr_net * a_GFIv1 - 0.001)

    # Working confluence pointers (will be updated iteration by iteration)
    curr_R = ROW_confluence.copy()
    curr_C = COL_confluence.copy()

    idx_channel = list(zip(*np.where(channel > 0)))

    for k in range(max_iter):
        updated = 0
        for r, c in idx_channel:
            nr = curr_R[r, c]
            nc = curr_C[r, c]
            if np.isnan(nr):
                continue
            nr, nc = int(nr), int(nc)

            A_next  = float(FlowAcc[nr, nc])
            hr_next = (((A_next + 1.0) * (cellsize ** 2)) / 1_000_000.0) ** n

            H_to_next = float(DEM[r, c]) - float(DEM[nr, nc])
            if H_to_next <= 0:
                H_to_next = 0.001

            WD_potential = max(0.0, hr_next * a_GFIv1 - H_to_next)

            if WD_potential > WD_net[r, c]:
                Ariver_net[r, c]    = A_next
                DEM_river_net[r, c] = float(DEM[nr, nc])
                WD_net[r, c]        = WD_potential
                # Follow the confluence chain one level further
                curr_R[r, c] = ROW_confluence[nr, nc]
                curr_C[r, c] = COL_confluence[nr, nc]
                updated += 1

        if feedback:
            feedback.pushInfo(f"    Iteration {k+1}/{max_iter}: {updated} pixels updated")
        if updated == 0:
            break   # converged early

    # Map updated channel hydraulics back to hillslope pixels
    valid_idx = ~np.isnan(ROW_channel)
    R_idx     = ROW_channel[valid_idx].astype(int)
    C_idx     = COL_channel[valid_idx].astype(int)

    Ariver_v2 = np.full((rows, cols), np.nan, dtype=np.float32)
    H_v2      = np.full((rows, cols), np.nan, dtype=np.float32)

    Ariver_v2[valid_idx]   = Ariver_net[R_idx, C_idx]
    Ariver_v2[channel > 0] = Ariver_net[channel > 0]

    H_v2[valid_idx]   = DEM[valid_idx] - DEM_river_net[R_idx, C_idx]
    H_v2[channel > 0] = DEM[channel > 0] - DEM_river_net[channel > 0]
    H_v2[H_v2 <= 0]   = 0.001

    Hr_v2 = (((Ariver_v2 + 1.0) * (cellsize ** 2)) / 1_000_000.0) ** n

    with np.errstate(divide='ignore', invalid='ignore'):
        GFI_v2 = np.where(
            (Hr_v2 > 0) & (H_v2 > 0),
            np.log(Hr_v2 / H_v2).astype(np.float32),
            np.nan
        ).astype(np.float32)

    if feedback:
        finite = GFI_v2[np.isfinite(GFI_v2)]
        feedback.pushInfo(f"    GFI 2.0 range: [{finite.min():.3f}, {finite.max():.3f}]")

    return GFI_v2, Hr_v2, H_v2, Ariver_v2


# =============================================================================
# INPUT VALIDATION
# =============================================================================

def validate_gfi_inputs(dem, flowacc, fdir, flood_map, slope, feedback=None):
    """
    Pre-flight checks on all raster inputs.
    Returns: (is_valid: bool, message: str)
    """
    issues = []

    if np.all(np.isnan(dem)):
        issues.append("DEM is entirely NaN.")
    elif np.nanmax(dem) == np.nanmin(dem):
        issues.append("DEM has no variation (flat/constant raster).")

    if np.nanmax(flowacc) <= 0:
        issues.append("Flow Accumulation has no positive values.")

    if len(fdir[(fdir != 0) & ~np.isnan(fdir)]) == 0:
        issues.append("Flow Direction raster has no valid values.")

    if flood_map is not None:
        uv = np.unique(flood_map[~np.isnan(flood_map)])
        if not np.all(np.isin(uv, [0, 1])):
            issues.append("Flood Map must be strictly binary (0 / 1).")
        if np.all(flood_map == 0) or np.all(np.isnan(flood_map)):
            issues.append("Flood Map contains no flooded pixels.")

    if np.all(slope == 0) or np.all(np.isnan(slope)):
        issues.append("Slope raster is entirely zero or NaN.")

    if not (dem.shape == flowacc.shape == fdir.shape == slope.shape):
        issues.append(
            f"Shape mismatch — DEM:{dem.shape} ACC:{flowacc.shape} "
            f"DIR:{fdir.shape} SLP:{slope.shape}")

    if issues:
        msg = "Input validation failed:\n" + "\n".join(f"  • {i}" for i in issues)
        if feedback:
            feedback.pushInfo(msg)
        return False, msg

    if feedback:
        feedback.pushInfo("  All inputs passed validation.")
    return True, "OK"
