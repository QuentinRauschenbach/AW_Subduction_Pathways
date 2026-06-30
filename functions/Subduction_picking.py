import numpy as np
import xarray as xr
import logging
import scipy
import matplotlib.pyplot as plt
import matplotlib as mpl
import cmocean as cmo
import pandas as pd
from matplotlib.colors import to_rgba

import sys 
sys.path.append('/albedo/work/user/quraus001/FESOM_particles/analysis/functions/')
from ocean_helper import classify_water_mass

def find_first_subduction(drifter_idx, ds, MIN_STEPS, STEPS_PER_DAY):
    """Find first timestep of MIN_DAYS consecutive subduction conditions."""
    temp      = ds.temp.isel(drifter=drifter_idx).rolling(time=STEPS_PER_DAY, center=True).mean()
    salt      = ds.salt.isel(drifter=drifter_idx).rolling(time=STEPS_PER_DAY, center=True).mean()
    temp_surf = ds.temp_surface.isel(drifter=drifter_idx).rolling(time=STEPS_PER_DAY, center=True).mean()
    salt_surf = ds.salt_surface.isel(drifter=drifter_idx).rolling(time=STEPS_PER_DAY, center=True).mean()

    labels_drifter = classify_water_mass(temp, salt, return_numeric=False)
    labels_surface = classify_water_mass(temp_surf, salt_surf, return_numeric=False)

    mask = (
        ((labels_drifter == 'AW') | (labels_drifter == 'Mixed')) &
        ((labels_surface == 'PW') | (labels_surface == 'MW'))
    )

    count = 0
    for t in range(len(mask)):
        count = count + 1 if mask[t] else 0
        if count >= MIN_STEPS:
            return t - MIN_STEPS + 1

    return None

# NOt used at the moment
def find_last_ML(drifter_idx, ds, MIN_STEPS, STEPS_PER_DAY):
    """Find first timestep of the last MIN_DAYS in of the drifter in the mixed layer."""
    depth     = ds.depth.isel(drifter=drifter_idx).rolling(time=STEPS_PER_DAY, center=True).mean()
    mld       = ds.mld.isel(drifter=drifter_idx).rolling(time=STEPS_PER_DAY, center=True).mean()
    
    # Create a boolean mask where True means the drifter is within the mixed layer
    mask = depth <= mld

    count = 0
    # Loop backward through the time series
    for t in range(len(mask) - 1, -1, -1):
        # xarray's rolling mean leaves NaNs at the edges; we treat NaNs as False
        is_in_ml = bool(mask[t].values) if not mask[t].isnull() else False
        
        count = count + 1 if is_in_ml else 0
        
        # Once we find the required consecutive steps moving backward,
        # 't' is the starting index of this block.
        if count >= MIN_STEPS:
            return t

    return None

def pick_drop_idx(pw_time_idx, ds, drifter_idx, STEPS_PER_DAY, logger,
                  SEARCH_BUFFER_DAYS  = 90,     # days to still search after pw_time_idx
                    ROLLING_DAYS        = 30,     # days for rolling mean
                    ASCENT_PENALTY      = 5.0,    # penalty weight for upward motion
                    SIGMA_DAYS          = 60,      # Damping in time
                    SIGNIFICANCE_MARGIN = 10.0,   # m — particle must not resurface within this of start depth
                    #MIN_NET_GAIN        = 20.0,   # m — minimum net depth gain (currently unused in gate)
                    SETTLE_OFFSET_DAYS  = 15,     # days after center before measuring settled depth
                    SETTLE_WINDOW_DAYS  = 60,     # max days after center for settle window
                    W_SEARCH_HALF_DAYS  = 7      # half-width (days) for daily-w peak search):
):
    """
    Find the start and center of the subduction drop for a single drifter.

    Strategy:
      1. Compute forward/backward cumulative scores (with ascent penalty) on
         a 30-day rolling-smoothed depth signal.
      2. Combine into a geometric-mean score, damp toward pw_time_idx with a
         Gaussian, and find the best candidate peak.
      3. Walk from the damped peak → undamped peak → second-best damped peak
         until the event passes a post-descent significance check.
      4. Refine the center to the highest daily vertical velocity near the peak,
         then return (mean_drop, idx_start).

    Parameters
    ----------
    pw_time_idx : int
        Time index at which polar water is confirmed above the particle.
    ds : xr.Dataset
        Full dataset (must contain 'depth').
    drifter_idx : int
        Index along the drifter dimension.

    Returns
    -------
    mean_drop : int or None
        Time index of the drop center (highest daily w near peak), or None if
        no significant event was found.
    idx_start : int or None
        Time index of the start of the descent, or None if not significant.
    """

    # ------------------------------------------------------------------ #
    # Parameters                                                           #
    # ------------------------------------------------------------------ #

    spd = STEPS_PER_DAY  # alias for brevity

    ds_drifter  = ds.isel(drifter=drifter_idx)
    end_idx     = min(pw_time_idx + SEARCH_BUFFER_DAYS * spd, len(ds.time))
    search      = ds_drifter.isel(time=slice(0, end_idx))

    # ------------------------------------------------------------------ #
    # Smoothed depth signals                                               #
    # ------------------------------------------------------------------ #
    depth_smooth = (
        search.depth
        .rolling(time=ROLLING_DAYS * spd, center=True)
        .mean()
        .interpolate_na(dim="time", method="nearest", fill_value="extrapolate")
    )
    depth_daily = (
        search.depth
        .rolling(time=spd, center=True)
        .mean()
        .interpolate_na(dim="time", method="nearest", fill_value="extrapolate")
    )

    # Also keep a full-dataset smooth for the debug plot axis
    depth_rolling_full = ds_drifter.depth.rolling(time=ROLLING_DAYS * spd, center=True).mean()

    # ------------------------------------------------------------------ #
    # Forward / backward penalised velocity arrays                        #
    # ------------------------------------------------------------------ #
    w_fwd      = depth_smooth.diff("time").values          # positive = sinking
    w_fwd_day  = depth_daily.diff("time").values
    w_bwd      = w_fwd[::-1]

    def _penalise(w):
        return np.where(w > 0, w, w * ASCENT_PENALTY)

    pen_fwd = _penalise(w_fwd)
    pen_bwd = _penalise(w_bwd)

    # ------------------------------------------------------------------ #
    # Cumulative reset score                                               #
    # ------------------------------------------------------------------ #
    def _reset_cumsum(w):
        """Cumulative sum that resets to zero whenever it would go negative."""
        run_sum = np.cumsum(w)
        run_min = np.minimum.accumulate(run_sum)
        return run_sum - run_min

    score_fwd = _reset_cumsum(pen_fwd)
    score_bwd = _reset_cumsum(pen_bwd)[::-1]   # flip back to forward time

    combined        = np.sqrt(score_fwd * score_bwd)
    combined_damped = combined * _gaussian_weight(combined, pw_time_idx, SIGMA_DAYS, spd)

    # Backward score with early-time guard
    score_bwd_mod             = score_bwd.copy()
    score_bwd_mod[:(ROLLING_DAYS // 2) * spd] = np.nan

    # Forward score with early-time guard
    score_fwd_mod             = score_fwd.copy()
    score_fwd_mod[-(ROLLING_DAYS // 2) * spd:] = np.nan

    # ------------------------------------------------------------------ #
    # Helper: find true peak of the hump containing idx (1/e threshold)   #
    # ------------------------------------------------------------------ #
    def _hump_peak(score, anchor_idx, threshold_frac=0.37):
        threshold    = score[anchor_idx] * threshold_frac
        labels, _    = scipy.ndimage.label(score >= threshold)
        basin        = np.where(labels == labels[anchor_idx])[0]
        return basin[np.argmax(score[basin])]

    # ------------------------------------------------------------------ #
    # Helper: closest peak in daily-w within ±W_SEARCH_HALF_DAYS of idx  #
    # ------------------------------------------------------------------ #
    def _refine_to_w_peak(idx):
        half  = W_SEARCH_HALF_DAYS * spd
        start = max(0, idx - half)
        end   = min(len(w_fwd_day), idx + half)
        return start + np.argmax(w_fwd_day[start:end])
    
    def _refine_to_depth_min(idx):
        half  = W_SEARCH_HALF_DAYS * spd
        start = max(0, idx)# - half)
        end   = min(len(depth_daily.values), idx + half)
        return start + np.argmin(depth_daily.values[start:end])
    
    def _refine_to_depth_max(idx):
        """Finds deepest point (end of drop) using a centered window."""
        half  = W_SEARCH_HALF_DAYS * spd
        start = max(0, idx - half)  # Centered window handles early/late score peaks
        end   = min(len(depth_daily.values), idx + half)
        return start + np.argmax(depth_daily.values[start:end])
    

    # ------------------------------------------------------------------ #
    # Helper: mask hump of failed_idx and return next-best peak           #
    # ------------------------------------------------------------------ #
    def _next_best_peak(score, failed_idx):
        threshold      = score[failed_idx] * 0.1
        labels, _      = scipy.ndimage.label(score >= threshold)
        new_score      = score.copy()
        failed_label   = labels[failed_idx]
        if failed_label > 0:
            new_score[labels == failed_label] = 0.0
        else:
            new_score[failed_idx] = 0.0
        return np.nanargmax(new_score), new_score

    # ------------------------------------------------------------------ #
    # Helper: post-descent significance check                              #
    # ------------------------------------------------------------------ #
    def _is_significant(idx_start, idx_center):
        settle_start = idx_center + SETTLE_OFFSET_DAYS * spd
        settle_end   = min(
            len(ds_drifter.time) - 1,
            max(pw_time_idx, idx_center + SETTLE_WINDOW_DAYS * spd)
        )
        # Settle window is empty — can't assess significance
        if settle_start >= settle_end:
            logger.debug(f"Drifter {drifter_idx}: empty settle window, marking as not significant")
            return False, np.nan, np.nan
    
        post   = ds_drifter.depth.isel(time=slice(settle_start, settle_end))
        z0     = ds_drifter.depth.isel(time=idx_start).item()
        net    = post.median().item() - z0
        sigma  = post.std().item()
        z_min  = post.min().item()

        is_stable       = net > 2.0 * sigma
        did_not_resurface = z_min > z0 + SIGNIFICANCE_MARGIN

        return is_stable and did_not_resurface, net, sigma

    # ------------------------------------------------------------------ #
    # Helper: warn if particle resurfaces after the picked drop           #
    # ------------------------------------------------------------------ #
    def _check_resurface_warning(mean_drop, idx_start, depth_at_drop, depth_before_drop):
        gap_start = mean_drop + 15 * spd
        if gap_start >= pw_time_idx:
            return ""
        interlude = depth_smooth.isel(time=slice(gap_start, pw_time_idx))
        if interlude.min() < depth_at_drop * 0.9:
            logger.warning(f"Drifter {drifter_idx}: particle resurfaces after drop — consider another peak")
            return "resurface warning"
        if interlude.min() < depth_before_drop:
            logger.warning(f"Drifter {drifter_idx}: particle resurfaces to pre-drop depth — strongly consider another peak")
            return "high resurface warning"

    # ------------------------------------------------------------------ #
    # Debug plot setup                                                     #
    # ------------------------------------------------------------------ #
    if logger.isEnabledFor(logging.DEBUG):
        fig, ax = debug_plot(drifter_idx, ds, pw_time_idx)
        ax.plot(ds.time, depth_rolling_full, color="black", linestyle="--")
        ax_twin = ax.twinx()
    
        ax_twin.fill_between(search.time.values[:-1], score_bwd,       facecolor=to_rgba("darkgreen", alpha=0.1), edgecolor=to_rgba("darkgreen", alpha=1), label="backward score")
        ax_twin.fill_between(search.time.values[:-1], score_fwd,       facecolor=to_rgba("red", alpha=0.1),       edgecolor=to_rgba("red", alpha=1),       label="forward score")
        ax_twin.fill_between(search.time.values[:-1], combined,        facecolor=to_rgba("dodgerblue", alpha=0.1), edgecolor=to_rgba("dodgerblue", alpha=1), label="combined score")
        ax_twin.fill_between(search.time.values[:-1], combined_damped, facecolor=to_rgba("darkblue", alpha=0.1),  edgecolor=to_rgba("darkblue", alpha=1),  label="combined score damped")
        ax_twin.fill_between(search.time.values[:-1], w_fwd_day * 10,  facecolor=to_rgba("magenta", alpha=0.1),   edgecolor=to_rgba("magenta", alpha=1),   label="daily w (×10)")
        
        ax_twin.set_ylim(0, np.nanmax(score_fwd) * 1.2)

    try:
        # ------------------------------------------------------------------ #
        # Peak-finding cascade                                                 #
        # ------------------------------------------------------------------ #

        # --- Attempt 1: damped score ---
        mid_damped      = np.nanargmax(combined_damped)
        idx_center      = _hump_peak(combined, mid_damped)
        idx_start       = _hump_peak(score_bwd_mod, idx_center)
        idx_stop        = _hump_peak(score_fwd_mod, idx_center)
        mean_drop       = _refine_to_w_peak(idx_center)
        idx_start       = _refine_to_depth_min(idx_start)
        idx_stop        = _refine_to_depth_max(idx_stop)
        depth_at_drop   = depth_smooth.values[mean_drop]
        depth_before    = depth_smooth.values[idx_start]

        if logger.isEnabledFor(logging.DEBUG):
            ax.axvline(ds.time.values[mid_damped],  color="darkblue",     linestyle=":",  label="0. damped peak", linewidth=2)
            ax.axvline(ds.time.values[idx_center],  color="dodgerblue",   linestyle="-",  label="1. center",      linewidth=2)
            ax.axvline(ds.time.values[idx_start],   color="darkgreen",    linestyle="--", label="1. start",       linewidth=2)
            ax.axvline(ds.time.values[idx_stop],    color="firebrick",    linestyle="--", label="1. stop",        linewidth=2)
            ax.axvline(ds.time.values[mean_drop],   color="deepskyblue",  linestyle="-.", label="1. w-peak",      linewidth=2)

        sig, net, sigma = _is_significant(idx_start, mean_drop)
        if sig:
            logger.debug(f"Found drop for Drifter {drifter_idx} during 1. attempt")
            resurf_str = _check_resurface_warning(mean_drop, idx_start, depth_at_drop, depth_before)
            return mean_drop, idx_start, idx_stop, resurf_str

        # --- Attempt 2: undamped score ---
        logger.info(f"Drifter {drifter_idx}: attempt 1 not significant, trying undamped score")
        idx_center    = np.nanargmax(combined)
        idx_start     = _hump_peak(score_bwd_mod, idx_center)
        idx_stop      = _hump_peak(score_fwd_mod, idx_center)
        mean_drop     = _refine_to_w_peak(idx_center)
        idx_start     = _refine_to_depth_min(idx_start)
        idx_stop      = _refine_to_depth_max(idx_stop)
        depth_at_drop = depth_smooth.values[mean_drop]
        depth_before  = depth_smooth.values[idx_start]

        if logger.isEnabledFor(logging.DEBUG):
            ax.axvline(ds.time.values[idx_center], color="deepskyblue", linestyle="-",  label="2. center", linewidth=2)
            ax.axvline(ds.time.values[idx_start],  color="limegreen",  linestyle="--", label="2. start",  linewidth=2)
            ax.axvline(ds.time.values[idx_stop],   color="salmon", linestyle="--", label="2. stop",  linewidth=2)

        sig, net, sigma = _is_significant(idx_start, mean_drop)
        if sig:
            logger.debug(f"Found drop for Drifter {drifter_idx} during 2. attempt")
            resurf_str = _check_resurface_warning(mean_drop, idx_start, depth_at_drop, depth_before)
            return mean_drop, idx_start, idx_stop, resurf_str

        # --- Attempt 3: second-best damped peak ---
        logger.info(f"Drifter {drifter_idx}: attempt 2 not significant, trying second-best damped peak")
        mid_damped, _  = _next_best_peak(combined_damped, mid_damped)
        idx_center     = _hump_peak(combined, mid_damped)
        idx_start      = _hump_peak(score_bwd_mod, idx_center)
        idx_stop       = _hump_peak(score_fwd_mod, idx_center)
        mean_drop      = _refine_to_w_peak(idx_center)
        idx_start      = _refine_to_depth_min(idx_start)
        idx_stop       = _refine_to_depth_max(idx_stop)
        depth_at_drop  = depth_smooth.values[mean_drop]
        depth_before   = depth_smooth.values[idx_start]

        if logger.isEnabledFor(logging.DEBUG):
            ax.axvline(ds.time.values[idx_center], color="blue",      linestyle="-",  label="3. center", linewidth=2)
            ax.axvline(ds.time.values[idx_start],  color="mediumseagreen", linestyle="--", label="3. start",  linewidth=2)
            ax.axvline(ds.time.values[idx_stop],  color="peru", linestyle="--", label="3. stop",  linewidth=2)
            fig.legend()

        sig, net, sigma = _is_significant(idx_start, mean_drop)
        if sig:
            logger.debug(f"Found drop for Drifter {drifter_idx} during 3. attempt")
            resurf_str = _check_resurface_warning(mean_drop, idx_start, depth_at_drop, depth_before)
            return mean_drop, idx_start, idx_stop, resurf_str

        # --- All attempts failed ---
        logger.info(f"Drifter {drifter_idx}: no significant drop found — particle likely already deep enough")
        if logger.isEnabledFor(logging.DEBUG):
            ax.set_title(f"Drifter {drifter_idx} — no significant drop", color="red")

        return None, None, None, None

    finally:
        if logger.isEnabledFor(logging.DEBUG):
            fig.legend()


def _gaussian_weight(score, center_idx, sigma_days, steps_per_day):
    """Gaussian damping weight centred on center_idx."""
    t            = np.arange(len(score))
    dist_days    = np.abs(t - center_idx) / steps_per_day
    return np.exp(-(dist_days ** 2) / (2 * sigma_days ** 2))

def debug_plot(drifter_idx, ds, pw_time_idx, limit=None):

    fig, ax = plt.subplots(1, figsize=(16,5))
    ax.set_title(f"Drifter {drifter_idx}")
    
    ax.plot(ds.time, ds.depth.isel(drifter=drifter_idx), color="black", linewidth=2)
    T = 6
    if 'temp_collumn' in ds.data_vars:
        sc = ax.contourf(ds.time, ds.depth_levels, ds.temp_column.isel(drifter=drifter_idx), levels = np.arange(-2,7), cmap = cmo.cm.thermal, extend="max")
        label = "Water column temperature (°C)"
    else:
        levels = np.arange(-2,7)  # adjust to your data range and step size
        cmap = cmo.cm.thermal  # or any other cmocean colormap
        norm = mpl.colors.BoundaryNorm(levels, ncolors=cmap.N)
        sc = ax.scatter(ds.time.values[::T], (np.ones(ds.time.shape)*2)[::T], c=ds.temp_surface.isel(drifter=drifter_idx).values[::T], cmap=cmap, norm=norm, marker="s")
        label = "Surface layer temperature (°C)"

    ax.set_ylim(np.nanmax(ds.depth.isel(drifter=drifter_idx)), 0)
    ax.set_ylabel("Depth (m)")
    fig.colorbar(sc, extend='max', label=label)
    ax.margins(x=0)
    ax.axvline(x=ds.time.values[pw_time_idx], color="black", linestyle="--", label="PW or MW encounter", linewidth=3)
    ax.legend(loc="lower right")

    if limit is not None:
        ax.axhline(y=limit)

    return fig, ax


def _extreme(vals):
    """Return the value with largest absolute magnitude."""
    v_max = np.nanmax(vals)
    v_min = np.nanmin(vals)
    return v_min if abs(v_min) > abs(v_max) else v_max

def _props_at_point(drifter_idx, t, ds, STEPS_PER_DAY):
    """Extract particle properties at a single timestep."""

    Ro = ds.Ro.values[drifter_idx, t-STEPS_PER_DAY:t+STEPS_PER_DAY] if 'Ro' in ds else None
    strain = ds.strain.values[drifter_idx, t-STEPS_PER_DAY:t+STEPS_PER_DAY] if 'strain' in ds else None
    ow = ds.okubo_weiss.values[drifter_idx, t-STEPS_PER_DAY:t+STEPS_PER_DAY] if 'okubo_weiss' in ds else None
    w = ds.w.values[drifter_idx, t-STEPS_PER_DAY:t+STEPS_PER_DAY] if 'w' in ds else None
    return {
        'lat':   ds.lat.values[drifter_idx, t],
        'lon':   ds.lon.values[drifter_idx, t],
        'depth': ds.depth.values[drifter_idx, t],
        'Ro_mean':   np.nanmean(abs(Ro)) if Ro is not None else np.nan,
        'vorticity': ds.vort.values[drifter_idx, t],
        'Ro_90':  np.nanpercentile(abs(Ro), 90) if Ro is not None else np.nan,
        'Ro_max': _extreme(Ro) if Ro is not None else np.nan,
        'strain_90': np.nanpercentile(abs(strain), 90) if strain is not None else np.nan,
        'strain_max': _extreme(strain) if strain is not None else np.nan,
        'okubo_weiss_mean': np.nanmean(ow) if ow is not None else np.nan,
        'okubo_weiss_max': _extreme(ow) if ow is not None else np.nan,
        'w_mean': np.nanmean(w) if w is not None else np.nan,
        'w_max': np.nanmax(w) if w is not None else np.nan,

    }

def _props_in_window(drifter_idx, t_start, t_end, ds, STEPS_PER_DAY):
    """Extract aggregated particle properties over a time window [t_start, t_end]."""
    t_start = max(0, t_start)
    t_end   = min(ds.time.shape[0], t_end)

    # Window too small — fall back to daily window around midpoint
    if t_end - t_start < STEPS_PER_DAY:
        mid     = (t_start + t_end) // 2
        t_start = max(0, mid - STEPS_PER_DAY // 2)
        t_end   = min(ds.time.shape[0], mid + STEPS_PER_DAY // 2)
        print(f"Drifter {drifter_idx}: window too small, using daily window [{t_start}, {t_end}]")

    def _get(var):
        if var in ds:
            return ds[var].values[drifter_idx, t_start:t_end]
        return None

    vort   = _get('vort')
    Ro     = _get('Ro')
    strain = _get('strain')
    ow     = _get('okubo_weiss')
    N2     = _get('N2')
    Ri     = _get('Ri')
    depth  = _get('depth')

    props = {
        'mean_vorticity':   np.nanmean(vort),
        'Ro_mean':          np.nanmean(abs(Ro)),
        'Ro_90':            np.nanpercentile(abs(Ro), 90),
        'Ro_max':           _extreme(Ro),
        'strain_mean':      np.nanmean(abs(strain)),
        'strain_max':       _extreme(strain),
        'strain_90':        np.nanpercentile(abs(strain), 90),
        'okubo_weiss_mean': np.nanmean(ow),
        'okubo_weiss_max':  _extreme(ow),
        'depth_mean':       np.nanmean(depth),
        'depth_min':        np.nanmin(depth),
        'depth_max':        np.nanmax(depth),
    }

    if N2 is not None:
        props['N2_mean'] = np.nanmean(N2)
    if Ri is not None:
        props['Ri_mean'] = np.nanmean(Ri)
        props['Ri_min']  = np.nanmin(Ri)

    return props

def extract_properties(drifter_idx, t0, ds, STEPS_PER_DAY, start_idx=None, drop_idx=None, stop_idx=None, category='subduction', resurf_warning_str=""):
    """Build a result record for one drifter.
    
    If start_idx/drop_idx/stop_idx are None, drop/start/descent fields are filled with NaN.
    """
    record = {
        'drifter_idx':    drifter_idx,
        'category':       category,
        'warning':        resurf_warning_str,
        'time':           t0,
        'time_drop':      int(drop_idx)  if drop_idx  is not None else np.nan,
        'time_start':     int(start_idx) if start_idx is not None else np.nan,
        'time_stop':      int(stop_idx)  if stop_idx  is not None else np.nan,
        'datetime':       ds.time.values[t0],
        'datetime_drop':  ds.time.values[int(drop_idx)]  if drop_idx  is not None else pd.NaT,
        'datetime_start': ds.time.values[int(start_idx)] if start_idx is not None else pd.NaT,
        'datetime_stop' : ds.time.values[int(stop_idx)]  if stop_idx  is not None else pd.NaT,

    }

    # Point properties at subduction confirmation and at drop center
    for k, v in _props_at_point(drifter_idx, t0,       ds, STEPS_PER_DAY).items():
        record[k] = v
    #for k, v in _props_at_point(drifter_idx, drop_idx, ds).items():
    #    record[f'{k}_drop'] = v
    #for k, v in _props_at_point(drifter_idx, start_idx, ds).items():
    #    record[f'{k}_start'] = v

    # Window properties: start → drop center (the actual descent)
    #for k, v in _props_in_window(drifter_idx, start_idx, drop_idx, ds, STEPS_PER_DAY).items():
    #    record[f'{k}_descent'] = v

    # Drop/start point properties
    for k, v in (_props_at_point(drifter_idx, int(drop_idx),  ds, STEPS_PER_DAY) if drop_idx  is not None else {}).items():
        record[f'{k}_drop']  = v
    for k, v in (_props_at_point(drifter_idx, int(start_idx), ds, STEPS_PER_DAY) if start_idx is not None else {}).items():
        record[f'{k}_start'] = v
    for k, v in (_props_at_point(drifter_idx, int(stop_idx), ds, STEPS_PER_DAY) if stop_idx is not None else {}).items():
        record[f'{k}_stop'] = v

    # Window/descent properties
    for k, v in (_props_in_window(drifter_idx, int(start_idx), int(drop_idx), ds, STEPS_PER_DAY) if (start_idx is not None and drop_idx is not None) else {}).items():
        record[f'{k}_descent'] = v
    

    return record

