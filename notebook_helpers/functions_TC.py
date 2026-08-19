import pandas as pd
import numpy as np
import numpy
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from argovisHelpers import helpers as avh
from argovisHelpers import analysis as ava
import math
import gsw
from dateutil import parser
from scipy.interpolate import interp1d
from math import radians, sin, cos, sqrt, atan2
from notebook_helpers.functions import plot_maps
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="scipy")
warnings.filterwarnings("ignore", message="Internal error in root finding")

# ── Dataset key mappings ──────────────────────────────────────────────────────
TEMPERATURE_KEY = {
    'argo':      'temperature',
    'cchdo':     'temperature',
    'easyocean': 'ctd_temperature',
}
SALINITY_KEY = {
    'argo':      'salinity',
    'cchdo':     'salinity',
    'easyocean': 'ctd_salinity',
}


# ─────────────────────────────────────────────────────────────────────────────
# Utility helpers
# ─────────────────────────────────────────────────────────────────────────────
        
def radius_to_box(center, radius_km):
    """Convert a center + radius (km) to a bounding box [[lon_min, lat_min], [lon_max, lat_max]]."""
    lon, lat = center
    dlat = radius_km / 111.0
    dlon = radius_km / (111.0 * math.cos(math.radians(lat)))
    return [[lon - dlon, lat - dlat], [lon + dlon, lat + dlat]]


def haversine_km(lon1, lat1, lon2, lat2):
    """Great-circle distance in km between two (lon, lat) points."""
    R = 6371.0
    phi1, phi2 = radians(lat1), radians(lat2)
    dphi    = radians(lat2 - lat1)
    dlambda = radians(lon2 - lon1)
    a = sin(dphi / 2)**2 + cos(phi1) * cos(phi2) * sin(dlambda / 2)**2
    return 2 * R * atan2(sqrt(a), sqrt(1 - a))


def _arr(val):
    """Safely convert Profile.getvar() output to a plain float ndarray."""
    return val.filled(np.nan) if hasattr(val, 'filled') else np.asarray(val, dtype=float)


def infer_dataset(p):
    """Infer dataset key ('argo', 'cchdo', 'easyocean') from rawdata or variable names."""
    try:
        source_str = p._rawdata['source'][0]['source'][0]
        if 'argo'      in source_str: return 'argo'
        if 'cchdo'     in source_str: return 'cchdo'
        if 'easyocean' in source_str: return 'easyocean'
    except Exception:
        pass
    names = p.variable_names() if callable(p.variable_names) else p.variable_names
    if 'ctd_temperature' in names: return 'easyocean'
    if 'temperature'     in names: return 'argo'
    return None


def get_flat_profiles(obj):
    """Recursively flatten nested lists/dicts of Profile objects into a flat list."""
    found = []
    if "Profile" in str(type(obj)): return [obj]
    if isinstance(obj, list):
        for item in obj: found.extend(get_flat_profiles(item))
    elif isinstance(obj, dict):
        for val in obj.values(): found.extend(get_flat_profiles(val))
    return found


def interpolate_to_grid(p, varname, grid):
    """Interpolate a profile variable onto a common pressure grid."""
    pres = p.getvar('pressure')
    vals = p.getvar(varname)
    if pres is None or vals is None:
        return np.full(len(grid), np.nan)
    pres  = np.array(pres, dtype=float)
    vals  = np.array(vals, dtype=float)
    valid = ~np.isnan(pres) & ~np.isnan(vals)
    if valid.sum() < 2:
        return np.full(len(grid), np.nan)
    order = np.argsort(pres[valid])
    f = interp1d(pres[valid][order], vals[valid][order],
                 bounds_error=False, fill_value=np.nan)
    return f(np.array(grid, dtype=float))


def _scalar_mld(p, key='MLD_potdensity'):
    """Extract a scalar MLD value from a profile."""
    v = p.getvar(key)
    return float(v[0]) if v is not None and len(v) > 0 else np.nan


def extract_profile_data(p, grid):
    """
    Extract T, S, MLD, and potential density from a Profile object,
    interpolating variables onto the common pressure grid.
    Uses interpolate_to_grid for each variable individually since
    derived variables (SA, PD) are added via setvar and are not
    visible to ava.interpolate_all.
    """
    mld_potd = _scalar_mld(p, 'MLD_potdensity')
    mld_temp = _scalar_mld(p, 'MLD_temperature')
    mld_dens = _scalar_mld(p, 'MLD_density')

    T = interpolate_to_grid(p, 'temperature', grid)
    if np.all(np.isnan(T)):
        T = interpolate_to_grid(p, 'ctd_temperature', grid)

    return {
        'profile_id':        p.id,
        'timestamp':         p.timestamp,
        'longitude':         p.longitude,
        'latitude':          p.latitude,
        'pressure':          np.array(grid, dtype=float),
        'temperature':       T,
        'salinity':          interpolate_to_grid(p, 'absolute_salinity',  grid),
        'potential_density': interpolate_to_grid(p, 'potential_density',  grid),
        'MLD':               mld_potd,
        'MLD_temperature':   mld_temp,
        'MLD_density':       mld_dens,
        'MLD_potdensity':    mld_potd,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Track data
# ─────────────────────────────────────────────────────────────────────────────

def get_yearly_storm_tracks(year, APIKEY):
    """
    Query TC track data year-by-year in monthly chunks.
    Returns a dictionary keyed by storm ID. Each value has:
      'name'    : storm name (str)
      'id'      : storm ID (str)
      'regions' : list of track points, each with 'center', 'time', 'wind_class'
    'wind_class' is the IBTrACS storm classification string (e.g. 'TD', 'TS',
    'HU', 'TY', 'ST', 'EX') — wind speed is not available from this endpoint.
    """
    temp_grouped = {}

    for month in range(1, 13):
        start_date = f"{year}-{month:02d}-01T00:00:00Z"
        end_year   = year if month < 12 else year + 1
        end_month  = month + 1 if month < 12 else 1
        end_date   = f"{end_year}-{end_month:02d}-01T00:00:00Z"

        print(f"--- Fetching TC data: {year}-{month:02d} ---")
        try:
            raw_data = avh.query(
                'tc',
                options={'startDate': start_date, 'endDate': end_date},
                verbose=False
            )
            if not raw_data:
                continue
            for point in raw_data:
                storm_id = (point['metadata'][0] if point.get('metadata')
                            else point['_id'].split('_')[0])

                # ── Storm name — try common field names in the track point ────
                # The canonical name lives in the TC metadata document; the
                # track point may not carry it in all API versions. We fall
                # back to a post-loop metadata query for any still-UNNAMED
                # storms below.
                name = ''
                for _key in ('name', 'NAME', 'storm_name'):
                    _val = point.get(_key, '')
                    if _val and str(_val).strip().upper() not in ('', 'NOT_NAMED', 'UNNAMED'):
                        name = str(_val).strip().upper()
                        break

                # ── Storm class — ordinal intensity proxy ────────────────────
                # Wind speed is not available from the Argovis TC endpoint.
                # 'class' carries the IBTrACS classification string instead
                # (TD, TS, HU, TY, TC, ST, EX, …).
                wind_class = point.get('class', '')

                if storm_id not in temp_grouped:
                    temp_grouped[storm_id] = {
                        'name':    name if name else 'UNNAMED',
                        'regions': [],
                        'id':      storm_id,
                    }
                elif name and temp_grouped[storm_id]['name'] == 'UNNAMED':
                    temp_grouped[storm_id]['name'] = name

                if not any(p['time'] == point['timestamp']
                           for p in temp_grouped[storm_id]['regions']):
                    temp_grouped[storm_id]['regions'].append({
                        'center':     point['geolocation']['coordinates'],
                        'time':       point['timestamp'],
                        'wind_class': wind_class,
                    })
        except Exception as e:
            print(f"Skipping {year}-{month:02d} due to error: {e}")

    # ── Post-process: try the tc/meta endpoint for any still-UNNAMED storms ──
    unnamed_ids = [sid for sid in temp_grouped
                   if temp_grouped[sid]['name'] == 'UNNAMED']
    if unnamed_ids:
        print(f"Fetching names for {len(unnamed_ids)} unnamed storm(s)…")
        for sid in unnamed_ids:
            try:
                meta = avh.query('tc/meta', options={'id': sid}, verbose=False)
                if meta and isinstance(meta, list) and len(meta) > 0:
                    mname = meta[0].get('name', '')
                    if mname and str(mname).strip().upper() not in ('', 'NOT_NAMED'):
                        temp_grouped[sid]['name'] = str(mname).strip().upper()
            except Exception:
                pass

    for sid in temp_grouped:
        temp_grouped[sid]['regions'].sort(key=lambda x: x['time'])

    print(f"Successfully organized {len(temp_grouped)} tracks.")
    return temp_grouped
# ─────────────────────────────────────────────────────────────────────────────
# Derived oceanographic properties
# ─────────────────────────────────────────────────────────────────────────────

def compute_derived_properties(p, temperature_key=None, salinity_key=None, ds=None):
    if temperature_key is None: temperature_key = TEMPERATURE_KEY
    if salinity_key    is None: salinity_key    = SALINITY_KEY
    if ds is None:              ds              = infer_dataset(p)

    SA = gsw.conversions.SA_from_SP(
        p.getvar(salinity_key[ds]), p.getvar('pressure'),
        p.longitude, p.latitude)
    p.setvar('absolute_salinity', SA)

    CT = gsw.conversions.CT_from_t(
        SA, p.getvar(temperature_key[ds]), p.getvar('pressure'))
    p.setvar('conservative_temperature', CT)

    sigma0 = gsw.sigma0(SA, CT) + 1000
    p.setvar('potential_density', sigma0)

    rho_insitu = gsw.rho(SA, CT, p.getvar('pressure'))
    p.setvar('density', rho_insitu)

    pressure    = p.getvar('pressure')
    temperature = p.getvar(temperature_key[ds])
    density     = p.getvar('density')
    pot_density = p.getvar('potential_density')

    # ── MLD estimates — isolated so a failure NaNs only the MLD, not the profile ──
    # ava.MLD_estimate fits a PCHIP spline and calls pchip.solve() to find where the
    # profile crosses (value at reference_pressure) + threshold_delta. Two failure modes:
    #   1. The profile never samples reference_pressure (shallowest level deeper than it),
    #      so an MLD referenced to that depth is undefined -> return None.
    #   2. scipy's _ppoly.real_roots raises "Internal error in root finding" on certain
    #      degenerate spline segments -> catch and return None rather than dropping the
    #      whole profile.
    REFERENCE_PRESSURE = 10
    _pres        = _arr(pressure)
    _finite      = _pres[np.isfinite(_pres)]
    _reaches_ref = _finite.size > 0 and np.nanmin(_finite) <= REFERENCE_PRESSURE
    if not _reaches_ref:
        _shallow = f"{np.nanmin(_finite):.1f}" if _finite.size else "n/a"
        print(f"    {p.id} | MLD skipped: shallowest level {_shallow} dbar "
              f"is below the {REFERENCE_PRESSURE}-dbar reference")

    def _safe_mld(field, thr):
        if field is None or not _reaches_ref:
            return [None]
        try:
            mld, _ = ava.MLD_estimate(pressure, field,
                                      threshold_delta=thr,
                                      reference_pressure=REFERENCE_PRESSURE)
            return [mld]
        except Exception as e:
            print(f"    {p.id} | MLD (thr={thr}) failed: {type(e).__name__}: {e}")
            return [None]

    p.setvar('MLD_temperature', _safe_mld(temperature, -0.2))
    p.setvar('MLD_density',     _safe_mld(density,      0.125))
    p.setvar('MLD_potdensity',  _safe_mld(pot_density,  0.03))
    p.setvar('MLD', p.getvar('MLD_potdensity'))

    # ── N² and oscillation period — isolated so failures don't affect MLD/PD ──
    try:
        SA_arr   = _arr(SA)
        CT_arr   = _arr(CT)
        pres_arr = _arr(pressure)
        valid    = ~np.isnan(SA_arr) & ~np.isnan(CT_arr) & ~np.isnan(pres_arr)
        SA_v, CT_v, pres_v = SA_arr[valid], CT_arr[valid], pres_arr[valid]
        order = np.argsort(pres_v)
        N2, p_mid = gsw.Nsquared(SA_v[order], CT_v[order], pres_v[order], lat=p.latitude)
        good = ~np.isnan(N2)
        p.setvar('N2', N2[good])
        p.setvar('N2_pressure', p_mid[good])
        with np.errstate(divide='ignore', invalid='ignore'):
            period = np.where(N2[good] > 0, (2 * np.pi / np.sqrt(N2[good])) / 60, np.nan)
        p.setvar('oscillation_period', period)
    except Exception as e:
        p.setvar('N2', None)
        p.setvar('N2_pressure', None)
        p.setvar('oscillation_period', None)
        print(f"    {p.id} | N2 failed: {type(e).__name__}: {e}")

# ─────────────────────────────────────────────────────────────────────────────
# Deduplication
# ─────────────────────────────────────────────────────────────────────────────

def deduplicate_storm_results(storm_results, track_points,
                               max_pair_dist_km=None, radius_km=None):
    """
    Remove profiles that appear in multiple overlapping track-point windows,
    and remove spatially redundant profiles within each window.

    Parameters
    ----------
    max_pair_dist_km : float or None
        Within each step, profiles of the SAME timing type (both before or
        both after TC passage) that are closer than this are considered
        redundant — keep only the one nearest the track point.
        Defaults to radius_km.
        Before/after pairs are NEVER removed by this criterion.
    radius_km : float or None
        Hard spatial cutoff applied after assignment — profiles beyond this
        distance from their assigned track point are removed.
    """
    if max_pair_dist_km is None:
        max_pair_dist_km = radius_km

    # ── Pass 1: assign each profile to its nearest track point ───────────────
    best = {}  # profile_id -> (dist_km, dt_s, step_idx)
    for step_idx, step_data in enumerate(storm_results):
        if not step_data:
            continue
        pt = track_points[step_idx]
        tc_lon, tc_lat = pt['center']
        tc_time = parser.parse(pt['time']).replace(tzinfo=None)
        for p in get_flat_profiles(step_data):
            try:
                key  = p.id
                dist = haversine_km(tc_lon, tc_lat, p.longitude, p.latitude)
                dt   = abs((pd.to_datetime(p.timestamp)
                              .to_pydatetime()
                              .replace(tzinfo=None) - tc_time).total_seconds())
                if key not in best or (dist, dt) < (best[key][0], best[key][1]):
                    best[key] = (dist, dt, step_idx)
            except Exception:
                continue

    # ── Pass 2: rebuild flat lists keeping only the best-assigned copy ────────
    deduplicated = [[] for _ in storm_results]
    for step_idx, step_data in enumerate(storm_results):
        if not step_data:
            continue
        for p in get_flat_profiles(step_data):
            try:
                if best.get(p.id, (None, None, -1))[2] == step_idx:
                    deduplicated[step_idx].append(p)
            except Exception:
                continue

    # ── Pass 2.5: remove profiles beyond true radius (box approximation fix) ──
    if radius_km is not None:
        for step_idx, profiles in enumerate(deduplicated):
            if not profiles:
                continue
            tc_lon, tc_lat = track_points[step_idx]['center']
            deduplicated[step_idx] = [
                p for p in profiles
                if haversine_km(tc_lon, tc_lat, p.longitude, p.latitude) <= radius_km
            ]

    # ── Pass 3: spatial deduplication within each step ───────────────────────
    # Runs separately on before and after groups so that a before/after pair
    # that happens to be spatially close is never removed.
    if max_pair_dist_km is not None:
        for step_idx, profiles in enumerate(deduplicated):
            if not profiles:
                continue
            tc_lon, tc_lat = track_points[step_idx]['center']
            tc_time = parser.parse(track_points[step_idx]['time']).replace(tzinfo=None)

            before_profiles, after_profiles = [], []
            for p in profiles:
                try:
                    t = pd.to_datetime(p.timestamp).to_pydatetime().replace(tzinfo=None)
                    (before_profiles if t < tc_time else after_profiles).append(p)
                except Exception:
                    before_profiles.append(p)

            def _spatial_dedup(group):
                """Keep one profile per spatial cluster (closest to track point wins)."""
                sorted_group = sorted(
                    group,
                    key=lambda p: haversine_km(tc_lon, tc_lat, p.longitude, p.latitude)
                )
                kept = []
                for candidate in sorted_group:
                    too_close = any(
                        haversine_km(candidate.longitude, candidate.latitude,
                                     k.longitude, k.latitude) < max_pair_dist_km
                        for k in kept
                    )
                    if not too_close:
                        kept.append(candidate)
                return kept

            deduplicated[step_idx] = (
                _spatial_dedup(before_profiles) + _spatial_dedup(after_profiles)
            )

    return deduplicated


# ─────────────────────────────────────────────────────────────────────────────
# Pairing and delta computation
# ─────────────────────────────────────────────────────────────────────────────

def build_pair_dataframe(all_insitu_deduplicated, all_result_track_pts,
                         interpolation_levels, max_pair_dist_km=None, radius_km=None):
    """
    Globally match before/after profile pairs across all track points per storm.
    Each profile appears in at most one pair.

    Pairing criteria:
      1. Both profiles must be within radius_km of the track point.
      2. Distance between the two profiles must be <= max_pair_dist_km
         (defaults to radius_km; can be set independently).
      3. Among valid candidates, keep the pair with minimum spatial distance
         between the two profiles.
      4. Tiebreak: minimum time gap between the two profiles (|t_after - t_before|).

    Parameters
    ----------
    max_pair_dist_km : float or None
        Maximum allowed distance between the before and after profile.
        Independent of radius_km:
          - radius_km controls how far each profile can be from the track point.
          - max_pair_dist_km controls how far the two profiles can be from each other.
        Defaults to radius_km if not set.

    Returns
    -------
    all_paired : nested dict
        {
          storm_id: {
            step_idx: {
              'track_lon', 'track_lat', 'track_time',
              'pair_dist_km', 'pair_dt_hours',
              dataset_name: {
                'before':             profile data dict,
                'after':              profile data dict,
                'after_minus_before': delta dict,
              }
            }
          }
        }
        Each profile data dict contains:
          profile_id, timestamp, longitude, latitude,
          pressure (common grid), temperature, salinity,
          potential_density, MLD, MLD_temperature, MLD_density, MLD_potdensity.
        The delta dict contains:
          pressure, temperature, salinity, potential_density (arrays), MLD (scalar).

    df_pairs : pd.DataFrame
        Flat table of scalar values for mapping. One row per assigned pair with:
          storm_id, step_idx, track_lon, track_lat, track_time, dataset,
          profile_before_id, profile_after_id, time_before, time_after,
          pair_dist_km, pair_dt_hours,
          MLD_before, MLD_after, delta_MLD.
    """
    if max_pair_dist_km is None:
        max_pair_dist_km = radius_km

    grid       = np.array(interpolation_levels, dtype=float)
    all_paired = {}
    rows       = []

    for storm_id, steps in all_insitu_deduplicated.items():
        track_pts = all_result_track_pts[storm_id]
        all_paired[storm_id] = {}

        # ── Collect all before/after profiles across all steps ────────────────
        all_before, all_after = [], []
        for step, pt in zip(steps, track_pts):
            tc_time = parser.parse(pt['time']).replace(tzinfo=None)
            for p in step:
                try:
                    t = pd.to_datetime(p.timestamp).to_pydatetime().replace(tzinfo=None)
                    (all_before if t < tc_time else all_after).append((p, t))
                except Exception:
                    continue

        if not all_before or not all_after:
            print(f"  {storm_id}: no before/after profiles — skipping")
            continue

        print(f"  {storm_id}: {len(all_before)} before, {len(all_after)} after profiles")

        # ── Generate all valid (track_point, before, after) candidates ────────
        candidates = []
        for step_idx, (step, pt) in enumerate(zip(steps, track_pts)):
            if not step:
                continue
            tc_lon, tc_lat = pt['center']

            before_near = [
                (p, t) for p, t in all_before
                if haversine_km(tc_lon, tc_lat, p.longitude, p.latitude) <= radius_km
            ]
            after_near = [
                (p, t) for p, t in all_after
                if haversine_km(tc_lon, tc_lat, p.longitude, p.latitude) <= radius_km
            ]

            for p_b, t_b in before_near:
                for p_a, t_a in after_near:
                    dist = haversine_km(p_b.longitude, p_b.latitude,
                                        p_a.longitude, p_a.latitude)
                    if max_pair_dist_km is not None and dist > max_pair_dist_km:
                        continue
                    dt = abs((t_a - t_b).total_seconds())
                    candidates.append((dist, dt, step_idx, p_b, t_b, p_a, t_a))

        # ── Greedy assignment: best score first, each entity used once ────────
        candidates.sort(key=lambda x: (x[0], x[1]))  # sort by (dist, dt)

        used_steps  = set()
        used_before = set()
        used_after  = set()
        assigned    = []

        for dist, dt, step_idx, p_b, t_b, p_a, t_a in candidates:
            if step_idx in used_steps:  continue
            if p_b.id   in used_before: continue
            if p_a.id   in used_after:  continue
            used_steps.add(step_idx)
            used_before.add(p_b.id)
            used_after.add(p_a.id)
            assigned.append((dist, dt, step_idx, p_b, t_b, p_a, t_a))

        print(f"  {storm_id}: {len(assigned)} pairs assigned")

        # ── Build nested output and flat rows ─────────────────────────────────
        for dist, dt, step_idx, p_b, t_b, p_a, t_a in assigned:
            pt             = track_pts[step_idx]
            tc_lon, tc_lat = pt['center']
            tc_time        = parser.parse(pt['time']).replace(tzinfo=None)
            ds             = infer_dataset(p_b)

            before_data = extract_profile_data(p_b, grid)
            after_data  = extract_profile_data(p_a, grid)

            mld_b = before_data['MLD']
            mld_a = after_data['MLD']

            delta_data = {
                'pressure':          grid.copy(),
                'temperature':       after_data['temperature']       - before_data['temperature'],
                'salinity':          after_data['salinity']          - before_data['salinity'],
                'potential_density': after_data['potential_density'] - before_data['potential_density'],
                'MLD': mld_a - mld_b
                       if not (np.isnan(mld_b) or np.isnan(mld_a)) else np.nan,
            }

            all_paired[storm_id][step_idx] = {
                'track_lon':     tc_lon,
                'track_lat':     tc_lat,
                'track_time':    tc_time,
                'pair_dist_km':  round(dist, 2),
                'pair_dt_hours': round(dt / 3600, 2),
                'pair_dt_days':  round(dt / 86400, 2),
                ds: {
                    'before':             before_data,
                    'after':              after_data,
                    'after_minus_before': delta_data,
                }
            }

            rows.append({
                'storm_id':          storm_id,
                'step_idx':          step_idx,
                'track_lon':         tc_lon,
                'track_lat':         tc_lat,
                'track_time':        tc_time,
                'dataset':           ds,
                'profile_before_id': p_b.id,
                'profile_after_id':  p_a.id,
                'time_before':       t_b,
                'time_after':        t_a,
                'pair_dist_km':      round(dist, 2),
                'pair_dt_hours':     round(dt / 3600, 2),
                'pair_dt_days':      round(dt / 86400, 2),
                'MLD_before':        mld_b,
                'MLD_after':         mld_a,
                'delta_MLD':         delta_data['MLD'],
            })

    df_pairs = pd.DataFrame(rows) if rows else pd.DataFrame()
    if not df_pairs.empty:
        print(f"\nBuilt dataframe: {len(df_pairs)} pairs across "
              f"{df_pairs['storm_id'].nunique()} storms")
    else:
        print("No pairs found.")

    return all_paired, df_pairs


# ─────────────────────────────────────────────────────────────────────────────
# Visualization
# ─────────────────────────────────────────────────────────────────────────────

def plot_storm_insitu(storm_label, track_points, all_profiles_raw,
                      radius_km, delta_days,
                      all_paired_storm=None,
                      time_linestyles=None, depth_ylim=500):
    """
    Two maps per storm in a single figure:
      1. All profiles returned by the API (pre-deduplication) on a stock imagery
         background (yellow markers). The TC track is colored following NHC/JTWC
         conventions by storm class (Super Typhoon → Hurricane/Typhoon →
         Tropical Storm → Depression → Extratropical → Disturbance).
      2. Only the selected before/after pairs, colored by days relative to TC
         passage (RdBu colormap; blue = before, red = after). Markers carry a
         black edge so light-shaded profiles near zero time offset remain visible.
    """
    if time_linestyles is None:
        time_linestyles = ['-', '--', ':', '-.']

    from matplotlib.collections import LineCollection
    from matplotlib.patches import Patch

    # ── NHC/JTWC color conventions ────────────────────────────────────────────
    CLASS_COLOR = {
        'ST': '#CC00CC',  # Super Typhoon        — magenta
        'HU': '#CC0000',  # Hurricane            — red
        'TY': '#CC0000',  # Typhoon              — red
        'TC': '#CC0000',  # Tropical Cyclone     — red
        'TS': '#FFDD00',  # Tropical Storm       — yellow
        'SS': '#FFA500',  # Subtropical Storm    — orange
        'TD': '#00AA00',  # Tropical Depression  — green
        'SD': '#0066FF',  # Subtropical Dep.     — blue
        'EX': '#777777',  # Extratropical        — gray
        'ET': '#777777',
        '':   '#333333',  # Disturbance/unknown  — dark gray
        'DB': '#333333',
        'LO': '#333333',
        'WV': '#333333',
        'XX': '#333333',
    }
    CLASS_NAME = {
        'ST': 'Super Typhoon',
        'HU': 'Hurricane',
        'TY': 'Typhoon',
        'TC': 'Tropical Cyclone',
        'TS': 'Tropical Storm',
        'SS': 'Subtropical Storm',
        'TD': 'Tropical Depression',
        'SD': 'Subtropical Depression',
        'EX': 'Extratropical',
        'ET': 'Extratropical',
        '':   'Disturbance',
        'DB': 'Disturbance',
        'LO': 'Low',
        'WV': 'Wave',
        'XX': 'Unknown',
    }
    # Display order in legend (strongest → weakest)
    CLASS_LEGEND_ORDER = ['ST', 'HU', 'TY', 'TC', 'TS', 'SS', 'TD', 'SD',
                          'EX', 'ET', '', 'DB', 'LO', 'WV', 'XX']

    # ── Collect ALL raw profiles for Map 1 ───────────────────────────────────
    all_raw_lons, all_raw_lats = [], []
    for step_data in all_profiles_raw:
        for p in get_flat_profiles(step_data):
            try:
                all_raw_lons.append(p.longitude)
                all_raw_lats.append(p.latitude)
            except Exception:
                continue

    # ── Collect SELECTED PAIR profiles for Map 2 ─────────────────────────────
    pair_lons, pair_lats, pair_time_deltas = [], [], []
    if all_paired_storm:
        for step_idx, entry in all_paired_storm.items():
            ds_keys = [k for k in entry if k not in (
                'track_lon', 'track_lat', 'track_time',
                'pair_dist_km', 'pair_dt_hours', 'pair_dt_days')]
            if not ds_keys:
                continue
            ds = ds_keys[0]
            tc_time = entry['track_time']
            if hasattr(tc_time, 'tzinfo') and tc_time.tzinfo is not None:
                tc_time = tc_time.replace(tzinfo=None)
            for timing in ('before', 'after'):
                pdata = entry[ds].get(timing)
                if pdata is None:
                    continue
                try:
                    d_time  = pd.to_datetime(pdata['timestamp'])\
                                .to_pydatetime().replace(tzinfo=None)
                    delta_t = (d_time - tc_time).total_seconds() / 86400
                    pair_lons.append(pdata['longitude'])
                    pair_lats.append(pdata['latitude'])
                    pair_time_deltas.append(delta_t)
                except Exception:
                    continue

    if not all_raw_lons and not pair_lons:
        print(f"{storm_label}: no colocated profiles found.")
        return

    track_lons_flat = [pt['center'][0] for pt in track_points]
    track_lats_flat = [pt['center'][1] for pt in track_points]
    all_lons_flat   = track_lons_flat + all_raw_lons + pair_lons
    all_lats_flat   = track_lats_flat + all_raw_lats + pair_lats

    pad = max(radius_km / 111, 2.0)
    extent = [
        min(all_lons_flat) - pad,
        max(all_lons_flat) + pad,
        max(min(all_lats_flat) - pad, -90),
        min(max(all_lats_flat) + pad,  90),
    ]
    lon_span = extent[1] - extent[0]
    lat_span = extent[3] - extent[2]
    fig_w    = 12
    fig_h    = max(4, fig_w * lat_span / lon_span)

    fig, (ax1, ax2) = plt.subplots(
        2, 1,
        figsize=(fig_w, fig_h * 2),
        subplot_kw={'projection': ccrs.PlateCarree()},
        constrained_layout=True,
    )

    # ── Map 1: all raw profiles + NHC-style colored track ────────────────────
    ax1.stock_img()
    ax1.coastlines(resolution='50m', color='black', linewidth=0.5)
    ax1.set_extent(extent, crs=ccrs.PlateCarree())
    gl1 = ax1.gridlines(draw_labels=True, linestyle=':', alpha=0.5)
    gl1.top_labels   = False
    gl1.right_labels = False

    classes    = [pt.get('wind_class', '') for pt in track_points]
    pts_xy     = np.column_stack([track_lons_flat, track_lats_flat])
    segments   = [pts_xy[i:i+2] for i in range(len(pts_xy) - 1)]
    seg_colors = [CLASS_COLOR.get(classes[i], '#333333') for i in range(len(segments))]

    lc = LineCollection(segments, colors=seg_colors, linewidth=2.5,
                        transform=ccrs.PlateCarree(), zorder=4)
    ax1.add_collection(lc)
    ax1.scatter(track_lons_flat, track_lats_flat,
                color='black', s=15, zorder=5, transform=ccrs.PlateCarree())

    # Legend: only classes that actually appear in this track
    seen_classes  = set(classes)
    seen_labels   = set()
    legend_handles = []
    for cls in CLASS_LEGEND_ORDER:
        if cls in seen_classes:
            label = CLASS_NAME.get(cls, cls)
            if label not in seen_labels:
                legend_handles.append(
                    Patch(facecolor=CLASS_COLOR.get(cls, '#333333'), label=label)
                )
                seen_labels.add(label)

    if all_raw_lons:
        ax1.scatter(
            all_raw_lons, all_raw_lats,
            color='yellow', s=35, alpha=0.75,
            edgecolors='gray', linewidths=0.4,
            transform=ccrs.PlateCarree(), zorder=5,
        )
        legend_handles.append(
            Patch(facecolor='yellow', edgecolor='gray',
                  label=f'All returned profiles (n={len(all_raw_lons)})')
        )

    if legend_handles:
        ax1.legend(handles=legend_handles, loc='lower left', fontsize=8)

    ax1.set_title(f"{storm_label} — all returned profiles")

    # ── Map 2: selected before/after pairs only ───────────────────────────────
    ax2.set_extent(extent, crs=ccrs.PlateCarree())
    ax2.add_feature(cfeature.LAND, facecolor='lightgray')
    ax2.add_feature(cfeature.COASTLINE)
    gl2 = ax2.gridlines(draw_labels=True, linestyle=':', alpha=0.5)
    gl2.top_labels   = False
    gl2.right_labels = False

    ax2.plot(track_lons_flat, track_lats_flat,
             'k--', linewidth=2, label='Storm Track', transform=ccrs.PlateCarree())
    ax2.scatter(track_lons_flat, track_lats_flat,
                color='black', s=15, transform=ccrs.PlateCarree(), zorder=4)

    if pair_lons:
        sc = ax2.scatter(
            pair_lons, pair_lats,
            c=pair_time_deltas, cmap='RdBu',
            vmin=-delta_days, vmax=delta_days,
            alpha=0.85, s=60,
            edgecolors='black', linewidths=0.6,
            transform=ccrs.PlateCarree(), zorder=5,
        )
        cbar = plt.colorbar(sc, ax=ax2, orientation='horizontal', pad=0.05, shrink=0.7)
        cbar.set_label('Days relative to storm passage  (− before, + after)')
        n_pairs = len(pair_lons) // 2
        ax2.set_title(
            f"{storm_label} — selected before/after pairs "
            f"(n={n_pairs}; blue = before, red = after)"
        )
    else:
        ax2.set_title(f"{storm_label} — no pairs selected")

    plt.show()


def plot_one_storm(storm_id, all_insitu_deduplicated, all_result_track_pts,
                   radius_km, delta_days,
                   all_insitu_results=None, all_paired=None):
    """Plot a single storm by ID."""
    if storm_id not in all_insitu_deduplicated:
        print(f"Storm ID '{storm_id}' not found. "
              f"Available: {sorted(all_insitu_deduplicated.keys())}")
        return
    raw_profiles = (all_insitu_results.get(storm_id)
                    if all_insitu_results is not None
                    else all_insitu_deduplicated[storm_id])
    plot_storm_insitu(
        storm_id,
        all_result_track_pts[storm_id],
        raw_profiles,
        radius_km=radius_km,
        delta_days=delta_days,
        all_paired_storm=all_paired.get(storm_id) if all_paired else None,
    )


def plot_storm_profile_pairs(storm_label, track_points, storm_profiles, max_plots=50,
                             time_linestyles=None, depth_ylim=500, all_paired=None,
                             show_mld=True):
    """
    For each paired track point, plot temperature (left), salinity (center),
    and potential density (right).
    Encoding: same color per pair across all four elements.
      - Before profile:  dashed (--), thicker
      - Before MLD:      dotted (:),  thicker  — same color as before profile
      - After profile:   solid  (-),  standard
      - After MLD:       solid  (-),  standard — same color as after profile
    Thicker 'before' lines remain visible when they coincide with 'after' values.
    """
    if time_linestyles is None:
        time_linestyles = ['-', '--', ':', '-.']
    PAIR_COLORS = plt.cm.tab10.colors
    pair_counter = 0
    for i, step_data in enumerate(storm_profiles):
        if not step_data:
            continue
        track_time = parser.parse(track_points[i]['time']).replace(tzinfo=None)
        track_lon  = track_points[i]['center'][0]
        track_lat  = track_points[i]['center'][1]
        if all_paired is not None:
            entry = all_paired.get(storm_label, {}).get(i)
            if entry is None:
                continue
            ds          = [k for k in entry if k not in (
                           'track_lon', 'track_lat', 'track_time',
                           'pair_dist_km', 'pair_dt_hours', 'pair_dt_days')][0]
            before_data = [entry[ds]['before']]
            after_data  = [entry[ds]['after']]
            pair_dist   = entry['pair_dist_km']
            pair_dt     = entry['pair_dt_hours']
        else:
            before_data, after_data = [], []
            for p in get_flat_profiles(step_data):
                try:
                    if not (p.hasvar('temperature') and p.hasvar('pressure')):
                        continue
                    d_time = pd.to_datetime(p.timestamp).to_pydatetime().replace(tzinfo=None)
                    (before_data if d_time < track_time else after_data).append(p)
                except Exception:
                    continue
            pair_dist = pair_dt = None
        if not before_data or not after_data:
            continue
        pair_counter += 1
        if pair_counter > max_plots:
            break
        color = PAIR_COLORS[pair_counter % len(PAIR_COLORS)]

        # ── Profile dates for the legend (Before / After entries) ─────────────
        def _date_str(d):
            ts = d.get('timestamp') if isinstance(d, dict) else getattr(d, 'timestamp', None)
            if ts is None:
                return ''
            try:
                return pd.to_datetime(ts).strftime('%Y-%m-%d')
            except Exception:
                return ''
        _b_date = _date_str(before_data[0])
        _a_date = _date_str(after_data[0])
        before_label = f"Before ({_b_date})" if _b_date else "Before"
        after_label  = f"After ({_a_date})"  if _a_date  else "After"

        fig, (ax_T, ax_S, ax_PD) = plt.subplots(1, 3, figsize=(14, 7), sharey=True)
        def get_vals(d, varname):
            if isinstance(d, dict):
                pres = d['pressure']
                vals = d.get(varname)
            else:
                pres = np.array(d.getvar('pressure'), dtype=float)
                vals = d.getvar(varname)
                if vals is None:
                    vals = d.getvar('ctd_' + varname)
                if vals is not None:
                    vals = np.array(vals, dtype=float)
            return pres, vals
        def plot_pair(ax, varname, xlabel):
            # before: dashed and thicker so it stays visible if it coincides with after
            for d, ls, lw, lbl in [
                (before_data[0], '--', 2.5, before_label),
                (after_data[0],  '-',  2.0, after_label),
            ]:
                pres, vals = get_vals(d, varname)
                if vals is None or np.all(np.isnan(vals)):
                    continue
                ax.plot(vals, pres, color=color, linestyle=ls,
                        linewidth=lw, label=lbl)
            ax.set_xlabel(xlabel)
            if depth_ylim is not None:
                ax.set_ylim(depth_ylim, 0)
            else:
                ax.invert_yaxis()
            handles, _ = ax.get_legend_handles_labels()
            if handles:
                ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)
        plot_pair(ax_T,  'temperature',       'Temperature (°C)')
        plot_pair(ax_S,  'salinity',          'Absolute Salinity (g/kg)')
        plot_pair(ax_PD, 'potential_density', 'Potential Density (kg/m³)')
        # ── MLD horizontal markers ────────────────────────────────────────────
        # Same pair color as the profiles.
        # Before MLD: dotted (:), thicker — matches before profile style.
        # After  MLD: solid  (-), standard — matches after profile style.
        if show_mld:
            def get_mld(d):
                if isinstance(d, dict):
                    return d.get('MLD', np.nan)
                return float(_scalar_mld(d) or np.nan)
            mld_b = get_mld(before_data[0])
            mld_a = get_mld(after_data[0])
            for ax in (ax_T, ax_S, ax_PD):
                if not np.isnan(mld_b):
                    ax.axhline(mld_b, color=color, linestyle=':',
                               linewidth=2.5, alpha=0.9,
                               label=f'MLD before ({mld_b:.0f} dbar)')
                if not np.isnan(mld_a):
                    ax.axhline(mld_a, color=color, linestyle='-',
                               linewidth=1.5, alpha=0.7,
                               label=f'MLD after ({mld_a:.0f} dbar)')
                ax.legend(fontsize=8)
        ax_T.set_ylabel('Pressure (dbar)')
        title = (f"{storm_label} | {track_time.strftime('%Y-%m-%d %H:%M')} | "
                 f"{track_lat:.2f}°N, {track_lon:.2f}°E")
        if pair_dist is not None:
            title += f" | pair dist {pair_dist} km | Δt {pair_dt:.1f} h"
        fig.suptitle(title, fontsize=10)
        plt.tight_layout()
        plt.show()

def plot_mean_across_storms(all_paired, depth_ylim=500):
    """
    Plot mean delta profiles across all storms and pairs.
    Two panels: ΔT (left) and ΔS (right).
    """
    all_delta_T, all_delta_S = [], []
    all_pressure = None

    for storm_id, step_dict in all_paired.items():
        for step_idx, entry in step_dict.items():
            ds = [k for k in entry if k not in (
                  'track_lon', 'track_lat', 'track_time',
                  'pair_dist_km', 'pair_dt_hours', 'pair_dt_days')][0]

            delta = entry[ds]['after_minus_before']

            if all_pressure is None:
                all_pressure = delta['pressure']

            all_delta_T.append(delta['temperature'])
            all_delta_S.append(delta['salinity'])

    if not all_delta_T:
        print("No pairs found in all_paired.")
        return

    mean_dT = np.nanmean(np.stack(all_delta_T, axis=0), axis=0)
    mean_dS = np.nanmean(np.stack(all_delta_S, axis=0), axis=0)

    n_pairs  = sum(len(sd) for sd in all_paired.values())
    n_storms = len(all_paired)

    fig, (ax_dT, ax_dS) = plt.subplots(1, 2, figsize=(10, 7), sharey=True,
                                        constrained_layout=True)

    ax_dT.plot(mean_dT, all_pressure, color='black', linewidth=2.5)
    ax_dT.axvline(0, color='gray', linewidth=0.8, linestyle='--')
    ax_dT.set_xlabel('ΔT (°C)')

    ax_dS.plot(mean_dS, all_pressure, color='black', linewidth=2.5)
    ax_dS.axvline(0, color='gray', linewidth=0.8, linestyle='--')
    ax_dS.set_xlabel('ΔS (g/kg)')

    for ax in (ax_dT, ax_dS):
        if depth_ylim is not None:
            ax.set_ylim(depth_ylim, 0)
        else:
            ax.invert_yaxis()
        ax.grid(True, alpha=0.3)

    ax_dT.set_ylabel('Pressure (dbar)')
    fig.suptitle(
        f"Mean ΔT and ΔS across all TCs — {n_pairs} pairs from {n_storms} storms\n"
        f"({', '.join(all_paired.keys())})",
        fontsize=11
    )
    plt.show()