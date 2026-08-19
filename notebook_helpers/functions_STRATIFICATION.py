"""
functions_STRATIFICATION.py
────────────────────────────
All helper functions for the Ocean Stratification & MLD notebook (ATOC 1060).

Student-facing wrappers (call these in the notebook):
    query_floats            — fetch profiles from Argovis
    compute_derived         — compute potential density + MLD estimates
    plot_float_map          — map of float locations
    plot_profiles           — clean 3-panel T / S / PD plot (no MLD markers)
    plot_tspd_with_mld      — 3-panel plot with student MLD markers
    plot_mld_vs_latitude    — MLD vs latitude (student estimates only)
    print_mld_table         — summary table of student estimates

Low-level functions (used internally or in the advanced KK notebook):
    plot_profiles_with_mld  — single-panel profile plot with optional MLD markers
    compute_derived_properties — full derived-quantity computation for one profile
    compute_dT_MLD_surf     — ΔT between surface and student MLD
    build_mld_summary_table — styled table comparing student vs algorithmic MLD
    inspect_float_at_depths — print T and S at two target depths
    n2_to_period            — convert N² to oscillation period (minutes)
"""

# ── Imports ───────────────────────────────────────────────────────────────────

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import gsw
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from argovisHelpers import helpers as avh
from argovisHelpers import analysis as ava
import warnings
warnings.filterwarnings('ignore')


# ─────────────────────────────────────────────────────────────────────────────
# Internal utilities
# ─────────────────────────────────────────────────────────────────────────────

def _arr(val):
    """Safely convert Profile.getvar() output to a plain float ndarray."""
    return val.filled(np.nan) if hasattr(val, 'filled') else np.asarray(val, dtype=float)


def _lat_label(lat):
    """Format latitude as e.g. '45.2°N' or '60.1°S'."""
    return f"{abs(lat):.1f}°{'N' if lat >= 0 else 'S'}"


def _get_profiles(platform_profiles):
    """Return a flat list of all profiles from the nested platform_profiles dict."""
    profiles = []
    for ds, regions_list in platform_profiles.items():
        for region_list in regions_list:
            for time_list in region_list:
                profiles.extend(time_list)
    return profiles


# ─────────────────────────────────────────────────────────────────────────────
# 1. Data retrieval
# ─────────────────────────────────────────────────────────────────────────────
    
def query_floats(float_ids, apikey='guest',
                 apiroot='https://argovis-api.colorado.edu/'):
    """
    Fetch Argo profiles from Argovis and return them as a platform_profiles dict.
    Parameters
    ----------
    float_ids : list of str
        Profile IDs in 'PLATFORM_CYCLE' format, e.g. ['7902209_042', ...].
    apikey    : str   (default 'guest')
    apiroot   : str
    Returns
    -------
    platform_profiles : nested dict
        Compatible with all functions in this file.
    """
    from notebook_helpers.functions import traverse_query, traverse_interpolate
    INTERP       = list(range(0, 2001, 5))
    target       = [{'id': fid} for fid in float_ids]
    time_periods = [{'startDate': '2020-01-01T00:00:00Z',
                     'endDate':   '2026-12-31T00:00:00Z'}]
    collections  = {'argo': 'temperature,1,salinity,1'}  # ',1' = Argo QC flag 1 (good data only)
    raw = traverse_query(target, time_periods, collections,
                         avh.queryProfile, apikey=apikey, apiroot=apiroot)
    platform_profiles = traverse_interpolate(raw, INTERP)
    loaded = _get_profiles(platform_profiles)
    print(f'Fetched {len(loaded)}/{len(float_ids)} profiles:')
    for p in loaded:
        print(f'  ✓  {p.rawdata["_id"]:22s}  '
              f'lat = {p.latitude:+6.2f}°   lon = {p.longitude:+7.2f}°')
    return platform_profiles


# ─────────────────────────────────────────────────────────────────────────────
# 2. Derived quantities
# ─────────────────────────────────────────────────────────────────────────────

def compute_derived(platform_profiles):
    """
    Compute potential density, in-situ density, algorithmic MLD (3 methods),
    N², and oscillation period for every profile in platform_profiles.

    Calls compute_derived_properties() on each profile internally.
    """
    temperature_key = {'argo': 'temperature'}
    salinity_key    = {'argo': 'salinity'}
    for ds, regions_list in platform_profiles.items():
        for region_list in regions_list:
            for time_list in region_list:
                for p in time_list:
                    compute_derived_properties(
                        p, temperature_key, salinity_key, ds)
    print('✓ Derived properties computed for all profiles.')


def compute_derived_properties(p, temperature_key, salinity_key, ds):
    """
    Compute and store on profile p:
      absolute_salinity, conservative_temperature, potential_density,
      density, MLD_temperature, MLD_density, MLD_potdensity, MLD,
      N², N2_pressure, oscillation_period.

    Parameters
    ----------
    p               : Profile object
    temperature_key : dict  e.g. {'argo': 'temperature'}
    salinity_key    : dict  e.g. {'argo': 'salinity'}
    ds              : str   dataset key, e.g. 'argo'
    """
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

    # MLD from temperature (ΔT = 0.2 °C from 10 dbar reference)
    if temperature is not None:
        mld_temp, _ = ava.MLD_estimate(pressure, temperature,
                                       threshold_delta=-0.2,
                                       reference_pressure=10)
        p.setvar('MLD_temperature', [mld_temp])
    else:
        p.setvar('MLD_temperature', [None])

    # MLD from in-situ density (Δρ = 0.125 kg m⁻³)
    if density is not None:
        mld_dens, _ = ava.MLD_estimate(pressure, density,
                                       threshold_delta=0.125,
                                       reference_pressure=10)
        p.setvar('MLD_density', [mld_dens])
    else:
        p.setvar('MLD_density', [None])

    # MLD from potential density (Δσθ = 0.03 kg m⁻³)
    if pot_density is not None:
        mld_potdens, _ = ava.MLD_estimate(pressure, pot_density,
                                          threshold_delta=0.03,
                                          reference_pressure=10)
        p.setvar('MLD_potdensity', [mld_potdens])
    else:
        p.setvar('MLD_potdensity', [None])

    p.setvar('MLD', p.getvar('MLD_potdensity'))

    # N² and oscillation period
    SA_arr   = _arr(SA)
    CT_arr   = _arr(CT)
    pres_arr = _arr(pressure)
    valid    = ~np.isnan(SA_arr) & ~np.isnan(CT_arr) & ~np.isnan(pres_arr)
    SA_v, CT_v, pres_v = SA_arr[valid], CT_arr[valid], pres_arr[valid]
    order = np.argsort(pres_v)
    N2, p_mid = gsw.Nsquared(SA_v[order], CT_v[order], pres_v[order],
                              lat=p.latitude)
    good = ~np.isnan(N2)
    p.setvar('N2', N2[good])
    p.setvar('N2_pressure', p_mid[good])
    with np.errstate(divide='ignore', invalid='ignore'):
        period = np.where(N2[good] > 0,
                          (2 * np.pi / np.sqrt(N2[good])) / 60, np.nan)
    p.setvar('oscillation_period', period)


# ─────────────────────────────────────────────────────────────────────────────
# 3. Map
# ─────────────────────────────────────────────────────────────────────────────

def plot_float_map(platform_profiles, colors, label=''):
    """
    Plot float locations on a global Robinson-projection map.
    Dots are color-coded to match the profile plots.
    """
    profiles = _get_profiles(platform_profiles)
    fig, ax  = plt.subplots(
        figsize=(13, 5),
        subplot_kw={'projection': ccrs.Robinson(central_longitude=180)}
    )
    ax.stock_img()
    ax.add_feature(cfeature.COASTLINE, linewidth=0.6)
    ax.add_feature(cfeature.BORDERS,   linewidth=0.3, linestyle=':')
    gl = ax.gridlines(draw_labels=True, linewidth=0.4, color='gray',
                      alpha=0.5, x_inline=False, y_inline=False)
    gl.top_labels   = False
    gl.right_labels = False

    handles = []
    for p, c in zip(profiles, colors):
        ax.scatter(p.longitude, p.latitude, color=c, s=120,
                   zorder=5, transform=ccrs.PlateCarree())
        ax.annotate(_lat_label(p.latitude),
                    xy=(p.longitude, p.latitude), xytext=(6, 4),
                    textcoords='offset points', fontsize=8, color=c,
                    transform=ccrs.PlateCarree())
        handles.append(mpatches.Patch(
            color=c,
            label=f"{p.rawdata['_id']}  ({_lat_label(p.latitude)})"))

    ax.legend(handles=handles, loc='lower left', fontsize=8)
    ax.set_title(f'Float Locations — {label}' if label else 'Float Locations',
                 fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.show()


# ─────────────────────────────────────────────────────────────────────────────
# 4. Single-panel profile plot (low-level, used by the wrappers below)
# ─────────────────────────────────────────────────────────────────────────────

def plot_profiles_with_mld(platform_profiles, profiles_with_MLD,
                            varname='temperature', mld_key='MLD',
                            ylim_max=400, title=None, figsize=(12, 7),
                            colors=None, ax=None, pressure_key='pressure',
                            mld_marker='s', mld_markersize=10,
                            mld_marker_alpha=1.0,
                            helper_mld_key=None, helper_mld_marker='*'):
    """
    Plot all profiles on a single axes panel with optional MLD markers.

    Parameters
    ----------
    platform_profiles  : nested dict  — output of query_floats()
    profiles_with_MLD  : list of dict — each dict has 'id' and optionally
                         the key given by mld_key
    varname            : str  — variable to plot
                         ('temperature', 'salinity', 'potential_density', 'N2', …)
    mld_key            : str or None  — dict key for student MLD depth;
                         None = no student markers
    ylim_max           : float — max depth on y-axis (dbar)
    title              : str or None — panel title (auto-generated if None)
    figsize            : tuple — used only when ax is None
    colors             : list of str or None
    ax                 : matplotlib Axes or None
    pressure_key       : str — 'pressure' or 'N2_pressure' for mid-point grid
    mld_marker         : str — marker style for student MLD (default 's')
    mld_markersize     : float
    mld_marker_alpha   : float — 0 = transparent, 1 = opaque
    helper_mld_key     : str or None — Profile attribute for algorithmic MLD
    helper_mld_marker  : str — marker style for algorithmic MLD (default '*')

    Returns
    -------
    all_handles : list  — legend handles
    """
    if colors is None:
        colors = ['#e41a1c', '#377eb8', '#4daf4a', '#ff7f00', '#984ea3']

    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=figsize)

    # Build profile lookup by ID
    profile_lookup = {}
    for ds, regions_list in platform_profiles.items():
        for region_list in regions_list:
            for time_list in region_list:
                for p in time_list:
                    profile_lookup[p.rawdata['_id']] = p

    all_handles = []

    for i, entry in enumerate(profiles_with_MLD):
        profile_id = entry['id']
        profile    = profile_lookup.get(profile_id)
        if not profile:
            continue

        color = colors[i % len(colors)]
        var   = _arr(profile.getvar(varname))
        pres  = _arr(profile.getvar(pressure_key))

        if var.ndim == 0 or pres.ndim == 0:
            continue
        valid = ~np.isnan(var) & ~np.isnan(pres)
        if not valid.any():
            continue

        line, = ax.plot(var[valid], pres[valid], color=color, linewidth=1.5,
                        label=f"{profile_id} ({_lat_label(profile.latitude)})")
        all_handles.append(line)

        order = np.argsort(pres[valid])

        # Student MLD marker
        if mld_key and entry.get(mld_key):
            mld_depth  = float(entry[mld_key])
            val_at_mld = np.interp(mld_depth, pres[valid][order], var[valid][order])
            ax.plot(val_at_mld, mld_depth, mld_marker,
                    color=color, markersize=mld_markersize,
                    alpha=mld_marker_alpha,
                    markeredgewidth=1.5, markeredgecolor='black', zorder=10)

        # Algorithmic MLD marker
        if helper_mld_key:
            helper_mld = profile.getvar(helper_mld_key)
            if helper_mld is not None and helper_mld[0] is not None:
                helper_depth   = float(helper_mld[0])
                val_at_helper  = np.interp(helper_depth,
                                           pres[valid][order], var[valid][order])
                ax.plot(val_at_helper, helper_depth, helper_mld_marker,
                        color=color, markersize=12, markeredgewidth=1.2,
                        markeredgecolor='black', zorder=9)

    # Legend entries for markers
    if mld_key and any(entry.get(mld_key) for entry in profiles_with_MLD):
        all_handles.append(plt.Line2D([], [], marker=mld_marker,
                                       color='gray', markersize=mld_markersize,
                                       alpha=mld_marker_alpha,
                                       markeredgecolor='black', linestyle='None',
                                       label='Your MLD estimate'))
    if helper_mld_key:
        lbl = '★ Temperature MLD' if helper_mld_marker == '*' else '● Pot. density MLD'
        all_handles.append(plt.Line2D([], [], marker=helper_mld_marker,
                                       color='gray', markersize=12,
                                       markeredgecolor='black', linestyle='None',
                                       label=lbl))

    if title is None:
        title = varname.replace('_', ' ').title()
    ax.set_xlabel(title, fontsize=12)
    ax.set_ylabel('Pressure (dbar)', fontsize=12)
    ax.set_title(title, fontsize=13)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(ylim_max, 0)

    if standalone:
        ax.legend(handles=all_handles, fontsize=8, loc='lower left')
        plt.tight_layout()
        plt.show()

    return all_handles


# ─────────────────────────────────────────────────────────────────────────────
# 5. Three-panel profile plots
# ─────────────────────────────────────────────────────────────────────────────

def plot_profiles(platform_profiles, colors, max_depth=700, label=''):
    """
    Three-panel plot: Temperature | Salinity | Potential Density.
    Each line = one float / latitude.  No MLD markers.
    Use this to visually identify the mixed layer before entering estimates.
    """
    empty = [{'id': p.rawdata['_id']} for p in _get_profiles(platform_profiles)]

    fig, (ax_t, ax_s, ax_pd) = plt.subplots(1, 3, figsize=(14, 10.5), sharey=True)

    handles = plot_profiles_with_mld(
        platform_profiles, empty, varname='temperature', mld_key=None,
        title='Temperature (°C)', colors=colors,
        ax=ax_t, ylim_max=max_depth)
    plot_profiles_with_mld(
        platform_profiles, empty, varname='salinity', mld_key=None,
        title='Salinity (psu)', colors=colors,
        ax=ax_s, ylim_max=max_depth)
    plot_profiles_with_mld(
        platform_profiles, empty, varname='potential_density', mld_key=None,
        title='Potential Density (kg m⁻³)', colors=colors,
        ax=ax_pd, ylim_max=max_depth)

    ax_t.legend(handles=handles, fontsize=9, loc='lower left')
    ax_s.set_ylabel('')
    ax_pd.set_ylabel('')

    suptitle = 'Temperature, Salinity, and Potential Density'
    if label:
        suptitle += f' — {label}'
    suptitle += '\n(each line = one float / latitude; pressure ≈ depth in meters)'
    fig.suptitle(suptitle, fontsize=12, fontweight='bold')
    plt.tight_layout()
    plt.show()


def plot_tspd_with_mld(platform_profiles, profiles_with_MLD, colors,
                        max_depth=700, label=''):
    """
    Three-panel T / S / PD plot with student MLD markers:
      ■  square on Temperature panel   = your MLD_T  estimate
      ■  square on Potential Density panel = your MLD_PD estimate

    profiles_with_MLD must be a list of dicts with keys:
        'id', 'MLD_T', 'MLD_PD'
    """
    fig, (ax_t, ax_s, ax_pd) = plt.subplots(1, 3, figsize=(14, 10.5), sharey=True)

    handles = plot_profiles_with_mld(
        platform_profiles, profiles_with_MLD,
        varname='temperature', mld_key='MLD_T',
        title='Temperature (°C)', colors=colors,
        ax=ax_t, ylim_max=max_depth)
    plot_profiles_with_mld(
        platform_profiles, profiles_with_MLD,
        varname='salinity', mld_key=None,
        title='Salinity (psu)', colors=colors,
        ax=ax_s, ylim_max=max_depth)
    plot_profiles_with_mld(
        platform_profiles, profiles_with_MLD,
        varname='potential_density', mld_key='MLD_PD',
        title='Potential Density (kg m⁻³)', colors=colors,
        ax=ax_pd, ylim_max=max_depth)

    ax_t.legend(handles=handles, fontsize=9, loc='lower left')
    ax_s.set_ylabel('')
    ax_pd.set_ylabel('')

    suptitle = 'T, S, PD with Your MLD Estimates'
    if label:
        suptitle += f' — {label}'
    fig.suptitle(suptitle, fontsize=12, fontweight='bold')
    plt.tight_layout()
    plt.show()


# ─────────────────────────────────────────────────────────────────────────────
# 6. MLD vs. latitude
# ─────────────────────────────────────────────────────────────────────────────

def plot_mld_vs_latitude(platform_profiles, profiles_with_MLD, colors, label=''):
    """
    Scatter plot of Mixed Layer Depth vs. latitude — student estimates only.

    profiles_with_MLD : list of dicts with keys 'id', 'MLD_T', 'MLD_PD'.
    Squares (■) = MLD_T.   Circles (●) = MLD_PD.
    """
    lat_map = {p.rawdata['_id']: p.latitude
               for p in _get_profiles(platform_profiles)}

    fig, ax = plt.subplots(figsize=(8, 5))
    for entry, c in zip(profiles_with_MLD, colors):
        lat = lat_map.get(entry['id'])
        if lat is None:
            continue
        if entry.get('MLD_T') is not None:
            ax.scatter(lat, float(entry['MLD_T']),
                       color=c, s=130, marker='s', zorder=5)
        if entry.get('MLD_PD') is not None:
            ax.scatter(lat, float(entry['MLD_PD']),
                       color=c, s=130, marker='o', zorder=5)

    ax.legend(handles=[
        Line2D([0],[0], marker='s', color='gray', markerfacecolor='gray',
               markersize=9, linestyle='None', label='MLD$_T$ — your estimate'),
        Line2D([0],[0], marker='o', color='gray', markerfacecolor='gray',
               markersize=9, linestyle='None', label='MLD$_{PD}$ — your estimate'),
    ], fontsize=9, loc='upper right')

    ax.invert_yaxis()
    ax.set_xlabel('Latitude  (°N — negative = South)', fontsize=12)
    ax.set_ylabel('Mixed Layer Depth (dbar)', fontsize=12)
    ax.axvline(0, color='k', linewidth=0.8, linestyle='--', alpha=0.4)
    ax.grid(True, linestyle=':', alpha=0.5)
    ax.set_title(f'Mixed Layer Depth vs. Latitude — {label}' if label else
                 'Mixed Layer Depth vs. Latitude',
                 fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.show()


# ─────────────────────────────────────────────────────────────────────────────
# 7. Summary table
# ─────────────────────────────────────────────────────────────────────────────

def print_mld_table(platform_profiles, profiles_with_MLD, label=''):
    """
    Return a DataFrame summarising student MLD_T and MLD_PD estimates.

    profiles_with_MLD : list of dicts with keys 'id', 'MLD_T', 'MLD_PD'.
    """
    lat_map = {p.rawdata['_id']: p.latitude
               for p in _get_profiles(platform_profiles)}
    rows = []
    for entry in profiles_with_MLD:
        lat = lat_map.get(entry['id'], float('nan'))
        rows.append({
            'Float ID'           : entry['id'],
            'Latitude'           : _lat_label(lat),
            'Your MLD_T (dbar)'  : entry.get('MLD_T'),
            'Your MLD_PD (dbar)' : entry.get('MLD_PD'),
        })
    if label:
        print(f'\n── MLD Summary: {label} ──')
    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────────────────
# 8. Advanced / KK-notebook functions (kept for compatibility)
# ─────────────────────────────────────────────────────────────────────────────

def compute_dT_MLD_surf(profiles_with_MLD, platform_profiles):
    """
    Compute ΔT between the surface and the student-estimated MLD for each profile.

    profiles_with_MLD : list of dicts with keys 'id' and 'MLD' (single depth).

    Returns a list of dicts with 'id', 'MLD', 'dT_MLD_Surf'.
    """
    ds_lookup_id  = {}
    ds_lookup_lat = {}
    for ds, regions_list in platform_profiles.items():
        for region_list in regions_list:
            for time_list in region_list:
                for p in time_list:
                    ds_lookup_id[p.rawdata['_id']] = p
                    ds_lookup_lat[round(p.latitude, 3)] = p

    results = []
    for profile in profiles_with_MLD:
        if isinstance(profile, dict):
            pid       = profile['id']
            mld_depth = float(profile['MLD'])
            p_obj     = ds_lookup_id.get(pid)
        else:
            lat, mld_depth = profile
            p_obj = ds_lookup_lat.get(round(lat, 3))
            pid   = p_obj.rawdata['_id'] if p_obj else None

        dT = None
        if p_obj is not None:
            temps     = _arr(p_obj.getvar('temperature'))
            pressures = _arr(p_obj.getvar('pressure'))
            valid     = ~np.isnan(temps) & ~np.isnan(pressures)
            vt, vp    = temps[valid], pressures[valid]
            if len(vp) > 1:
                idx    = np.argsort(vp)
                t_surf = vt[idx][0]
                t_mld  = np.interp(mld_depth, vp[idx], vt[idx])
                dT     = round(float(t_mld - t_surf), 3)

        results.append({'id': pid, 'MLD': mld_depth, 'dT_MLD_Surf': dT})
    return results


def build_mld_summary_table(platform_profiles, profiles_with_MLD_student):
    """
    Build and display a styled summary table comparing student vs algorithmic MLD.

    profiles_with_MLD_student : list of dicts with 'id', 'MLD', 'dT_MLD_Surf'.

    Returns the unstyled DataFrame sorted by latitude.
    """
    from IPython.display import display as ipy_display

    ds_lookup = {p['id']: p for p in profiles_with_MLD_student}
    rows = []
    for ds, regions_list in platform_profiles.items():
        for region_list in regions_list:
            for time_list in region_list:
                for p in time_list:
                    pid   = p.rawdata['_id']
                    entry = ds_lookup.get(pid, {})
                    mld_t  = p.getvar('MLD_temperature')
                    mld_pd = p.getvar('MLD_potdensity')
                    sal    = p.getvar('salinity')
                    if sal is not None and mld_pd is not None and mld_pd[0] is not None:
                        n_ml      = max(1, int(np.round(float(mld_pd[0]) / 5)) + 1)
                        ml_sal    = round(float(np.nanmean(sal[:n_ml])), 2)
                    else:
                        ml_sal = np.nan
                    rows.append({
                        'Profile ID'          : pid,
                        'Latitude'            : round(p.latitude, 3),
                        'Longitude'           : round(p.longitude, 3),
                        'MLD (m)\nStudent'    : entry.get('MLD'),
                        'MLD (m)\nTemperature': round(float(mld_t[0]),  1)
                                                if mld_t  and mld_t[0]  is not None else np.nan,
                        'ΔT\n(°C)'            : entry.get('dT_MLD_Surf'),
                        'MLD (m)\nPot. Density': round(float(mld_pd[0]), 1)
                                                 if mld_pd and mld_pd[0] is not None else np.nan,
                        'ML Salinity\n(psu)'  : ml_sal,
                    })
    df = pd.DataFrame(rows).sort_values('Latitude')
    styled = df.style.set_properties(**{
        'max-width'   : '110px',
        'white-space' : 'pre-wrap',
        'text-align'  : 'center',
    })
    ipy_display(styled)
    print("\nCompare your visual MLD estimates (Student) with the algorithmic estimates.")
    print("Since you estimated MLD from the temperature profile, the 'Temperature' column")
    print("is your most direct comparison. Differences with the potential-density-based")
    print("estimates reveal the role of salinity in setting the mixed-layer structure.")
    print("\nML Salinity is the mean salinity from the surface down to MLD_potdensity.")
    return df


def inspect_float_at_depths(selected_float, depth_1, depth_2, g=9.81):
    """
    Print temperature and salinity at two target depths for a given profile.

    Parameters
    ----------
    selected_float : Profile object
    depth_1, depth_2 : float  — target depths in dbar
    g : float  — gravitational acceleration (m s⁻²), default 9.81

    Returns
    -------
    dict with keys idx_1, idx_2, p_1, p_2, t_1, s_1, t_2, s_2,
                   pressure_arr, temp_arr, sal_arr
    """
    pressure_arr = _arr(selected_float.getvar('pressure'))
    temp_arr     = _arr(selected_float.getvar('temperature'))
    sal_arr      = _arr(selected_float.getvar('salinity'))

    idx_1 = np.argmin(np.abs(pressure_arr - depth_1))
    idx_2 = np.argmin(np.abs(pressure_arr - depth_2))

    print(f"{'Depth':>10}  {'Pressure':>10}  {'T (°C)':>10}  {'S (psu)':>10}")
    print("-" * 45)
    print(f"{'Depth 1':>10}  {pressure_arr[idx_1]:>10.1f}  "
          f"{temp_arr[idx_1]:>10.4f}  {sal_arr[idx_1]:>10.4f}")
    print(f"{'Depth 2':>10}  {pressure_arr[idx_2]:>10.1f}  "
          f"{temp_arr[idx_2]:>10.4f}  {sal_arr[idx_2]:>10.4f}")

    return dict(idx_1=idx_1, idx_2=idx_2,
                p_1=pressure_arr[idx_1], p_2=pressure_arr[idx_2],
                t_1=temp_arr[idx_1],     s_1=sal_arr[idx_1],
                t_2=temp_arr[idx_2],     s_2=sal_arr[idx_2],
                pressure_arr=pressure_arr,
                temp_arr=temp_arr,
                sal_arr=sal_arr)


def n2_to_period(n2):
    """Convert N² (s⁻²) to oscillation period in minutes."""
    return (2 * np.pi / np.sqrt(n2)) / 60 if n2 is not None and n2 > 0 else np.nan