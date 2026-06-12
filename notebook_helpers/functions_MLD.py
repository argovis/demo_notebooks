import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.stats import sem
from argovisHelpers import helpers as avh
from argovisHelpers import analysis as ava
import gsw
from notebook_helpers.functions import compare_profiles

def _arr(val):
    """Safely convert Profile.getvar() output to a plain float ndarray."""
    return val.filled(np.nan) if hasattr(val, 'filled') else np.asarray(val, dtype=float)


def plot_profiles_with_mld(platform_profiles, profiles_with_MLD, varname='temperature',
                           mld_key='MLD', ylim_max=400, title=None, figsize=(12, 7),
                           colors=None, ax=None, pressure_key='pressure',
                           mld_marker='s', mld_markersize=10, mld_marker_alpha=1.0,
                           helper_mld_key=None, helper_mld_marker='*'):
    """
    Plot all profiles on a single panel with optional MLD markers.

    Parameters
    ----------
    platform_profiles : dict
        The nested platform_profiles structure.
    profiles_with_MLD : list of dict
        Each dict has 'id' and optionally the key given by mld_key.
    varname : str
        Variable to plot (e.g. 'temperature', 'salinity', 'potential_density', 'N2').
    mld_key : str or None
        Key in profiles_with_MLD for the student MLD depth. If None, no student
        markers are plotted.
    ylim_max : float
        Maximum depth shown on y-axis (dbar).
    title : str or None
        Plot title. If None, auto-generated from varname.
    figsize : tuple
        Figure size (only used if ax is None).
    colors : list of str or None
        List of colors, one per profile. If None, uses a default palette.
    ax : matplotlib Axes or None
        If provided, plot on this axes. If None, create a new figure.
    pressure_key : str
        Name of the pressure variable to use (default 'pressure').
        Use 'N2_pressure' for variables on the mid-point grid.
    mld_marker : str
        Marker style for student MLD markers. Default 's' (square).
    mld_markersize : float
        Size of the student MLD marker. Default 10.
    mld_marker_alpha : float
        Transparency of the student MLD marker face (0=transparent, 1=opaque).
        Default 1.0. Reduce to ~0.5 to reveal a helper marker behind it.
    helper_mld_key : str or None
        Attribute name on the Profile object for a helper MLD depth (e.g.
        'MLD_temperature', 'MLD_potdensity'). If None, no helper markers plotted.
    helper_mld_marker : str
        Marker style for helper MLD markers. Default '*' (star).

    Returns
    -------
    all_handles : list
        Legend handles for the plotted profiles and MLD markers.
    """
    if colors is None:
        colors = ['#e41a1c', '#377eb8', '#4daf4a', '#ff7f00', '#984ea3']

    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=figsize)

    # Build profile lookup
    profile_lookup = {}
    for ds, regions_list in platform_profiles.items():
        for region_list in regions_list:
            for time_list in region_list:
                for p in time_list:
                    profile_lookup[p.rawdata['_id']] = p

    all_handles = []

    for i, entry in enumerate(profiles_with_MLD):
        profile_id = entry['id']
        profile = profile_lookup.get(profile_id)
        if not profile:
            continue

        color = colors[i % len(colors)]

        var  = _arr(profile.getvar(varname))
        pres = _arr(profile.getvar(pressure_key))

        if var.ndim == 0 or pres.ndim == 0:
            continue
        valid = ~np.isnan(var) & ~np.isnan(pres)
        if not valid.any():
            continue

        line, = ax.plot(var[valid], pres[valid], color=color, linewidth=1.5,
                        label=f"{profile_id} ({profile.latitude:.1f}°N)")
        all_handles.append(line)

        order = np.argsort(pres[valid])

        # Student MLD marker (■ square by default)
        if mld_key and entry.get(mld_key):
            mld_depth = float(entry[mld_key])
            val_at_mld = np.interp(mld_depth, pres[valid][order], var[valid][order])
            ax.plot(val_at_mld, mld_depth, mld_marker,
                    color=color, markersize=mld_markersize,
                    alpha=mld_marker_alpha,
                    markeredgewidth=1.5,
                    markeredgecolor='black', zorder=10)

        # Helper MLD marker (★ or ● depending on caller)
        if helper_mld_key:
            helper_mld = profile.getvar(helper_mld_key)
            if helper_mld is not None and helper_mld[0] is not None:
                helper_depth = float(helper_mld[0])
                val_at_helper = np.interp(helper_depth, pres[valid][order], var[valid][order])
                ax.plot(val_at_helper, helper_depth, helper_mld_marker,
                        color=color, markersize=12, markeredgewidth=1.2,
                        markeredgecolor='black', zorder=9)

    # Legend entries
    if mld_key and any(entry.get(mld_key) for entry in profiles_with_MLD):
        all_handles.append(plt.Line2D([], [], marker=mld_marker,
                                       color='gray', markersize=mld_markersize,
                                       alpha=mld_marker_alpha,
                                       markeredgecolor='black', linestyle='None',
                                       label='Your MLD estimate (■)'))

    if helper_mld_key:
        label = '★ Temperature MLD' if helper_mld_marker == '*' else '● Pot. density MLD'
        all_handles.append(plt.Line2D([], [], marker=helper_mld_marker,
                                       color='gray', markersize=12,
                                       markeredgecolor='black', linestyle='None',
                                       label=label))

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


def compute_derived_properties(p, temperature_key, salinity_key, ds):
    """Compute SA, CT, potential density, MLD (3 methods), N², and oscillation period."""
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

    # MLD from in-situ density (Δρ = 0.125 kg/m³)
    if density is not None:
        mld_dens, _ = ava.MLD_estimate(pressure, density,
                                       threshold_delta=0.125,
                                       reference_pressure=10)
        p.setvar('MLD_density', [mld_dens])
    else:
        p.setvar('MLD_density', [None])

    # MLD from potential density (Δσθ = 0.03 kg/m³)
    if pot_density is not None:
        mld_potdens, _ = ava.MLD_estimate(pressure, pot_density,
                                          threshold_delta=0.03,
                                          reference_pressure=10)
        p.setvar('MLD_potdensity', [mld_potdens])
    else:
        p.setvar('MLD_potdensity', [None])

    p.setvar('MLD', p.getvar('MLD_potdensity'))

    # N² and oscillation period
    SA_arr  = _arr(SA)
    CT_arr  = _arr(CT)
    pres_arr = _arr(pressure)
    valid = ~np.isnan(SA_arr) & ~np.isnan(CT_arr) & ~np.isnan(pres_arr)
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
    
def compute_dT_MLD_surf(profiles_with_MLD, platform_profiles):
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
    Build a styled summary table comparing student vs algorithmic MLD estimates.

    Parameters
    ----------
    platform_profiles : dict
        The nested platform_profiles structure.
    profiles_with_MLD_student : list of dict
        Each dict has 'id', 'MLD', and optionally 'dT_MLD_Surf'.

    Returns
    -------
    df_table : pd.DataFrame
        Unstyled DataFrame (sorted by latitude).
    """
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
                        mld_depth_val = float(mld_pd[0])
                        n_ml_levels   = max(1, int(np.round(mld_depth_val / 5)) + 1)
                        ml_sal_val    = round(float(np.nanmean(sal[:n_ml_levels])), 2)
                    else:
                        ml_sal_val = np.nan
                    rows.append({
                        'Profile ID':            pid,
                        'Latitude':              round(p.latitude, 3),
                        'Longitude':             round(p.longitude, 3),
                        'MLD (m)\nStudent':      entry.get('MLD'),
                        'MLD (m)\nTemperature':  round(float(mld_t[0]),  1) if mld_t  and mld_t[0]  is not None else np.nan,
                        'ΔT\n(°C)':             entry.get('dT_MLD_Surf'),
                        'MLD (m)\nPot. Density': round(float(mld_pd[0]), 1) if mld_pd and mld_pd[0] is not None else np.nan,
                        'ML Salinity\n(psu)':    ml_sal_val,
                    })
    df_table = pd.DataFrame(rows).sort_values('Latitude')
    styled = df_table.style.set_properties(**{
        'max-width': '110px',
        'white-space': 'pre-wrap',
        'text-align': 'center',
    })
    display(styled)
    print("\nCompare your visual MLD estimates (Student) with the algorithmic estimates.")
    print("Since you estimated MLD from the temperature profile, the 'Temperature' column")
    print("is your most direct comparison. Differences with the potential-density-based estimates")
    print("reveal the role of salinity in setting the mixed-layer structure.")
    print("\nML Salinity is the mean salinity over all interpolated levels from the surface")
    print("down to the potential-density mixed layer depth.")
    return df_table


def plot_mld_vs_latitude(platform_profiles, profiles_with_MLD_student):
    """
    Plot MLD vs latitude comparing student estimates with algorithmic methods.

    Parameters
    ----------
    platform_profiles : dict
        The nested platform_profiles structure.
    profiles_with_MLD_student : list of dict
        Each dict has 'id' and 'MLD'.
    """
    # Student data
    df_student = pd.DataFrame(profiles_with_MLD_student)
    df_student['MLD'] = pd.to_numeric(df_student['MLD'])
    df_student['lat'] = [
        next(r.latitude for ds, rl in platform_profiles.items()
             for r2 in rl for tl in r2 for r in tl
             if r.rawdata['_id'] == p['id'])
        for p in profiles_with_MLD_student
    ]
    df_student = df_student.sort_values('lat')

    # Helper data from platform_profiles
    helper_rows = []
    for ds, regions_list in platform_profiles.items():
        for region_list in regions_list:
            for time_list in region_list:
                for p in time_list:
                    mld_t  = p.getvar('MLD_temperature')
                    mld_pd = p.getvar('MLD_potdensity')
                    helper_rows.append({
                        'latitude':    p.latitude,
                        'temperature': float(mld_t[0])  if mld_t  and mld_t[0]  is not None else np.nan,
                        'potdensity':  float(mld_pd[0]) if mld_pd and mld_pd[0] is not None else np.nan,
                    })
    df_helper = pd.DataFrame(helper_rows).sort_values('latitude')

    plt.figure(figsize=(8, 6))
    plt.plot(df_helper['latitude'], df_helper['temperature'],
             'o--', color='lightsalmon', markersize=6, linewidth=1,
             label='Helper – Temperature (ΔT = 0.2 °C)')
    plt.plot(df_helper['latitude'], df_helper['potdensity'],
             's--', color='lightsteelblue', markersize=6, linewidth=1,
             label='Helper – Pot. Density (Δσθ = 0.03)')
    plt.plot(df_student['lat'], df_student['MLD'], 'D-', color='black',
             markerfacecolor='tomato', markersize=9, linewidth=2,
             label='Your MLD')
    plt.xlabel('Latitude')
    plt.ylabel('MLD (dbar)')
    plt.title('Mixed Layer Depth vs Latitude')
    plt.gca().invert_yaxis()
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend(fontsize=9)
    plt.show()


def inspect_float_at_depths(selected_float, depth_1, depth_2, g=9.81):
    """
    Print T and S at two target depths for a given float.

    Parameters
    ----------
    selected_float : Profile object
        Float profile with pressure, temperature, and salinity data.
    depth_1, depth_2 : float
        Target depths in dbar.
    g : float
        Gravitational acceleration (m/s²).

    Returns
    -------
    results : dict
        Dictionary with indices, pressures, T, S at the two depths.
    """
    pressure_arr = _arr(selected_float.getvar('pressure'))
    temp_arr = _arr(selected_float.getvar('temperature'))
    sal_arr = _arr(selected_float.getvar('salinity'))

    idx_1 = np.argmin(np.abs(pressure_arr - depth_1))
    idx_2 = np.argmin(np.abs(pressure_arr - depth_2))

    p_1 = pressure_arr[idx_1]
    p_2 = pressure_arr[idx_2]
    t_1 = temp_arr[idx_1]
    s_1 = sal_arr[idx_1]
    t_2 = temp_arr[idx_2]
    s_2 = sal_arr[idx_2]

    print(f"{'Depth':>10}  {'Pressure':>10}  {'T (°C)':>10}  {'S (psu)':>10}")
    print("-" * 45)
    print(f"{'Depth 1':>10}  {p_1:>10.1f}  {t_1:>10.4f}  {s_1:>10.4f}")
    print(f"{'Depth 2':>10}  {p_2:>10.1f}  {t_2:>10.4f}  {s_2:>10.4f}")

    return {
        'idx_1': idx_1, 'idx_2': idx_2,
        'p_1': p_1, 'p_2': p_2,
        't_1': t_1, 's_1': s_1,
        't_2': t_2, 's_2': s_2,
        'pressure_arr': pressure_arr,
        'temp_arr': temp_arr,
        'sal_arr': sal_arr
    }


def n2_to_period(n2):
    """Convert N² (s⁻²) to oscillation period in minutes."""
    return (2 * np.pi / np.sqrt(n2)) / 60 if n2 is not None and n2 > 0 else np.nan

