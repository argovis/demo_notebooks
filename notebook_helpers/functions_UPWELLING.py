import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.stats import sem
from scipy.interpolate import RegularGridInterpolator
from argovisHelpers import helpers as avh
from argovisHelpers import analysis as ava


def _arr(val):
    """Safely convert Profile.getvar() output to a plain float ndarray."""
    return val.filled(np.nan) if hasattr(val, 'filled') else np.asarray(val, dtype=float)


def _smooth_grid(lon_arr, depth_arr, data_2d, n_lon=200):
    """Interpolate a (lon, depth) matrix onto a finer longitude grid."""
    valid = ~np.all(np.isnan(data_2d), axis=1)
    if valid.sum() < 2:
        return lon_arr, data_2d
    lon_v, data_v = lon_arr[valid], data_2d[valid, :]
    lon_fine = np.linspace(lon_v.min(), lon_v.max(), n_lon)
    interp = RegularGridInterpolator(
        (lon_v, depth_arr), data_v,
        method='linear', bounds_error=False, fill_value=np.nan
    )
    lg, dg = np.meshgrid(lon_fine, depth_arr, indexing='ij')
    return lon_fine, interp((lg, dg))


def plot_section(ax, combined, depth_mask, island_lon,
                 inset_lon_min, inset_lon_max, inset_depth_min, inset_depth_max,
                 label='', levels=np.arange(8, 30, 2), iso_levels=[15, 25],
                 xlim=None, ylim=(500, 0),
                 interpolated_data=None, depth_arr=None, lon_fine=None,
                 domain_cap=None, smooth=False):
    """
    Panel (a): mean temperature cross-section.

    Three modes (checked in order):
      1. interpolated_data/depth_arr/lon_fine passed explicitly  (legacy Fig 3)
      2. smooth=True  — auto-interpolate from the xarray combined grid
      3. Neither       — use xarray's built-in plotting  (Fig 4 default)
    domain_cap: if set, drop longitude bins above this value before plotting.
    """
    if interpolated_data is not None:
        # Legacy explicit path (kept for backwards compat)
        cf = ax.contourf(lon_fine, depth_arr, interpolated_data.T,
                         levels=levels, cmap='RdYlBu_r', extend='both')
        ax.contour(lon_fine, depth_arr, interpolated_data.T,
                   levels=levels, colors='black', linewidths=0.75, alpha=0.6)
        cs = ax.contour(lon_fine, depth_arr, interpolated_data.T,
                        levels=iso_levels, colors='white', linewidths=2.5)
        plt.colorbar(cf, ax=ax, label='Temperature (°C)')
    elif smooth:
        # Auto-smooth path
        if depth_mask is None:
            temp_da = combined.temperature
        elif isinstance(depth_mask, np.ndarray) and depth_mask.dtype == bool:
            temp_da = combined.temperature.isel(level=depth_mask)
        else:
            temp_da = combined.temperature.isel(level=depth_mask)

        lon_arr = temp_da.longitude.values
        d_arr   = temp_da.level.values
        t_mat   = temp_da.values  # (lon, depth)

        if domain_cap is not None:
            keep = lon_arr <= domain_cap
            lon_arr, t_mat = lon_arr[keep], t_mat[keep, :]

        lf, t_fine = _smooth_grid(lon_arr, d_arr, t_mat)
        cf = ax.contourf(lf, d_arr, t_fine.T,
                         levels=levels, cmap='RdYlBu_r', extend='both')
        ax.contour(lf, d_arr, t_fine.T,
                   levels=levels, colors='black', linewidths=0.75, alpha=0.6)
        cs = ax.contour(lf, d_arr, t_fine.T,
                        levels=iso_levels, colors='white', linewidths=2.5)
        plt.colorbar(cf, ax=ax, label='Temperature (°C)')
        if xlim is None:
            valid = ~np.all(np.isnan(t_mat), axis=1)
            xlim = (lon_arr[valid].min(), lon_arr[valid].max())
    else:
        temp_plot = combined.temperature.isel(level=depth_mask)
        temp_plot.plot.contourf(x='longitude', y='level', levels=levels,
                                cmap='RdYlBu_r', yincrease=False, ax=ax, add_colorbar=True)
        temp_plot.plot.contour(x='longitude', y='level', levels=levels,
                               colors='black', linewidths=0.75, alpha=0.6,
                               yincrease=False, ax=ax)
        cs = temp_plot.plot.contour(x='longitude', y='level', levels=iso_levels,
                                    colors='white', linewidths=2.5, yincrease=False, ax=ax)

    plt.clabel(cs, inline=True, fontsize=11, fmt='%1.0f°C')
    ax.axvline(x=island_lon, color='black', linewidth=2, zorder=5)
    ax.add_patch(Rectangle(
        xy=(inset_lon_min, inset_depth_min),
        width=inset_lon_max - inset_lon_min,
        height=inset_depth_max - inset_depth_min,
        linewidth=1.5, edgecolor='white', facecolor='none', linestyle='--', zorder=5
    ))
    ax.text(inset_lon_min + 0.15, inset_depth_max - 25, 'Inset',
            color='white', fontsize=10, zorder=6)
    if xlim:
        ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_ylabel('Depth (m)')
    ax.set_xlabel('Longitude (°E)')
    ax.set_title(f'{label} mean temperature')


def plot_anomaly_zoom(ax, combined, depth_mask,
                      inset_lon_min, inset_lon_max, inset_depth_min, inset_depth_max,
                      island_lon=None, label='', iso_levels=[15, 25],
                      anom_levels=np.arange(-0.8, 0.9, 0.1),
                      interpolated_anom=None, interpolated_temp=None,
                      depth_arr=None, lon_fine=None,
                      domain_cap=None, smooth=False):
    """
    Panel (b): zonal anomaly, zoomed on the inset region.

    Three modes (same logic as plot_section):
      1. interpolated_anom/interpolated_temp passed explicitly  (legacy)
      2. smooth=True  — auto-interpolate
      3. Neither       — xarray plotting (Fig 4 default)
    """
    if interpolated_anom is not None:
        cf = ax.contourf(lon_fine, depth_arr, interpolated_anom.T,
                         levels=anom_levels, cmap='RdBu_r', extend='both')
        ax.contour(lon_fine, depth_arr, interpolated_anom.T,
                   levels=anom_levels, colors='black', linewidths=0.75, alpha=0.6)
        cs = ax.contour(lon_fine, depth_arr, interpolated_temp.T,
                        levels=iso_levels, colors='white', linewidths=2.5)
        plt.clabel(cs, inline=True, fontsize=11, fmt='%1.0f°C')
        plt.colorbar(cf, ax=ax, label='Temperature anomaly (°C)')
    elif smooth:
        if depth_mask is None:
            temp_da = combined.temperature
        elif isinstance(depth_mask, np.ndarray) and depth_mask.dtype == bool:
            temp_da = combined.temperature.isel(level=depth_mask)
        else:
            temp_da = combined.temperature.isel(level=depth_mask)

        lon_arr = temp_da.longitude.values
        d_arr   = temp_da.level.values
        t_mat   = temp_da.values

        if domain_cap is not None:
            keep = lon_arr <= domain_cap
            lon_arr, t_mat = lon_arr[keep], t_mat[keep, :]

        valid   = ~np.all(np.isnan(t_mat), axis=1)
        zmean   = np.nanmean(t_mat[valid, :], axis=0)
        a_mat   = t_mat - zmean

        # Clip to inset region
        inset = (lon_arr >= inset_lon_min) & (lon_arr <= inset_lon_max) & valid
        lon_in, a_in, t_in = lon_arr[inset], a_mat[inset, :], t_mat[inset, :]

        lf, a_fine = _smooth_grid(lon_in, d_arr, a_in)
        _,  t_fine = _smooth_grid(lon_in, d_arr, t_in)

        cf = ax.contourf(lf, d_arr, a_fine.T,
                         levels=anom_levels, cmap='RdBu_r', extend='both')
        ax.contour(lf, d_arr, a_fine.T,
                   levels=anom_levels, colors='black', linewidths=0.75, alpha=0.6)
        cs = ax.contour(lf, d_arr, t_fine.T,
                        levels=iso_levels, colors='white', linewidths=2.5)
        plt.clabel(cs, inline=True, fontsize=11, fmt='%1.0f°C')
        plt.colorbar(cf, ax=ax, label='Temperature anomaly (°C)')
    else:
        temp_plot = combined.temperature.isel(level=depth_mask)
        anomaly   = temp_plot - temp_plot.mean(dim='longitude')
        zoom_anom = anomaly.sel(longitude=slice(inset_lon_min, inset_lon_max))
        zoom_temp = temp_plot.sel(longitude=slice(inset_lon_min, inset_lon_max))
        zoom_anom.plot.contourf(x='longitude', y='level', levels=anom_levels,
                                cmap='RdBu_r', yincrease=False, ax=ax,
                                extend='both', add_colorbar=True)
        zoom_anom.plot.contour(x='longitude', y='level', levels=anom_levels,
                               colors='black', linewidths=0.75, alpha=0.6,
                               yincrease=False, ax=ax)
        cs = zoom_temp.plot.contour(x='longitude', y='level', levels=iso_levels,
                                    colors='white', linewidths=2.5, yincrease=False, ax=ax)
        plt.clabel(cs, inline=True, fontsize=11, fmt='%1.0f°C')

    if island_lon:
        ax.axvline(x=island_lon, color='black', linewidth=2, zorder=5)
    ax.set_ylim(inset_depth_max, inset_depth_min)
    ax.set_ylabel('Depth (m)')
    ax.set_xlabel('Longitude (°E)')
    ax.set_title(f'{label} anomaly (zoom)')


def plot_profiles_comparison(ax, profiles, depth_mask, depth_plot,
                             upwelling_lon, west_of, label='',
                             xlim=(9, 30), ylim=(400, 0),
                             west_data=None, east_data=None):
    """
    Panel (c): mean profiles west of vs within upwelling, with 95% CI.
    Pass west_data/east_data as pre-computed arrays,
    or leave as None to compute from profiles.
    """
    if west_data is not None and east_data is not None:
        west_arr, east_arr = west_data, east_data
    else:
        uw_lon_min, uw_lon_max = upwelling_lon
        west_temps, east_temps = [], []
        for p in profiles:
            lon = p.rawdata['geolocation']['coordinates'][0] % 360
            lat = p.rawdata['geolocation']['coordinates'][1]
            if lon is None or lat is None:
                continue
            try:
                t = _arr(p.getvar('temperature'))[depth_mask]
            except Exception:
                continue
            if np.all(np.isnan(t)):
                continue
            if lon < west_of:
                west_temps.append(t)
            elif uw_lon_min <= lon <= uw_lon_max:
                east_temps.append(t)
        if not west_temps:
            raise ValueError(f"No profiles found west of {west_of}°E.")
        if not east_temps:
            raise ValueError(f"No profiles found in upwelling zone {upwelling_lon}°E.")
        west_arr = np.vstack(west_temps)
        east_arr = np.vstack(east_temps)

    for arr, color, lbl in [
        (west_arr, 'red',  'West of upwelling'),
        (east_arr, 'blue', 'Within upwelling'),
    ]:
        mean = np.nanmean(arr, axis=0)
        se   = 2 * sem(arr, axis=0, nan_policy='omit')
        ax.fill_betweenx(depth_plot, mean - se, mean + se,
                         color=color, alpha=0.2, zorder=3)
        ax.plot(mean, depth_plot, color=color, lw=2.5, label=lbl, zorder=4)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel('Temp (°C)')
    ax.set_ylabel('Depth (m)')
    ax.set_title(f'{label} profiles')
    ax.legend(loc='lower left')
    ax.grid(linestyle='--', alpha=0.5)