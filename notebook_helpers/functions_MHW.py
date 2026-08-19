import numpy as np
import xarray as xr
import gsw
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.patches import Rectangle
from argovisHelpers import analysis as ava

from notebook_helpers.functions_candidates_for_functions import remove_seasonal_cycle


def linear_detrend(da, time_dim='timestamp'):
    """Remove a least-squares linear trend along `time_dim`, NaN-safe, per other dimension."""
    t = ((da[time_dim] - da[time_dim][0]) / np.timedelta64(1, 'D')).astype(float)
    def _dt(y, x):
        m = np.isfinite(y)
        if m.sum() < 2:
            return y - np.nanmean(y)
        b = np.polyfit(x[m], y[m], 1)
        return y - np.polyval(b, x)
    return xr.apply_ufunc(
        _dt, da, t,
        input_core_dims=[[time_dim], [time_dim]],
        output_core_dims=[[time_dim]],
        vectorize=True,
    )


def compute_linear_trend(da, time_dim='timestamp', per='year'):
    """Least-squares linear trend of `da` along `time_dim`, NaN-safe, vectorized
    over every other dimension (e.g. longitude/latitude for a trend map) —
    returns the slope only (see `linear_detrend` for the detrended residual).

    Parameters
    ----------
    da       : xr.DataArray, must have a `time_dim` coordinate
    time_dim : str           name of the time dimension
    per      : 'year' or 'day' — units of the returned slope

    Returns
    -------
    xr.DataArray — same shape as `da` minus `time_dim`; slope in <da units> / `per`
    """
    t = ((da[time_dim] - da[time_dim][0]) / np.timedelta64(1, 'D')).astype(float)
    scale = 365.25 if per == 'year' else 1.0
    def _slope(y, x):
        m = np.isfinite(y)
        if m.sum() < 2:
            return np.nan
        b = np.polyfit(x[m], y[m], 1)
        return b[0] * scale
    return xr.apply_ufunc(
        _slope, da, t,
        input_core_dims=[[time_dim], [time_dim]],
        vectorize=True,
    )


def deseason_detrend(da, time_dim='timestamp'):
    """Deseasonalized + linearly-detrended anomaly (climatology removed, then trend removed)."""
    return linear_detrend(remove_seasonal_cycle(da, time_dim=time_dim), time_dim=time_dim)


def wind_speed_to_stress(ws, rho_air=1.225, Cd=1.3e-3):
    """Wind-stress magnitude from scalar wind speed: tau = rho_air * Cd * |U|^2  (N/m^2).
    NOTE: magnitude only — stress curl requires vector (u, v) winds."""
    return rho_air * Cd * ws**2


def compute_ohc_profile(CT, sigma0, depth_max=300, level_dim='level', cp=3991.867):
    """Upper-ocean heat content per unit area (J/m^2), integrated 0–depth_max.
    Inputs are (time × level) area-mean fields; level coordinate in dbar ≈ m."""
    CTs  = CT.sel({level_dim: slice(0, depth_max)})
    rhos = sigma0.sel({level_dim: slice(0, depth_max)})
    return cp * (rhos * CTs).integrate(level_dim)


def derive_TS_grid(T_prof, S_prof, lon0, lat0, level_dim='level'):
    """SA, CT, sigma0 on a (time × level) area-mean T/S grid via TEOS-10.

    Returns
    -------
    SA, CT, sigma0 : xr.DataArray   same shape and coords as T_prof
    """
    levels = T_prof[level_dim].values
    P2d    = np.broadcast_to(levels, T_prof.shape)
    SA     = gsw.SA_from_SP(S_prof.values, P2d, lon0, lat0)
    CT     = gsw.CT_from_t(SA, T_prof.values, P2d)
    sig0   = gsw.sigma0(SA, CT) + 1000
    wrap   = lambda a, name: xr.DataArray(a, coords=T_prof.coords, dims=T_prof.dims, name=name)
    return wrap(SA, 'absolute_salinity'), wrap(CT, 'conservative_temperature'), wrap(sig0, 'sigma0')


def mld_timeseries(sigma0, level_dim='level', time_dim='timestamp',
                   threshold_delta=0.03, reference_pressure=10):
    """MLD(t) from a (time × level) sigma0 field using ava.MLD_estimate per timestep."""
    levels = sigma0[level_dim].values
    mld = np.array([
        ava.MLD_estimate(levels, sigma0.isel({time_dim: i}).values,
                         threshold_delta=threshold_delta,
                         reference_pressure=reference_pressure)[0]
        for i in range(sigma0.sizes[time_dim])
    ], dtype=float)
    return xr.DataArray(mld, coords={time_dim: sigma0[time_dim]}, dims=[time_dim], name='MLD')


def n2_timeseries(SA, CT, lat0, level_dim='level', time_dim='timestamp'):
    """N²(time × level_mid) from (time × level) SA/CT via gsw.Nsquared per timestep."""
    levels = SA[level_dim].values
    cols   = []
    for i in range(SA.sizes[time_dim]):
        sa, ct = SA.isel({time_dim: i}).values, CT.isel({time_dim: i}).values
        n2, _  = gsw.Nsquared(sa, ct, levels, lat=lat0)
        cols.append(n2)
    p_mid_ref = gsw.Nsquared(SA.isel({time_dim: 0}).values,
                              CT.isel({time_dim: 0}).values, levels, lat=lat0)[1]
    return xr.DataArray(np.vstack(cols),
                        coords={time_dim: SA[time_dim], 'level_mid': p_mid_ref},
                        dims=[time_dim, 'level_mid'], name='N2')



def plot_sst_evolution_maps(sst_grid, map_dates, lat_band, lon_band,
                             title='SST evolution — 2013–2016 NEP MHW',
                             ncol=3, cmap='RdYlBu_r', cbar_label='SST (°C)',
                             symmetric_cbar=False):
    """Plot SST snapshots at selected dates with averaging-box overlay.

    Parameters
    ----------
    sst_grid       : xr.DataArray   (timestamp × latitude × longitude)
    map_dates      : list of str    ISO date strings, e.g. '2013-09-15'
    lat_band       : tuple          (lat_min, lat_max) of the averaging box overlay
    lon_band       : tuple          (lon_min, lon_max) of the averaging box overlay
    title          : str            figure suptitle
    ncol           : int            number of columns in the subplot grid
    cmap           : str            colormap name
    cbar_label     : str            colorbar label
    symmetric_cbar : bool           if True, set vmin/vmax symmetrically around 0
    """
    nrow = int(np.ceil(len(map_dates) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5 * ncol, 3.2 * nrow),
                              subplot_kw={'projection': ccrs.PlateCarree()}, squeeze=False)
    if symmetric_cbar:
        vmax = max(abs(float(sst_grid.min())), abs(float(sst_grid.max())))
        vmin = -vmax
    else:
        vmin, vmax = float(sst_grid.min()), float(sst_grid.max())
    for k, d in enumerate(map_dates):
        ax   = axes[k // ncol][k % ncol]
        snap = sst_grid.sel(timestamp=d, method='nearest').transpose('latitude', 'longitude')
        pc   = ax.pcolormesh(snap.longitude, snap.latitude, snap,
                             cmap=cmap, vmin=vmin, vmax=vmax,
                             transform=ccrs.PlateCarree(), shading='auto')
        ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=2)
        ax.coastlines(zorder=3)
        ax.add_patch(Rectangle((lon_band[0], lat_band[0]),
                               lon_band[1] - lon_band[0], lat_band[1] - lat_band[0],
                               fill=False, edgecolor='black', linewidth=2,
                               transform=ccrs.PlateCarree(), zorder=4))
        ax.set_title(str(snap.timestamp.values)[:10])
        gl = ax.gridlines(draw_labels=True, alpha=0.2)
        gl.top_labels = gl.right_labels = False
    for k in range(len(map_dates), nrow * ncol):
        axes[k // ncol][k % ncol].set_visible(False)
    fig.colorbar(pc, ax=axes, orientation='vertical', shrink=0.7, label=cbar_label)
    plt.suptitle(title, y=1.02, fontweight='bold')
    plt.show()


def plot_surface_timeseries(series, time_dim='timestamp', fontsize=13,
                            figsize=None, suptitle=None):
    """Two-row figure: total (top) and deseasonalised + detrended anomaly (bottom).

    Parameters
    ----------
    series   : dict   {label: xr.DataArray}  — each DataArray must have a `time_dim` coord
    time_dim : str    name of the time dimension
    fontsize : int    base font size for titles, axis labels, and tick labels
    figsize  : (width, height), optional — defaults to (5 * len(series), 8)
    suptitle : str, optional — overall figure title (e.g. a region name)
    """
    if figsize is None:
        figsize = (5 * len(series), 8)
    fig, axes = plt.subplots(2, len(series), figsize=figsize, sharex=True, squeeze=False)
    for j, (name, da) in enumerate(series.items()):
        da.plot(ax=axes[0][j], color='C0')
        axes[0][j].set_title(name, fontsize=fontsize)
        axes[0][j].set_xlabel('')
        anom = deseason_detrend(da, time_dim=time_dim)
        anom.plot(ax=axes[1][j], color='C3')
        axes[1][j].axhline(0, color='gray', lw=0.5)
        axes[1][j].set_title(f'{name} — anomaly', fontsize=fontsize)
        axes[1][j].set_xlabel('')

    # Adaptive year spacing: aim for ~7 x-axis tick labels regardless of the
    # time span, so a short window or a multi-decade one both stay legible.
    all_times = np.concatenate([np.asarray(da[time_dim].values) for da in series.values()])
    span_years = (all_times.max() - all_times.min()) / np.timedelta64(365, 'D')
    year_interval = max(1, int(np.ceil(span_years / 7)))

    for ax in axes.flat:
        ax.xaxis.set_major_locator(mdates.YearLocator(year_interval))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
        ax.xaxis.set_minor_locator(mdates.MonthLocator(bymonth=[4, 7, 10]))
        ax.tick_params(labelsize=fontsize - 1)
        ax.yaxis.label.set_size(fontsize - 1)
    axes[0][0].set_ylabel('total', fontsize=fontsize)
    axes[1][0].set_ylabel('deseason. + detrended', fontsize=fontsize)
    fig.autofmt_xdate()

    if suptitle is not None:
        plt.tight_layout(rect=[0, 0, 1, 0.93])
        fig.suptitle(suptitle, fontsize=fontsize + 2, fontweight='bold')
    else:
        plt.tight_layout()
    plt.show()


def _hovmoller_panel(ax, da, time_dim='timestamp', level_name='level', cmap='RdYlBu_r',
                     label='', overlay_mld=None, ylim=(300, 0)):
    """Single Hovmöller contourf panel with optional MLD time-series overlay."""
    pc = ax.contourf(da[time_dim], da[level_name], da.transpose(level_name, time_dim),
                     levels=21, cmap=cmap, extend='both')
    if overlay_mld is not None:
        ax.plot(overlay_mld[time_dim], overlay_mld, color='black', lw=1.8, label='MLD')
        ax.legend(loc='lower left', fontsize=8)
    ax.set_ylim(ylim[0], ylim[1])
    ax.set_ylabel('Depth (dbar)')
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax.xaxis.set_minor_locator(mdates.MonthLocator(bymonth=[4, 7, 10]))
    plt.colorbar(pc, ax=ax, label=label)


def plot_depth_time_hovmoller(fields, mld_ts, time_dim='timestamp', ylim=(300, 0)):
    """Plot depth–time Hovmöller panels (total + deseason./detrended anomaly) for multiple fields.

    Parameters
    ----------
    fields  : list of (label, DataArray, level_dim_name, cmap)
              MLD is overlaid only on panels whose level_dim_name == 'level'.
    mld_ts  : xr.DataArray   MLD time series to overlay on level-dimension panels
    ylim    : tuple           (depth_max, depth_min) — y-axis is plotted inverted
    """
    fig, axes = plt.subplots(len(fields), 2, figsize=(15, 4 * len(fields)))
    for i, (name, da, lvl, cmap) in enumerate(fields):
        mld_for = mld_ts if lvl == 'level' else None
        _hovmoller_panel(axes[i][0], da, time_dim=time_dim, level_name=lvl,
                         cmap=cmap, label=name, overlay_mld=mld_for, ylim=ylim)
        axes[i][0].set_title(f'{name} — total')
        anom = deseason_detrend(da)
        _hovmoller_panel(axes[i][1], anom, time_dim=time_dim, level_name=lvl,
                         cmap='RdBu_r', label=f'{name} anomaly',
                         overlay_mld=mld_for, ylim=ylim)
        axes[i][1].set_title(f'{name} — deseasonalized + detrended')
    fig.autofmt_xdate()
    plt.tight_layout()
    plt.show()