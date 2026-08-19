import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature


def plot_gcos_ohc_layers(ds, layers, labels=None, colors=None,
                         title='GCOS Earth Heat Inventory — global ocean heat content by layer',
                         ylabel='Ocean heat content anomaly (ZJ)'):
    """Plot one or more depth-layer OHC time series from the GCOS Earth Heat Inventory,
    with shaded uncertainty bands where a matching '<var>_uncertainty' variable exists.

    Parameters
    ----------
    ds     : xr.Dataset   opened GCOS_EHI_*.nc file (time-only dimension)
    layers : list of str  data variable names, e.g. 'ocean_heat_content_0-300m'
    labels : list of str, optional  legend labels (defaults to `layers`)
    colors : list, optional  color cycle (defaults to tab10)
    """
    if labels is None:
        labels = layers
    if colors is None:
        colors = list(plt.cm.tab10.colors)

    fig, ax = plt.subplots(figsize=(10, 6))
    for i, (var, label) in enumerate(zip(layers, labels)):
        da = ds[var]
        color = colors[i % len(colors)]
        ax.plot(da.time, da, color=color, linewidth=2, label=label)
        unc_name = f'{var}_uncertainty'
        if unc_name in ds:
            unc = ds[unc_name]
            ax.fill_between(da.time, da - unc, da + unc, color=color, alpha=0.15)
    ax.axhline(0, color='gray', linewidth=0.5)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(loc='upper left', fontsize=9)
    ax.grid(True, linestyle=':', linewidth=0.5)
    plt.tight_layout()
    plt.show()


def plot_regional_vs_global_ohc(regional_series, global_da, global_uncertainty=None,
                                region_colors=None, global_color='black',
                                global_label='Global (GCOS EHI)',
                                title='Regional (Argovis, on the fly) vs. global (GCOS EHI) OHC',
                                regional_ylabel='Regional OHC anomaly (J/m²)',
                                global_ylabel='Global OHC anomaly (ZJ)'):
    """Twin-axis comparison of regional on-the-fly OHC (Argovis, J/m²) against a
    single global OHC layer (GCOS Earth Heat Inventory, ZJ). Units and spatial
    scope differ by design — this compares the SHAPE/TIMING of the trend across
    regional and global scales, not the magnitude.

    Parameters
    ----------
    regional_series : dict  {region_label: xr.DataArray(time)}  e.g. from compute_ohc_profile
    global_da : xr.DataArray(time)  a GCOS layer, e.g. ds['ocean_heat_content_0-300m']
    global_uncertainty : xr.DataArray(time), optional  matching '<var>_uncertainty'
    """
    if region_colors is None:
        region_colors = list(plt.cm.tab10.colors)

    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax2 = ax1.twinx()

    for i, (label, da) in enumerate(regional_series.items()):
        time_dim = da.dims[0]
        ax1.plot(da[time_dim], da, color=region_colors[i % len(region_colors)],
                linewidth=1.5, label=label)
    ax1.set_ylabel(regional_ylabel)

    ax2.plot(global_da.time, global_da, color=global_color, linewidth=2.5,
             linestyle='--', label=global_label)
    if global_uncertainty is not None:
        ax2.fill_between(global_da.time, global_da - global_uncertainty,
                         global_da + global_uncertainty, color=global_color, alpha=0.12)
    ax2.set_ylabel(global_ylabel)

    # Zoom to the regional data's time span — otherwise a decades-long global
    # record (e.g. 1960-2020) squeezes a short regional window into a sliver.
    all_times = np.concatenate([np.asarray(da[da.dims[0]].values) for da in regional_series.values()])
    ax1.set_xlim(all_times.min(), all_times.max())

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left', fontsize=9)

    ax1.set_title(title)
    ax1.grid(True, linestyle=':', linewidth=0.5)
    plt.tight_layout()
    plt.show()


def plot_weekly_vs_annual_sst(weekly_series, annual_series,
                              title='Weekly SST vs. its annual mean'):
    """Overlay raw weekly SST (thin, light) with its annual mean (bold), one
    panel per region, to show how much clearer a long-term trend becomes once
    weekly noise is averaged away.

    Parameters
    ----------
    weekly_series, annual_series : dict  {region_label: xr.DataArray(time)}
        Must share the same keys. `annual_series` is typically built with
        `weekly_series[name].resample(timestamp='1YS').mean()`.
    """
    names = list(weekly_series.keys())
    fig, axes = plt.subplots(1, len(names), figsize=(6 * len(names), 5), squeeze=False)
    for j, name in enumerate(names):
        ax = axes[0][j]
        w, a = weekly_series[name], annual_series[name]
        w_dim, a_dim = w.dims[0], a.dims[0]
        ax.plot(w[w_dim], w, color='lightsteelblue', linewidth=0.8, alpha=0.9, label='Weekly (raw)')
        ax.plot(a[a_dim], a, color='C3', linewidth=2.5, marker='o', label='Annual mean')
        ax.set_title(name, fontsize=10)
        ax.set_ylabel('SST (°C)')
        ax.legend(fontsize=8, loc='best')
        ax.grid(True, linestyle=':', linewidth=0.5)
    plt.tight_layout(rect=[0, 0, 1, 0.92])
    fig.suptitle(title, fontweight='bold')
    plt.show()


def plot_ohc_trend_map(trend, title='OHC trend', cbar_label='OHC trend (J/m²/yr)',
                       cmap='RdBu_r'):
    """Single-panel map of a linear-trend field (e.g. from compute_linear_trend),
    on a diverging colormap centered at zero so warming/cooling read symmetrically.

    Parameters
    ----------
    trend : xr.DataArray (longitude x latitude)
    """
    vmax = float(np.nanmax(np.abs(trend.values)))
    fig, ax = plt.subplots(figsize=(9, 7), subplot_kw={'projection': ccrs.PlateCarree()})
    snap = trend.transpose('latitude', 'longitude')
    pc = ax.pcolormesh(snap.longitude, snap.latitude, snap,
                       cmap=cmap, vmin=-vmax, vmax=vmax,
                       transform=ccrs.PlateCarree(), shading='auto')
    ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=2)
    ax.coastlines(zorder=3)
    gl = ax.gridlines(draw_labels=True, alpha=0.3)
    gl.top_labels = gl.right_labels = False
    ax.set_title(title)
    plt.tight_layout()

    # Match the colorbar's height to the map's actual rendered extent. Cartopy
    # shrinks the map to its data aspect ratio *after* layout, but
    # fig.colorbar(..., shrink=...) sizes off the axes' pre-shrink box, so a
    # `shrink` value always ends up taller than the map itself. Forcing a draw
    # and reading the axes' real position side-steps that.
    fig.canvas.draw()
    pos = ax.get_position()
    cax = fig.add_axes([pos.x1 + 0.02, pos.y0, 0.02, pos.height])
    fig.colorbar(pc, cax=cax, label=cbar_label)
    plt.show()
