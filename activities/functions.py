import warnings
from datetime import datetime, timezone

import cartopy.crs as ccrs
import gsw
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy
import numpy as np

from argovisHelpers import helpers as avh


def traverse_query(space, time, collections, queryhelper, apikey, apiroot,
                   qsp_base={}, kwargs_base={}):
    """
    Iterate over space/time/collection combinations and collect API results.

    Parameters
    ----------
    space : list of dict
        Spatial query parameter dicts.
    time : list of dict
        Temporal query parameter dicts.
    collections : dict
        Maps route -> data value.
    queryhelper : callable
        Function that executes a single API call.
    apikey : str
        Argovis API key.
    apiroot : str
        Argovis API root URL.
    qsp_base : dict, optional
        Base query string parameters merged into every request.
    kwargs_base : dict, optional
        Base keyword arguments merged into every queryhelper call.
    """
    results = {}
    for route, data in collections.items():
        results[route] = []
        for r_i, r in enumerate(space):
            results[route].append([])
            for t_i, t in enumerate(time):
                if "/meta" in route:
                    qsp = qsp_base | r | t
                else:
                    qsp = qsp_base | r | t | {"data": data}
                kwargs = {'route': route, 'options': qsp, 'apikey': apikey, 'apiroot': apiroot} | kwargs_base
                results[route][r_i].append(queryhelper(**kwargs))
    return results


def plot_maps(lats, lons, title="Map", margin=10, fig_settings=None,
              labels=None, colors=None, markers=None, markersize=5,
              edgecolors='black', linewidths=0.4):
    """
    Plot scatter locations on a map.

    Encoding: COLOR -> dataset, MARKER -> time period.

    Parameters
    ----------
    lats, lons : list or list-of-lists
        Coordinates to plot. Pass a list-of-lists for multiple groups.
    markers : list of str, optional
        Marker cycle for groups. Defaults to ['o', '^', 'p', 'D', 'v', '*'].
    colors : list, optional
        Color cycle for groups. Defaults to tab10 colormap.
    """
    if fig_settings is None:
        fig_settings = {}
    if markers is None:
        markers = ['o', '^', 'p', 'D', 'v', '*']
    if colors is None:
        colors = list(plt.cm.tab10.colors)
    if len(lats) == 0:
        return
    ax = fig_settings.get('ax')
    show_plot = fig_settings.get('show_plot', ax is None)
    if ax is None:
        figsize = fig_settings.get('figsize', (10, 8))
        fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree()})
        ax.stock_img()
        ax.coastlines(resolution='50m', color='black', linewidth=0.5)
    lats_list = lats if isinstance(lats[0], (list, np.ndarray)) else [lats]
    lons_list = lons if isinstance(lons[0], (list, np.ndarray)) else [lons]
    all_lons_flat, all_lats_flat = [], []
    for i, (lat_group, lon_group) in enumerate(zip(lats_list, lons_list)):
        label  = labels[i] if labels and i < len(labels) else None
        color  = colors[i % len(colors)]
        marker = markers[i % len(markers)]
        ax.scatter(lon_group, lat_group, color=color, marker=marker, s=markersize**2.5,
                   edgecolors=edgecolors, linewidths=linewidths,
                   transform=ccrs.PlateCarree(), zorder=10, label=label)
        all_lons_flat.extend(lon_group)
        all_lats_flat.extend(lat_group)
    min_lon = np.nanmin(all_lons_flat) - margin
    max_lon = np.nanmax(all_lons_flat) + margin
    min_lat = max(np.nanmin(all_lats_flat) - margin, -90)
    max_lat = min(np.nanmax(all_lats_flat) + margin,  90)
    ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = True
    if labels:
        ax.legend(loc='lower left', fontsize=7, framealpha=0.7)
    ax.set_title(title)
    if show_plot:
        plt.show()
    return ax

def average_profiles_dict(profiles_sel, varname_lookup=None, use_nanmean=False):
    averaged = {}
    for ds, time_list in profiles_sel.items():
        if not time_list:
            averaged[ds] = []
            continue
        var_names = time_list[0].variable_names()
        pressure  = numpy.asarray(time_list[0].getvar('pressure'))
        n_levels  = len(pressure)
        means = {}
        for var in var_names:
            arrays = []
            for p in time_list:
                try:
                    raw_val = p.getvar(var)
                    if isinstance(raw_val, numpy.ma.MaskedArray):
                        val = raw_val.filled(numpy.nan).astype(float)
                    else:
                        val = numpy.asarray(raw_val, dtype=float)
                    if val.shape == (n_levels,):
                        arrays.append(val)
                except Exception:
                    continue
            if arrays:
                stacked = numpy.stack(arrays)
                if use_nanmean:
                    means[var] = numpy.nanmean(stacked, axis=0)
                else:
                    # Only include levels where all profiles have valid values
                    mask = numpy.any(numpy.isnan(stacked), axis=0)
                    result = numpy.mean(stacked, axis=0)
                    result[mask] = numpy.nan
                    means[var] = result
        keys   = list(means.keys())
        values = list(means.values())
        timestamps = [p.rawdata.get('timestamp') for p in time_list if p.rawdata.get('timestamp')]
        mid_ts = timestamps[len(timestamps) // 2] if timestamps else None
        if mid_ts is not None and not isinstance(mid_ts, str):
            mid_ts = mid_ts.strftime('%Y-%m-%dT%H:%M:%SZ')
        raw = {
            'data_info': [keys, [''] * len(keys), [[] for _ in keys]],
            'data'     : values,
            'geolocation': {'coordinates': [None, None]},
            'timestamp': mid_ts,
            'id'       : f'{ds}_mean',
        }
        averaged[ds] = [avh.Profile(raw)]
    return averaged

def compare_profiles(ts_profiles_sel, varname, dataset_colors, time_linestyles,
                     depth_ylim=None, varname_lookup=None, dataset=None,
                     region_label=None, datestring=None, mean_style=False,
                     ax=None, markersize=None, markevery=None, **plot_kwargs):
    """
    Plot vertical profiles for one or more datasets.
 
    Averaging is NOT done inside this function. To plot mean profiles, pre-average
    using average_profiles_dict() and pass the result here with mean_style=True.
 
    Parameters
    ----------
    ts_profiles_sel : dict
        Dict keyed by dataset name, each value a list of Profile objects.
    varname : str
        Variable name to plot on the x-axis.
    dataset_colors : dict
        Maps dataset name -> color string.
    time_linestyles : list of str
        Linestyle cycle for datasets.
    depth_ylim : int, float, or None
        Maximum pressure (dbar) on the y-axis. None plots full depth.
    varname_lookup : dict, optional
        Maps dataset name -> actual variable name in that dataset.
    dataset : str, optional
        Single dataset key to plot. If None, all keys are used.
    region_label : str, optional
        Region label appended to the plot title.
    datestring : str, optional
        Date/period string for the title. Auto-derived from timestamps if None.
    mean_style : bool
        If True, use mean-profile aesthetics (bold line + markers).
        If False, use individual-profile aesthetics (thin, transparent lines).
    ax : matplotlib Axes, optional
        Existing axes to draw into. If None, a new figure is created.
    markersize : int, float, or None
        Marker size for mean-profile markers. Default None (markers off).
        Only takes effect when markevery is also set.
    markevery : int or None
        Plot a marker every N levels. Default None (markers off).
        Set to 1 to show a marker at every level, or e.g. 5 for every 5th level.
        When set, markersize defaults to 3 if not specified.
    **plot_kwargs
        Override any of: color, linestyle, marker, linewidth, alpha,
        markeredgecolor, markeredgewidth.
    """
    datasets = [dataset] if dataset is not None else list(ts_profiles_sel.keys())
    created_fig = ax is None
    if created_fig:
        fig, ax = plt.subplots(figsize=(6, 8))
    if datestring is None:
        all_timestamps = [
            p.rawdata['timestamp'] for ds in datasets
            for p in ts_profiles_sel[ds] if p.rawdata.get('timestamp')
        ]
        datestring = f"{sorted(all_timestamps)[0][:4]}-{sorted(all_timestamps)[-1][:4]}" if all_timestamps else "unknown"
    legend_handles, any_data = [], False
    for i, ds in enumerate(datasets):
        v = varname_lookup.get(ds, varname) if varname_lookup else varname
        time_list = ts_profiles_sel[ds]
        if not time_list or not time_list[0].hasvar(v):
            continue
        ds_color  = plot_kwargs.get('color',           dataset_colors.get(ds, 'black'))
        ls        = plot_kwargs.get('linestyle',       time_linestyles[i % len(time_linestyles)])
        marker    = plot_kwargs.get('marker',          'o')
        linewidth = plot_kwargs.get('linewidth',       0.5)
        alpha     = plot_kwargs.get('alpha',           0.3)
        mec       = plot_kwargs.get('markeredgecolor', None)
        mew       = plot_kwargs.get('markeredgewidth', None)
        label     = f"{ds} - {datestring}" + (f" | {region_label}" if region_label is not None else "")
        if mean_style:
            p     = time_list[0]
            _var  = numpy.asarray(p.getvar(v),         dtype=float)
            _pres = numpy.asarray(p.getvar('pressure'), dtype=float)
            ax.plot(_var, _pres,
                    color=ds_color, linestyle=ls,
                    marker=marker if markevery is not None else None,
                    linewidth=1.5,
                    markersize=markersize if markersize is not None else 3,
                    markevery=markevery,
                    markeredgecolor=mec,
                    markeredgewidth=mew)
            legend_handles.append(plt.Line2D([], [], color=ds_color, linestyle=ls,
                                             marker=marker if markevery is not None else None,
                                             markersize=markersize if markersize is not None else 3,
                                             markeredgecolor=mec,
                                             markeredgewidth=mew,
                                             label=label))
        else:
            for p in time_list:
                v_row = p.getvar(v)
                z_row = p.getvar('pressure')
                if v_row is not None and len(v_row) == len(z_row) > 0:
                    ax.plot(v_row, z_row, color=ds_color, linestyle=ls,
                            marker=marker if markevery is not None else None,
                            markevery=markevery,
                            markersize=markersize if markersize is not None else 3,
                            markeredgecolor=mec,
                            markeredgewidth=mew,
                            linewidth=linewidth, alpha=alpha)
            legend_handles.append(plt.Line2D([], [], color=ds_color, linestyle=ls,
                                             marker=marker if markevery is not None else None,
                                             markersize=markersize if markersize is not None else 3,
                                             markeredgecolor=mec,
                                             markeredgewidth=mew,
                                             linewidth=1.5, label=label))
        any_data = True
    if not any_data:
        if created_fig:
            plt.close(fig)
        return []
    if depth_ylim is not None:
        ax.set_ylim(depth_ylim, 0)
    else:
        ymin, ymax = ax.get_ylim()
        if ymin < ymax:
            ax.set_ylim(ymax, ymin)
    ax.set_title(f"{'Mean' if mean_style else 'Profiles'} | {varname} | {datestring}"
                 + (f" | {region_label}" if region_label is not None else ""))
    ax.set_xlabel(varname)
    ax.set_ylabel("Pressure (dbar)")
    ax.grid(True, linestyle=":", linewidth=0.5)
    if created_fig:
        ax.legend(handles=legend_handles, loc='best', fontsize='small')
        plt.show()
    return legend_handles

def plot_profiles_TSdiagram(valuesT, valuesS, var_nameT, var_nameS,
                             values_to_color=None, title="", fig_settings=None):
    """T-S diagram, optionally colored by a third variable (depth, time, etc.)."""
    if fig_settings is None:
        fig_settings = {}
    ax        = fig_settings.get('ax')
    show_plot = False
    if ax is None:
        figsize   = fig_settings.get('figsize', (8, 8))
        fig, ax   = plt.subplots(figsize=figsize)
        show_plot = True

    color      = fig_settings.get('color',     'black')
    marker     = fig_settings.get('marker',    '.')
    markersize = fig_settings.get('markersize', 10)
    alpha      = fig_settings.get('alpha',      0.3)

    if values_to_color is not None:
        try:
            cbar_clim = fig_settings.get('cbar_clim',
                (np.concatenate(values_to_color).min(), np.concatenate(values_to_color).max()))
        except:
            cbar_clim = fig_settings.get('cbar_clim',
                (min(values_to_color), max(values_to_color)))
        if isinstance(cbar_clim[0], datetime):
            cbar_clim = (mdates.date2num(cbar_clim[0]), mdates.date2num(cbar_clim[1]))
        cbar_label = fig_settings.get('cbar_label', None)
        cbar_cmap  = fig_settings.get('cbar_cmap',  'viridis')

    for i in range(len(valuesT)):
        v_rowT = valuesT[i]
        v_rowS = valuesS[i]
        if values_to_color:
            z_row = values_to_color[i]
            if isinstance(z_row, np.ndarray):
                if len(v_rowT) != len(z_row) and len(z_row) != 1: continue
                if len(v_rowS) != len(z_row) and len(z_row) != 1: continue
            else:
                z_row = ([mdates.date2num(z_row)] if isinstance(z_row, datetime)
                         else [z_row]) * len(v_rowT)
        if v_rowT is None or len(v_rowT) == 0: continue
        if v_rowS is None or len(v_rowS) == 0: continue
        if len(v_rowT) != len(v_rowS): continue

        if not values_to_color:
            ax.plot(v_rowS, v_rowT, color=color, marker=marker,
                    markersize=markersize, linestyle='none', alpha=alpha)
        else:
            sc = ax.scatter(v_rowS, v_rowT, s=markersize, c=z_row,
                            alpha=alpha, cmap=cbar_cmap)
            sc.set_clim(cbar_clim[0], cbar_clim[1])

    if values_to_color:
        cbar = fig.colorbar(sc, ax=ax)
        if isinstance(values_to_color[i], datetime):
            cbar.ax.yaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
        if cbar_label is not None:
            cbar.set_label(cbar_label)

    ax.set_xlabel(var_nameS)
    ax.set_ylabel(var_nameT)
    ax.set_title(title)
    ax.grid(True, linestyle='--', alpha=0.5)
    if show_plot:
        plt.show()



##### OPTIONAL
def to_utc_naive(t):
    if t is None:
        return np.datetime64("NaT")
    if t.tzinfo is not None:
        return t.astimezone(timezone.utc).replace(tzinfo=None)
    return t  # already naive

def profiles_to_xarray(profiles, var_rename=None, dtype=np.float32):
    import numpy as np
    import xarray as xr

    if len(profiles) == 0:
        raise ValueError("Empty profile list.")

    pressure = np.asarray(profiles[0].getvar('pressure'))
    n_profiles = len(profiles)
    n_levels = len(pressure)

    varnames = profiles[0].variable_names()
    if var_rename is None:
        var_rename = {}

    data_vars = {}

    for var in varnames:
        try:
            test = profiles[0].getvar(var)
        except Exception:
            continue

        if test is None:
            continue

        test = np.asarray(test)

        out_name = var_rename.get(var, var.lower())

        # --- profile variable
        if test.shape == (n_levels,):
            arr = np.full((n_profiles, n_levels), np.nan, dtype=dtype)

            for i, p in enumerate(profiles):
                try:
                    arr[i, :] = p.getvar(var)
                except Exception:
                    continue

            data_vars[out_name] = (("index", "levels"), arr)

        # --- scalar variable
        elif test.size == 1:
            arr = np.full(n_profiles, np.nan, dtype=dtype)

            for i, p in enumerate(profiles):
                try:
                    val = p.getvar(var)
                    arr[i] = np.asarray(val).item()
                except Exception:
                    continue

            data_vars[out_name] = (("index",), arr)

        # --- ignore weird shapes (optional: extend later)
        else:
            continue

    lon = np.array([p.longitude for p in profiles])
    lat = np.array([p.latitude for p in profiles])
    #time = np.array([p.timestamp for p in profiles], dtype="datetime64[ns]")
    time = np.array(
            [to_utc_naive(p.timestamp) for p in profiles],
            dtype="datetime64[s]"
    
                    )
    ds = xr.Dataset(
        data_vars=data_vars,
        coords={
            "index": np.arange(n_profiles),
            "levels": pressure,
            "longitude": ("index", lon),
            "latitude": ("index", lat),
            "time": ("index", time),
        }
    )

    ds["levels"].attrs["units"] = "dbar"

    return ds