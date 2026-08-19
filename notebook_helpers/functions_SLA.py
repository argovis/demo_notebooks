import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

def tidy_grid(da):
    """Sort a queried gridded field (xr.DataArray on a regular lat/lon/time/level
    grid, e.g. as returned by avh.queryGrid) by every spatial/temporal dim it
    has, so downstream .sel(slice), .interp, and contourf all see
    monotonic-increasing coords.

    If the box crosses the antimeridian (e.g. a Pacific-centred query like
    120E-80W), the API returns longitude split into a negative branch
    (-180..-80) and a positive branch (120..180). Sorting those as-is puts
    the true centre of the box (the dateline) at the two edges of the plot
    and the unqueried gap in the middle. Detect that case — shifting the
    negative branch by +360 shrinks the total span — and shift before
    sorting so the coordinate is contiguous.
    """
    if 'longitude' in da.dims:
        lon = da.longitude.values
        raw_span = lon.max() - lon.min()
        shifted = np.where(lon < 0, lon + 360, lon)
        if shifted.max() - shifted.min() < raw_span:
            da = da.assign_coords(longitude=('longitude', shifted))
    dims = [d for d in ('timestamp', 'level', 'latitude', 'longitude') if d in da.dims]
    return da.sortby(dims)


def find_isotherm_depth(temp, iso=20.0, level_dim='level'):
    """Level (e.g. pressure in dbar, as for RG09 — not necessarily metres of depth;
    check the input's level_units) of the `iso`-degree isotherm, by linear
    interpolation between the two native levels that bracket the crossing —
    no fine pre-interpolation. Assumes temperature decreasing with depth and
    `level` ascending."""
    def _col(t, z):
        for k in range(len(t) - 1):
            if (t[k] - iso) >= 0 and (t[k + 1] - iso) < 0:
                f = (t[k] - iso) / (t[k] - t[k + 1])
                return z[k] + f * (z[k + 1] - z[k])
        return np.nan
    return xr.apply_ufunc(
        _col, temp, temp[level_dim],
        input_core_dims=[[level_dim], [level_dim]],
        vectorize=True,
    )
    
def plot_hovmoller(wind_anom, sst_anom_hov, sla_anom, d20_anom,
                   enso_index=None, figsize=(8, 10)):
    """
    Plot a 4-panel Hovmöller diagram: Wind | SST | SLA | D20.
    A narrow ENSO-state stripe is drawn on the left edge of each panel.
    """
    def _add_cb(mappable, ax, label):
        cb = plt.colorbar(mappable, ax=ax, orientation='horizontal', pad=0.08, label=label)
        plt.setp(cb.ax.get_xticklabels(), rotation=45, ha='right')
        return cb

    lon_min = sla_anom.longitude.min().item()
    lon_max = sla_anom.longitude.max().item()
    t_min   = sla_anom.timestamp.min()
    t_max   = sla_anom.timestamp.max()

    fig, (ax_wind, ax_sst, ax_sla, ax_d20) = plt.subplots(
        1, 4, figsize=figsize,
        gridspec_kw={'width_ratios': [2, 2, 2, 2]},
        sharey=True
    )

    # ── Panel 1: Wind anomaly ─────────────────────────────────────────────────
    c_wind = ax_wind.contourf(wind_anom.longitude, wind_anom.timestamp, wind_anom,
                               levels=np.linspace(-3, 3, 21),
                               cmap='RdBu_r', extend='both')
    _add_cb(c_wind, ax_wind, 'Zonal Wind Anomaly (m/s)')
    ax_wind.set_title('Wind anom', fontsize=12)
    ax_wind.set_xlabel('Longitude')
    ax_wind.set_ylabel('Time')

    # ── Panel 2: SST anomaly ──────────────────────────────────────────────────
    if sst_anom_hov is not None:
        c_sst = ax_sst.contourf(sst_anom_hov.longitude, sst_anom_hov.timestamp, sst_anom_hov,
                                 levels=np.linspace(-2.5, 2.5, 21),
                                 cmap='RdBu_r', extend='both')
        ax_sst.contour(sst_anom_hov.longitude, sst_anom_hov.timestamp, sst_anom_hov,
                       levels=[-0.5, 0.5], colors=['blue', 'red'],
                       linewidths=0.8, linestyles='--')
        _add_cb(c_sst, ax_sst, 'SST Anomaly (°C)')
    else:
        ax_sst.text(0.5, 0.5, 'SST not\navailable', ha='center', va='center',
                    transform=ax_sst.transAxes, fontsize=10, color='gray')
    ax_sst.set_title('SST anom', fontsize=12)
    ax_sst.set_xlabel('Longitude')

    # ── Panel 3: SLA anomaly ──────────────────────────────────────────────────
    c_sla = ax_sla.contourf(sla_anom.longitude, sla_anom.timestamp, sla_anom,
                             levels=np.linspace(-0.2, 0.2, 21),
                             cmap='RdBu_r', extend='both')
    _add_cb(c_sla, ax_sla, 'SLA (m)')
    ax_sla.set_title('SLA anom', fontsize=12)
    ax_sla.set_xlabel('Longitude')

    # ── Panel 4: D20 anomaly ──────────────────────────────────────────────────
    c_d20 = ax_d20.contourf(d20_anom.longitude, d20_anom.timestamp, d20_anom,
                             levels=np.linspace(-50, 50, 21),
                             cmap='RdBu_r', extend='both')
    levels_lines = np.arange(-60, 61, 20)
    cs_d20 = ax_d20.contour(d20_anom.longitude, d20_anom.timestamp, d20_anom,
                             levels=levels_lines, colors='black', linewidths=0.6, alpha=0.6)
    ax_d20.contour(d20_anom.longitude, d20_anom.timestamp, d20_anom,
                   levels=[0], colors='black', linewidths=1.8)
    ax_d20.clabel(cs_d20, inline=True, fontsize=8, fmt='%1.0f')
    _add_cb(c_d20, ax_d20, '$D_{20}$ Anomaly (dbar)')
    ax_d20.set_title('$D_{20}$ anom', fontsize=12)
    ax_d20.set_xlabel('Longitude')

    # ── Shared time axis (earliest at top) ────────────────────────────────────
    ax_wind.invert_yaxis()
    ax_wind.set_ylim(t_max, t_min)

    plt.suptitle('Equatorial Pacific Hovmöller Diagrams', fontsize=14, y=1.01)
    plt.tight_layout()

    # ── ENSO stripe on the left edge of each panel (MUST be last) ─────────────
    if enso_index is not None:
        vals = enso_index.values
        rgb  = np.array([
            [1.0, 0.0, 0.0] if v >= 0.5 else
            [0.0, 0.0, 1.0] if v <= -0.5 else
            [0.85, 0.85, 0.85]
            for v in vals
        ]).reshape(-1, 1, 3)

        for ax in [ax_wind, ax_sst, ax_sla, ax_d20]:
            pos = ax.get_position()                     # read AFTER all colorbars + tight_layout
            ax_stripe = fig.add_axes([pos.x0 - 0.01, pos.y0, 0.003, pos.height])
            ax_stripe.axis('off')
            ax_stripe.imshow(rgb, aspect='auto', origin='upper', extent=[0, 1, 0, 1])

        ax_wind.tick_params(axis='y', pad=5)

    plt.show()

def plot_temperature_sections(temp_eq, dates_to_compare,
                              depth_max=200, figsize=(12, 10)):
    """
    Plot longitude–depth temperature sections + 20°C isotherm for a set of
    snapshot dates, stacked vertically with a shared colorbar.

    Parameters
    ----------
    temp_eq          : xr.DataArray  – temperature averaged over the equatorial
                                       band (time × level × longitude)
    dates_to_compare : dict          – {label: date_str} e.g.
                                       {'Before': '2013-01-15', ...}
    depth_max        : int           – maximum level (dbar, RG09's native pressure
                                       coordinate — not metres of depth) to show
                                       (default 200)
    figsize          : tuple         – figure size (width, height)
    """

    fig, axes = plt.subplots(
        nrows=len(dates_to_compare), ncols=1,
        figsize=figsize, sharex=True, sharey=True
    )

    # Make iterable even for a single snapshot
    if len(dates_to_compare) == 1:
        axes = [axes]

    for ax, (label, date_str) in zip(axes, dates_to_compare.items()):

        section      = temp_eq.sel(timestamp=date_str, method='nearest')
        section_plot = section.transpose('level', 'longitude')

        # Temperature field
        cf = ax.contourf(section_plot.longitude, section_plot.level, section_plot,
                         levels=np.arange(8, 32, 1),
                         cmap='RdYlBu_r', extend='both')

        # 20°C isotherm
        cs = ax.contour(section_plot.longitude, section_plot.level, section_plot,
                        levels=[20], colors='black', linewidths=3)
        ax.clabel(cs, fmt='%1.0f°C', inline=True, fontsize=10)

        ax.set_title(f'Equatorial Pacific Section: {label}', fontsize=14)
        ax.set_ylabel('Pressure (dbar)')
        ax.invert_yaxis()
        ax.set_ylim(depth_max, 0)

    axes[-1].set_xlabel('Longitude (Degrees East)')

    # Shared colorbar
    cbar = fig.colorbar(cf, ax=axes,
                        orientation='vertical', pad=0.02, shrink=0.8)
    cbar.set_label('Temperature (°C)')

    plt.show()