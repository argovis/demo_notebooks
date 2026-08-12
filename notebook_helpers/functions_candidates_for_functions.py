import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from argovisHelpers import helpers as avh
import math
import xarray as xr
from dateutil import parser
from notebook_helpers.functions import plot_maps, compare_profiles
from argovisHelpers import analysis as ava
import gsw

def _arr(val):
    """Safely convert Profile.getvar() output to a plain float ndarray."""
    return val.filled(np.nan) if hasattr(val, 'filled') else np.asarray(val, dtype=float)


def _safe_interp(profiles, levels):
    """Interpolate profiles, skipping any with fewer than 2 valid levels."""
    out = []
    for p in profiles:
        try:
            out.append(ava.interpolate_all(p, levels))
        except ValueError:
            pass
    return out


def parse_profiles(nested_list):
    profiles = [p for time_slice in nested_list for p in time_slice]
    df = pd.DataFrame(profiles)
    df = df[['geolocation', 'timestamp', 'date']].copy() if 'date' in df.columns else df.iloc[:, [1, 2, 3]].copy()
    df.columns = ['lon', 'lat', 'time']
    df['lon']  = df['lon'] % 360
    df['time'] = pd.to_datetime(df['time'])
    df['year'], df['month'] = df['time'].dt.year, df['time'].dt.month
    return df


def build_grid(profiles, x_centers, x_half_width,
               y_centers=None, y_half_width=None,
               lat_min=None, lat_max=None,
               lon_min=None, lon_max=None,
               use_rawdata=False,
               profile_levels=None):
    """
    Bin profiles by longitude (and optionally latitude) and compute mean temperature.

    profile_levels : array-like, required when use_rawdata=True.
        Profiles are assumed already interpolated onto these levels.
    """
    def get_coords(p):
        if use_rawdata:
            lon = p.rawdata['geolocation']['coordinates'][0]
            lat = p.rawdata['geolocation']['coordinates'][1]
        else:
            lon = p['geolocation']['coordinates'][0]
            lat = p['geolocation']['coordinates'][1]
        return lon, lat

    def get_temperature(p):
        if use_rawdata:
            if profile_levels is not None:
                # Already interpolated onto these levels — just extract
                return np.array(_arr(p.getvar('temperature')), dtype=float)
            else:
                raise ValueError(
                    "build_grid: profile_levels is required when use_rawdata=True. "
                    "Pass the levels your profiles were interpolated onto."
                )
        else:
            raise ValueError(
                "build_grid: use_rawdata=False is no longer supported. "
                "Pass use_rawdata=True with profile_levels."
            )

    def get_levels():
        if use_rawdata and profile_levels is not None:
            return profile_levels
        else:
            raise ValueError(
                "build_grid: profile_levels is required when use_rawdata=True."
            )

    def get_bin_mean(x_center, y_center=None):
        x1, x2 = x_center - x_half_width, x_center + x_half_width
        y1 = y_center - y_half_width if y_center is not None else None
        y2 = y_center + y_half_width if y_center is not None else None
        x_is_lat = lon_min is not None
        temps = []
        for p in profiles:
            lon, lat = get_coords(p)
            if lon is None or lat is None:
                continue
            if lat_min is not None and not (lat_min <= lat <= lat_max):
                continue
            if lon_min is not None and not (lon_min <= lon % 360 <= lon_max):
                continue
            if x_is_lat:
                if not (x1 <= lat < x2):
                    continue
            else:
                if not (x1 <= lon % 360 < x2):
                    continue
            if y1 is not None and not (y1 <= lat < y2):
                continue
            try:
                t = get_temperature(p)
            except Exception as e:
                print(f"Error for profile at lon={lon:.2f}, lat={lat:.2f}: {e}")
                continue
            if t is None or np.all(np.isnan(t)):
                continue
            temps.append(np.array(t, dtype=float))
        if not temps:
            return None
        mean_temp = np.nanmean(np.vstack(temps), axis=0)
        lvls = get_levels()
        return xr.Dataset({'temperature': xr.DataArray(mean_temp, dims='level',
                                                        coords={'level': lvls})})

    def concat_grid(grid, centers, dim_name):
        non_empty = [ds for ds in grid.values() if ds is not None]
        if not non_empty:
            raise ValueError("build_grid: no profiles found in any bin.")
        template = non_empty[0]
        nan_fill = xr.full_like(template, fill_value=np.nan)
        return xr.concat(
            [grid[i] if grid[i] is not None else nan_fill
             for i in range(len(centers))],
            dim=xr.DataArray(centers, dims=dim_name)
        )

    if y_centers is not None:
        grid = {
            (i, j): get_bin_mean(xc, yc)
            for i, xc in enumerate(x_centers)
            for j, yc in enumerate(y_centers)
        }
        non_empty = [ds for ds in grid.values() if ds is not None]
        if not non_empty:
            raise ValueError("build_grid: no profiles found in any bin.")
        template = non_empty[0]
        nan_fill = xr.full_like(template, fill_value=np.nan)
        combined = xr.concat(
            [
                xr.concat(
                    [grid.get((i, j), nan_fill) if grid.get((i, j)) is not None else nan_fill
                     for j in range(len(y_centers))],
                    dim=xr.DataArray(y_centers, dims='latitude')
                )
                for i in range(len(x_centers))
            ],
            dim=xr.DataArray(x_centers, dims='longitude')
        )
    else:
        grid     = {i: get_bin_mean(xc) for i, xc in enumerate(x_centers)}
        combined = concat_grid(grid, x_centers, 'longitude')

    return grid, combined


def compute_n2_profile(p):
    """Compute N² profile using gsw.Nsquared instead of the custom helper."""
    SA  = _arr(p.getvar('absolute_salinity'))
    CT  = _arr(p.getvar('conservative_temperature'))
    pres = _arr(p.getvar('pressure'))

    valid = ~np.isnan(SA) & ~np.isnan(CT) & ~np.isnan(pres)
    SA_v, CT_v, pres_v = SA[valid], CT[valid], pres[valid]

    order = np.argsort(pres_v)
    SA_v, CT_v, pres_v = SA_v[order], CT_v[order], pres_v[order]

    N2, p_mid = gsw.Nsquared(SA_v, CT_v, pres_v, lat=p.latitude)

    good = ~np.isnan(N2)
    return p_mid[good], N2[good]


# ── Helper functions: spatial average and seasonal-cycle removal ────────────

def spatial_average(da, lat_band, lat_dim='latitude', lon_dim=None):
    """
    Return the cos-latitude-weighted mean of *da* over *lat_band* = (lat_min, lat_max).
    Only a latitude range can be specified — there is no equivalent lon_band
    argument. If *lon_dim* is given, longitude is also fully collapsed to a
    single averaged value (scalar output per time step); otherwise longitude
    is left untouched at its full queried resolution and only latitude is
    collapsed (Hovmöller output).

    Requires a gridded field: *da* must be an xr.DataArray on a regular
    lat/lon grid (e.g. as returned by tidy_grid), not a list of individual
    profiles.

    Parameters
    ----------
    da       : xr.DataArray  – gridded input field with a latitude coordinate
    lat_band : tuple         – the (lat_min, lat_max) *range of values* to average over
    lat_dim  : str           – the *name* of the latitude coordinate in da (default
                                'latitude'), not a range — lets this work even if a
                                dataset labels its latitude coordinate differently
    lon_dim  : str or None   – the *name* of the longitude coordinate to also collapse,
                                if given (same "coordinate name" role as lat_dim, not
                                a range — there is no lon_band equivalent)

    Returns
    -------
    xr.DataArray with the specified dimensions collapsed.
    """
    sub     = da.sel({lat_dim: slice(*lat_band)})
    weights = np.cos(np.deg2rad(sub[lat_dim]))
    weights.name = 'weights'
    dims = [lat_dim] if lon_dim is None else [lat_dim, lon_dim]
    return sub.weighted(weights).mean(dim=dims)


def remove_seasonal_cycle(da, time_dim='timestamp', smooth_window=None):
    """
    Remove the monthly climatological seasonal cycle from *da*.

    Steps
    -----
    1. Build the monthly climatology by averaging each calendar month.
    2. Subtract it from the original time series (anomaly).
    3. Optionally apply a running-mean smoother of width *smooth_window*
       (number of time steps, centred).

    Parameters
    ----------
    da            : xr.DataArray  – time series or time × space field
    time_dim      : str           – name of the time dimension (default 'timestamp')
    smooth_window : int or None   – rolling-mean window; None = no smoothing

    Returns
    -------
    xr.DataArray – anomaly (same shape as *da*, or smoothed if requested).
    """
    clim  = da.groupby(f'{time_dim}.month').mean(time_dim)
    anom  = da.groupby(f'{time_dim}.month') - clim
    if smooth_window is not None:
        anom = anom.rolling({time_dim: smooth_window}, center=True).mean()
    return anom