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

