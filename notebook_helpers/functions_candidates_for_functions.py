import warnings
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
    """
    Flatten `compression='minimal'` query results into a DataFrame of profile
    locations and times (no measurements), e.g. for sampling statistics.

    nested_list : list of time slices, each a list of minimal rows
                  [id, longitude, latitude, timestamp, ...]

    Returns a DataFrame with lon (0–360), lat, time, year, month.
    """
    profiles = [p for time_slice in nested_list for p in time_slice]
    df = pd.DataFrame(profiles).iloc[:, [1, 2, 3]].copy()
    df.columns = ['lon', 'lat', 'time']
    df['lon']  = df['lon'] % 360
    df['time'] = pd.to_datetime(df['time'])
    df['year'], df['month'] = df['time'].dt.year, df['time'].dt.month
    return df


def bin_profiles(profiles, levels, centers, half_width, varname='temperature',
                 along='longitude', lat_range=None, lon_range=None):
    """
    Average profiles in bins along longitude or latitude.

    This is a plain bin mean, NOT a mapped/gridded product: there is no
    weighting by distance or time, no mapping error, and no correction for
    uneven sampling in space, season, or year. Each bin is simply the mean of
    whichever profiles fall in it, so results reflect where and when floats
    happened to sample. For gridded fields use a mapped product (e.g. RG09,
    LocalGP).

    Bins are [center - half_width, center + half_width) and may overlap
    (e.g. 0.5°-wide bins every 0.25°, as in Karnauskas & Giglio 2022).

    Parameters
    ----------
    profiles   : list of Profile – already interpolated onto `levels`
                 (e.g. with ava.interpolate_all)
    levels     : array-like      – the levels the profiles were interpolated onto
    centers    : array-like      – bin centers along `along`; longitudes in 0–360
    half_width : float           – half the bin width, in degrees
    varname    : str             – variable to average (default 'temperature')
    along      : str             – 'longitude' or 'latitude'
    lat_range  : (min, max)      – optional; keep only profiles in this latitude range
    lon_range  : (min, max)      – optional; keep only profiles in this longitude
                                   range (0–360)

    Returns
    -------
    xr.Dataset with `<varname>` (along × level) and `n_profiles` (along).
    Empty bins are NaN, with n_profiles = 0.
    """
    if along not in ('longitude', 'latitude'):
        raise ValueError("along must be 'longitude' or 'latitude'")

    levels = np.asarray(levels, dtype=float)
    centers = np.asarray(centers, dtype=float)
    lon = np.array([p.longitude % 360 for p in profiles])
    lat = np.array([p.latitude for p in profiles])

    data = np.full((len(profiles), len(levels)), np.nan)
    for i, p in enumerate(profiles):
        v = p.getvar(varname)
        if v is None:
            continue
        if len(v) != len(levels):
            raise ValueError(f"profile {p.id}: '{varname}' has {len(v)} levels, expected "
                             f"{len(levels)} — interpolate profiles onto `levels` first")
        data[i] = v

    keep = ~np.all(np.isnan(data), axis=1)
    if lat_range is not None:
        keep &= (lat >= lat_range[0]) & (lat <= lat_range[1])
    if lon_range is not None:
        keep &= (lon >= lon_range[0]) & (lon <= lon_range[1])
    coord = lon if along == 'longitude' else lat

    mean = np.full((len(centers), len(levels)), np.nan)
    count = np.zeros(len(centers), dtype=int)
    for j, c in enumerate(centers):
        in_bin = keep & (coord >= c - half_width) & (coord < c + half_width)
        count[j] = in_bin.sum()
        if count[j]:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', RuntimeWarning)  # all-NaN levels
                mean[j] = np.nanmean(data[in_bin], axis=0)

    return xr.Dataset(
        {varname: ((along, 'level'), mean), 'n_profiles': (along, count)},
        coords={along: centers, 'level': levels},
        attrs={'method': f'bin mean of individual profiles, bin half-width {half_width}° '
                         '— not a mapped product'},
    )


# ── Helper functions: seasonal-cycle and trend removal ──────────────────────

def remove_seasonal_and_trend(da, time_dim='timestamp', seasonal='climatology',
                              n_harmonics=2, trend=None, fit_period=None,
                              smooth_window=None, return_components=False):
    """
    Remove the seasonal cycle and (optionally) a trend from a time series or a
    time × space field, returning the anomaly.

    The climatology (time mean + seasonal cycle) and the trend are fit together
    by least squares, separately at every grid point, NaN-safe. Fitting them
    jointly avoids the bias of doing them one after the other when the record
    is short or starts/ends in a particular season.

    Parameters
    ----------
    da            : xr.DataArray – series or field with a datetime `time_dim`
    time_dim      : str          – name of the time dimension (default 'timestamp')
    seasonal      : str or None  – 'climatology': one mean per calendar month
                                   'harmonic'   : annual harmonics (smoother; better
                                                  for short or gappy records)
                                   None         : remove the time mean only
    n_harmonics   : int          – number of harmonics if seasonal='harmonic'
                                   (1 = annual, 2 = annual + semi-annual, ...)
    trend         : str or None  – None, 'linear', or 'quadratic'
    fit_period    : (start, end) – optional; estimate climatology and trend from
                                   this period only (e.g. ('1993-01-01', '2012-12-31')),
                                   then remove them from the whole record
    smooth_window : int or None  – optional centred running mean of the anomaly,
                                   in time steps
    return_components : bool     – if True, return a Dataset with `anomaly`,
                                   `climatology` and `trend`

    Returns
    -------
    xr.DataArray anomaly (same shape as `da`), or an xr.Dataset of components.
    With seasonal='climatology' and trend=None this matches subtracting the
    monthly climatology.
    """
    if seasonal not in ('climatology', 'harmonic', None):
        raise ValueError("seasonal must be 'climatology', 'harmonic' or None")
    if trend not in (None, 'linear', 'quadratic'):
        raise ValueError("trend must be None, 'linear' or 'quadratic'")

    times = da[time_dim].to_index()
    if fit_period is None:
        in_fit = np.ones(len(times), dtype=bool)
    else:
        start, end = pd.Timestamp(fit_period[0]), pd.Timestamp(fit_period[1])
        in_fit = np.asarray((times >= start) & (times <= end))
        if not in_fit.any():
            raise ValueError(f'no time steps inside fit_period {fit_period}')

    # time in years, centred on the fit period (keeps the quadratic well-conditioned)
    t_fit = times[in_fit]
    t0 = t_fit.min() + (t_fit.max() - t_fit.min()) / 2
    t_yr = np.asarray((times - t0) / pd.Timedelta(days=365.25), dtype=float)

    # design matrix: climatology columns first, then trend columns
    if seasonal == 'climatology':
        clim_cols = [(times.month == m).astype(float) for m in range(1, 13)]
    else:
        clim_cols = [np.ones_like(t_yr)]
        if seasonal == 'harmonic':
            for k in range(1, n_harmonics + 1):
                clim_cols += [np.cos(2 * np.pi * k * t_yr), np.sin(2 * np.pi * k * t_yr)]
    degree = {None: 0, 'linear': 1, 'quadratic': 2}[trend]
    trend_cols = [t_yr ** d for d in range(1, degree + 1)]
    X = np.column_stack(clim_cols + trend_cols)
    n_clim = len(clim_cols)

    def _components(coef, used):
        # coef: (n_columns, n_series); rows of unused columns are ignored
        clim = X[:, :n_clim] @ np.where(used[:n_clim, None], coef[:n_clim], 0)
        clim[X[:, :n_clim][:, ~used[:n_clim]].any(axis=1)] = np.nan  # unsampled months
        return clim, X[:, n_clim:] @ coef[n_clim:]

    # one row per grid point, time last
    da_t = da.astype(float).transpose(..., time_dim)
    y = da_t.values.reshape(-1, len(times))
    clim = np.full(y.shape, np.nan)
    trnd = np.full(y.shape, np.nan)

    finite = np.isfinite(y[:, in_fit])
    complete = finite.all(axis=1)
    if complete.any():
        # all grid points without gaps share the same fit: solve them in one call
        used = X[in_fit].any(axis=0)
        if in_fit.sum() > used.sum():
            coef = np.zeros((X.shape[1], complete.sum()))
            coef[used] = np.linalg.lstsq(X[in_fit][:, used], y[complete][:, in_fit].T, rcond=None)[0]
            c, t = _components(coef, used)
            clim[complete], trnd[complete] = c.T, t.T
    for i in np.flatnonzero(~complete & finite.any(axis=1)):
        # grid points with gaps: fit each one on its own valid samples
        ok = in_fit & np.isfinite(y[i])
        used = X[ok].any(axis=0)
        if ok.sum() <= used.sum():
            continue                                  # not enough points to fit
        coef = np.zeros((X.shape[1], 1))
        coef[used, 0] = np.linalg.lstsq(X[ok][:, used], y[i, ok], rcond=None)[0]
        c, t = _components(coef, used)
        clim[i], trnd[i] = c[:, 0], t[:, 0]

    clim = da_t.copy(data=clim.reshape(da_t.shape)).transpose(*da.dims)
    trnd = da_t.copy(data=trnd.reshape(da_t.shape)).transpose(*da.dims)

    anomaly = da - clim - trnd
    if smooth_window is not None:
        anomaly = anomaly.rolling({time_dim: smooth_window}, center=True).mean()

    if return_components:
        return xr.Dataset({'anomaly': anomaly, 'climatology': clim, 'trend': trnd})
    return anomaly
