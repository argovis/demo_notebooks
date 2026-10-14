import os
import urllib.request
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

_AR6_REPO = 'https://raw.githubusercontent.com/IPCC-WG1/Chapter-9/main/Plotting_code_and_data'
AR6_URLS = {
    'Fig9-3c_data.nc': f'{_AR6_REPO}/Fig9_03_SST/Plotted_Data/Fig9-3c_data.nc',
    'Fig9-6b_data.nc': f'{_AR6_REPO}/Fig9_06_OHC/Plotted_Data/Fig9-6b_data.nc',
    'Fig9-6e_data.nc': f'{_AR6_REPO}/Fig9_06_OHC/Plotted_Data/Fig9-6e_data.nc',
}


def _ensure_ar6_file(path):
    """Download an AR6 Chapter 9 data file from the IPCC GitHub repo if `path` is missing."""
    if os.path.exists(path):
        return
    name = os.path.basename(path)
    if name not in AR6_URLS:
        raise FileNotFoundError(f"{path} not found, and no download URL is known for '{name}'.")
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    try:
        urllib.request.urlretrieve(AR6_URLS[name], path)
    except Exception as e:
        if os.path.exists(path):
            os.remove(path)
        raise RuntimeError(
            f"Could not download {name} from {AR6_URLS[name]} ({e}). "
            f"Download it manually and save it as {path}."
        ) from e


def load_ar6_grid(path, varname, lon_name='Longitude', lat_name='Latitude'):
    """Load an IPCC AR6 WG1 Chapter 9 "Plotted_Data" netCDF file
    (github.com/IPCC-WG1/Chapter-9), downloading it first if it is not already
    on disk, and return the named variable as a
    (latitude x longitude) DataArray on a -180..180 longitude convention,
    sorted for plotting — matching the convention `plot_ohc_trend_map` and
    our own Argovis-derived trend maps use.

    Parameters
    ----------
    path      : str  path to a Fig9-*_data.nc file
    varname   : str  the data variable to extract, e.g. 'SST_ChangeRate', 'OHC_Trend'
    lon_name, lat_name : str  names of the longitude/latitude coordinates in the file
                        (AR6 files use 'Longitude'/'Latitude', not lowercase)
    """
    _ensure_ar6_file(path)
    ds = xr.open_dataset(path)
    da = ds[varname].rename({lon_name: 'longitude', lat_name: 'latitude'})
    lon = da['longitude'].values
    da = da.assign_coords(longitude=np.where(lon > 180, lon - 360, lon))
    da = da.drop_duplicates('longitude').sortby('longitude').sortby('latitude')
    return da


def annual_mean_dropping_incomplete(da, time_dim='timestamp', min_frac=0.75):
    """Annual mean of `da`, dropping any calendar year with fewer samples than
    `min_frac` of the typical (median) per-year sample count. Guards against a
    partial year — e.g. right at a product's start/end cutoff — silently
    biasing the annual mean toward whichever season happens to be sampled.
    """
    counts = da.resample({time_dim: '1YS'}).count()
    means = da.resample({time_dim: '1YS'}).mean()
    valid_counts = counts.values[counts.values > 0]
    typical = np.median(valid_counts) if len(valid_counts) else 0
    return means.where(counts >= min_frac * typical, drop=True)


def plot_multi_region_timeseries(series, ylabel, title, colors=None):
    """Single-panel comparison of one time series per region, overlaid.

    Parameters
    ----------
    series : dict  {region_label: xr.DataArray(time)} — e.g. annual SST or OHC per region
    ylabel : str
    title  : str
    """
    if colors is None:
        colors = list(plt.cm.tab10.colors)
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    for i, (name, da) in enumerate(series.items()):
        time_dim = da.dims[0]
        ax.plot(da[time_dim], da, color=colors[i % len(colors)], marker='o', linewidth=2, label=name)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(loc='best', fontsize=9)
    ax.grid(True, linestyle=':', linewidth=0.5)
    plt.tight_layout()
    plt.show()


def plot_eei_components(ds, components, labels=None, colors=None,
                        title='Where does the excess energy go?',
                        ylabel='Cumulative heat content anomaly (ZJ)'):
    """Stacked-area chart of Earth-system heat-content components from the
    GCOS Earth Heat Inventory — the same atmosphere/land/cryosphere/ocean
    breakdown as Fig. 8 in von Schuckmann et al. (2023).

    Parameters
    ----------
    ds         : xr.Dataset  opened GCOS_EHI_*.nc file
    components : list of str  data variable names to stack (non-overlapping!),
                 e.g. ['atmospheric_heat_content', 'continental_heat_content',
                       'cryosphere_heat_content', 'ocean_heat_content_0-2000m',
                       'ocean_heat_content_below_2000m']
    labels     : list of str, optional
    """
    if labels is None:
        labels = components
    if colors is None:
        colors = list(plt.cm.tab10.colors)
    fig, ax = plt.subplots(figsize=(10, 6))
    values = [np.nan_to_num(ds[c].values) for c in components]
    ax.stackplot(ds.time, *values, labels=labels, colors=colors[:len(components)], alpha=0.85)
    ax.axhline(0, color='gray', linewidth=0.5)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(loc='upper left', fontsize=9)
    ax.grid(True, linestyle=':', linewidth=0.5)
    plt.tight_layout()
    plt.show()
