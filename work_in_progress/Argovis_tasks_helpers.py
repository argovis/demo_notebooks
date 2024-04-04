import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import random, numpy, pandas, matplotlib, math, copy
import scipy.interpolate

from argovisHelpers import helpers as avh
import xarray as xr

def get_route(collection_name):
    
    if collection_name == 'drifters':
        return 'https://argovisbeta01.colorado.edu/dapi/'
    else:
        return 'https://argovis-api.colorado.edu/'
    
def interpolate_profiles(profile, levels_varname, levels_new):
    # given a <profile>, a string for the name <levels_varname> of the variable containing the levels, and a list of desired values for vertical levels <levels_new>,
    # return a profile with profile.data levels at the desired vertical levels, with all available data interpolated to match
    # drop all QC and note `data_interpolated` in profile.data_warnings
    # please note that all the profile points will be used: if a selection by QC is needed, that should be done as part of the API query
    
    data_names = [levels_varname]
    interpolated_data = [levels_new]
    for key in profile['data_info'][0]:
        if '_argoqc' not in key and '_woceqc' not in key and key!=levels_varname:
            # create a list of data values and associated pressure and only keep the levels where both are valid values
            lvl = avh.data_inflate(profile)
            finites = [(level[levels_varname], level[key]) for level in lvl if level[levels_varname] is not None and level[key] is not None and not math.isnan(level[levels_varname]) and not math.isnan(level[key])]
            
            # store some values
            levels4interpolation = [x[0] for x in finites]
            data = [x[1] for x in finites]
            data_names.append(key)
            
            # interpolate avoiding extrapolation
            interpolated_data.append(scipy.interpolate.PchipInterpolator(levels4interpolation, data, extrapolate=False)(levels_new).tolist())
    
    interpolated_levels = list(zip(*interpolated_data))
    data = [{data_names[i]:d[i] for i in range(len(data_names))} for d in interpolated_levels]
    interpolated_profile = copy.deepcopy(profile) # don't mutate the original
    interpolated_profile['data'] = data
    if 'data_warnings' in interpolated_profile:
        interpolated_profile['data_warnings'].append('data_interpolated')
    else:
        interpolated_profile['data_warnings'] = ['data_interpolated']
    return interpolated_profile

def grids_to_xarray(grids,grids_meta):
    data_list        = []
    data_list_lev    = []
    data_list_lon    = []
    data_list_lat    = []
    data_list_tstamp = []
    for x in grids:
        for ix,x_lev in enumerate(grids_meta[0]['levels']):
            if ix <= len(x['data'][0])-1:
                data_list.append(x['data'][0][ix])
            else:
                data_list.append(numpy.nan)
                
            data_list_lev.append(x_lev)
            data_list_lon.append(x['geolocation']['coordinates'][0])
            data_list_lat.append(x['geolocation']['coordinates'][1])
            data_list_tstamp.append(x['timestamp'])
            
    bfr_lon = numpy.array(data_list_lon)
    bfr_lon[bfr_lon<20] = bfr_lon[bfr_lon<20]+360

    data_list_lon = bfr_lon.tolist()

    data_dict = {'data': data_list,'latitude':data_list_lat, 'longitude': data_list_lon,'levels': data_list_lev,'timestamp': data_list_tstamp} 

    data_df = pandas.DataFrame(data_dict)   
    df_rows = pandas.DataFrame(data_df).set_index(["latitude", "longitude","levels","timestamp"])
    return xr.Dataset.from_dataframe(df_rows)
    
def xarray_regional_mean(dxr, form='area'):
    # given an xarray dataset <dxr> with latitudes and longitudes as dimensions,
    # calculate the mean of all data variables, weighted by grid cell area
    weights = numpy.cos(numpy.deg2rad(dxr.latitude))
    weights.name = "weights"
    dxr_weighted = dxr.weighted(weights)
    
    if form =='area':
        return dxr_weighted.mean(("longitude", "latitude"))
    elif form == 'meridional':
        return dxr_weighted.mean(("latitude"))
    elif form == 'zonal':
        return dxr_weighted.mean(("longitude"))
    