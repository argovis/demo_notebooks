import random, pandas, matplotlib, math, copy
import scipy.interpolate
import numpy as np

import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.mpl.ticker as cticker

import dateutil

from argovisHelpers import helpers as avh
import xarray as xr

# get the correct route, given the collection name
def get_route(collection_name):
    
    if collection_name == 'drifters':
        return 'https://argovisbeta01.colorado.edu/dapi/'
    else:
        return 'https://argovis-api.colorado.edu/'

######## show list of variables for each collection in list
def show_variable_names_for_collections(collections_list,API_KEY,verbose=False):
    vars_list = []
    for icollection in collections_list:
        print('>>>>> '+icollection)
        try:
            bfr = avh.query(icollection+'/vocabulary', options={'parameter': 'data'}, verbose='true',apikey=API_KEY, apiroot=get_route(icollection)) 
            if verbose:
                print(bfr)
            vars_list.append(bfr)
        except:
            print('No data parameter for vocabulary query')
        try:
            bfr = avh.query(icollection+'/vocabulary', verbose='true',apikey=API_KEY, apiroot=get_route(icollection))
            if verbose:
                print(*bfr[0]['data'],sep=',')
            vars_list.append(bfr[0]['data'])
        except:
            print('Needs data parameter for vocabulary query')
    return vars_list

# interpolate profiles (this can be used also for grids/, as they are stored as profiles
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

# create an xarray from grids/ output (the output of grids_meta is also needed)
# varname indicates what variable to store in the xarray
def grids_to_xarray(grids,grids_meta,varname):
    data_list        = []
    data_list_lev    = []
    data_list_lon    = []
    data_list_lat    = []
    data_list_tstamp = []
    for x in grids:
        for ix,x_lev in enumerate(grids_meta[0]['levels']):
            if ix <= len(x['data'][x['data_info'][0].index(varname)])-1:
                data_list.append(x['data'][x['data_info'][0].index(varname)][ix])
            else:
                data_list.append(np.nan)
                
            data_list_lev.append(x_lev)
            data_list_lon.append(x['geolocation']['coordinates'][0])
            data_list_lat.append(x['geolocation']['coordinates'][1])
            data_list_tstamp.append(x['timestamp'])
            
    bfr_lon = np.array(data_list_lon)
    bfr_lon[bfr_lon<20] = bfr_lon[bfr_lon<20]+360

    data_list_lon = bfr_lon.tolist()

    data_dict = {'data': data_list,'latitude':data_list_lat, 'longitude': data_list_lon,'levels': data_list_lev,'timestamp': data_list_tstamp} 

    data_df = pandas.DataFrame(data_dict)   
    df_rows = pandas.DataFrame(data_df).set_index(["latitude", "longitude","levels","timestamp"])
    return xr.Dataset.from_dataframe(df_rows)

# area weigthed regional mean starting from an xarray
def xarray_regional_mean(dxr, form='area'):
    # given an xarray dataset <dxr> with latitudes and longitudes as dimensions,
    # calculate the horizontal average of all data variables, weighted by grid cell area
    weights = np.cos(np.deg2rad(dxr.latitude))
    weights.name = "weights"
    dxr_weighted = dxr.weighted(weights)
    
    if form =='area':
        return dxr_weighted.mean(("longitude", "latitude"))
    elif form == 'meridional':
        return dxr_weighted.mean(("latitude"))
    elif form == 'zonal':
        return dxr_weighted.mean(("longitude"))

# create a map from lists of longitudes and latitudes    
def map_lons_lats(lons,lats,dx=20,dy=20):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    lon_range = np.arange(np.floor(min(lons))-dx,np.ceil(max(lons))+dx,10)
    lat_range = np.arange(np.floor(min(lats))-dy,np.ceil(max(lats))+dy,10)
    
    if np.floor(min(lons))>=-180 and min(lon_range)<-180:
        lon_range=lon_range[lon_range>=-180]
    if np.floor(max(lons))<=180 and max(lon_range)>180:
        lon_range=lon_range[lon_range<=180]
    

    ax.set_extent([min(lon_range), max(lon_range), min(lat_range), max(lat_range)], crs=ccrs.PlateCarree())

    # Put a background image on for nice sea rendering.
    ax.stock_img()
    ax.plot(lons,lats,marker='.',linestyle='none',color='k')
    ax.set_xticks(lon_range, crs=ccrs.PlateCarree())
    lon_formatter = cticker.LongitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.set_yticks(lat_range, crs=ccrs.PlateCarree())
    lat_formatter = cticker.LatitudeFormatter()
    ax.yaxis.set_major_formatter(lat_formatter)

def list_values_for_parameter_to_api_query(selection_params,API_KEY):
    # list all the values for the parameter of interest
    vocab_params = {
        "parameter": selection_params['parameter_name'] 
    }
    for icollection in selection_params['collections']:
        print('>>>>>>>>>>>>>>> ' + icollection + ' <<<<<<<<<<<<<<<<<')
        param_vals = avh.query(icollection+'/vocabulary', options=vocab_params, apikey=API_KEY, apiroot=get_route(icollection))
        print(param_vals)
        # for woceline, we need to list the section_start_date as it is needed for querying the woceline occupation
        if selection_params['parameter_name'] == 'woceline':
            for i in param_vals:
                print('>>>>> ' + i + ' occupations:')
                i_vals = avh.query(icollection+'/meta', options={"woceline":i}, apikey=API_KEY, apiroot=get_route(icollection))
                for ii in i_vals:
                    for iii in ii['occupancies']:
                        print(iii['time_boundaries'])

# function to format the output of an api query, to make visualizing it easier
def format_api_output(api_output,selection_params,varname,index_collection,API_KEY=''):
# api_output is the output of an API query
# varname is the name of the variable we would like to save in the formatted output
# selection_params is a dictionary with info about the query (documentation to be included)
# index_collection is an integer pointing to the desired item in the list of collections within selection_params
    api_output_formatted = {}
    api_output_formatted['collection'] = selection_params['collections'][index_collection]
    api_output_formatted['varname']    = varname
    
    if api_output:
        # create a list with information for non gridded products
        if 'grids' not in selection_params['collections'][index_collection] and 'timeseries' not in selection_params['collections'][index_collection]:
            # if more than one variable is included in api_output, this will only store the one indicated in varname above
            api_output_formatted['_id']      =[x['_id'] for x in api_output]
            
            api_output_formatted['data']     =[x['data'][x['data_info'][0].index(varname)] for x in api_output]
            api_output_formatted['levels']   =[x['data'][x['data_info'][0].index(selection_params['varname_levels'][index_collection])] for x in api_output]
            api_output_formatted['timestamp']=[dateutil.parser.isoparse(x['timestamp']) for x in api_output]
            api_output_formatted['longitude']=[x['geolocation']['coordinates'][0] for x in api_output]
            api_output_formatted['latitude'] =[x['geolocation']['coordinates'][1] for x in api_output]

            api_output_formatted['data_units']=api_output[0]['data_info'][2][api_output[0]['data_info'][0].index(varname)][0]
            
            # if interp_levels are provided, then interpolate and create an xarray
            if 'interp_levels' in selection_params.keys():
                interpolated_profiles = []
                for idata in api_output:
                             interpolated_profiles.append(interpolate_profiles(profile=idata, levels_varname=selection_params['varname_levels'][index_collection], levels_new=selection_params['interp_levels']))
                        
                # if more than one variable is included in api_output, this will only store the one indicated in varname above
                d = [x['data'] for x in interpolated_profiles]
                d = [[level[varname] for level in x] for x in d]
                    
                # create xarray
                d_ind = np.array([list(range(1,len(api_output)+1))]*len(selection_params['interp_levels'])).transpose().tolist()
                d_lev = [selection_params['interp_levels']]*len(api_output)

                data_dict = {'data': np.array(d).flatten().tolist(),'levels': np.array(d_lev).flatten().tolist(),'index':np.array(d_ind).flatten().tolist()
                        } 
                data_df = pandas.DataFrame(data_dict)   
                df_rows = pandas.DataFrame(data_df).set_index(["levels","index"])
                xar     = xr.Dataset.from_dataframe(df_rows)
                xar.assign(longitude=(['index'],np.array(api_output_formatted['longitude'])))
                xar.assign(latitude=(['index'],np.array(api_output_formatted['latitude'])))
                xar.assign(timestamp=(['index'],np.array(api_output_formatted['timestamp'])))

                api_output_formatted['data_xarray'] = xar
        elif 'grids' in selection_params['collections'][index_collection]:
                # for grids, we only create data_xarray
                       
                for igm in api_output[0]['metadata']:
                    try:
                        # look in different meta_data documents until you find the one that has units for the variable of interest
                        grids_opt  = {
                            "id": igm
                            }
                        grids_meta = avh.query('grids/meta',options=grids_opt,verbose='true', apikey=API_KEY, apiroot=get_route(selection_params['collections'][index_collection]))
                        api_output_formatted['data_units'] = grids_meta[0]['data_info'][2][grids_meta[0]['data_info'][0].index(varname)][0] 
                        
                    except:
                        igm
#                         print(grids_meta['data_info'][0])
#                         print(grids_meta['data_info'][2])
                
                api_output_formatted['data_xarray'] = grids_to_xarray(api_output,grids_meta,varname)
                
    return api_output_formatted

# function to get a list of formatted api output
def get_api_output_formatted_list_1var_for_regions_and_timeranges(selection_params,API_KEY):
    # for each collection, region, and time range of interest, get api formatted output and store all in a list
    api_output_formatted_list = []

    for icl,icollection in enumerate(selection_params['collections']):
        for i,ireg in enumerate(selection_params['regions']):
            for istart,iend in zip(selection_params['startDate'],selection_params['endDate']):

                print('>>>>>>>>> '+icollection+' '+selection_params['varnames'][icl]+', '+selection_params['regions_tag'][i]+' '+istart[0:10]+' to '+iend[0:10])

                data_str = selection_params['varnames'][icl]+selection_params['varnames_qc'][icl]+selection_params['data_extra'][icl]
                iparam = {}
                iparam = {'data': data_str}
                if ireg:
                    iparam[selection_params['regions_type'][i]] = ireg
                if istart and 'glodap' not in icollection:
                    iparam['startDate'] = istart
                else:
                    istart = ''
                    if 'glodap' in icollection:
                        iparam['startDate'] = '1000-01-01T00:00:00.000Z'
                if iend and 'glodap' not in icollection:
                    iparam['endDate']   = iend
                else:
                    iend = ''
                    if 'glodap' in icollection:
                        iparam['endDate'] = '1000-01-02T00:00:00.000Z'
                api_output = avh.query(icollection, options=iparam, verbose='true',apikey=API_KEY, apiroot=get_route(icollection)) 
                
                api_output_formatted_all = {}
                for ivar in data_str.split(','):
                    if not ivar.isnumeric() and '~' not in ivar:
                        api_output_formatted = format_api_output(api_output=api_output,selection_params=selection_params,varname=selection_params['varnames'][icl],index_collection=icl,API_KEY=API_KEY) # please note that we specify varname as there may be more than one variable requested in other cases (in this specific function, we are focusing on comparing the same variable across datasets)

                        api_output_formatted_all[ivar]=api_output_formatted
                        # include some more info from selection_params
                        api_output_formatted_all[ivar]['region']   =ireg
                        api_output_formatted_all[ivar]['startDate']=istart
                        api_output_formatted_all[ivar]['endDate']  =iend

                        api_output_formatted_all[ivar]['region_type']=selection_params['regions_type'][i]
                        api_output_formatted_all[ivar]['region_tag']=selection_params['regions_tag'][i]
                        
                        api_output_formatted_all[ivar]['varname_title']=selection_params['varname_title']

                    
                    #print(api_output_formatted.keys())
                api_output_formatted_list.append(api_output_formatted_all)
    return api_output_formatted_list


def get_api_output_formatted_list_1var_for_parameter(selection_params,API_KEY):
    
    api_output_formatted_list = []
    
    for i,icollection in enumerate(selection_params['collections']):

        data_str = selection_params['varnames'][i]+selection_params['varnames_qc'][i]+selection_params['data_extra'][i]
        # we can only use the qc flags for datasets that support them
        if icollection not in ['argo','cchdo']:
            print('QC flags are used only for Argo and CCHDO')
            data_str = selection_params['varnames'][i]
            data_str_new = ''
            for ivar in data_str.split(','):
                if not ivar.isnumeric():
                    data_str_new = data_str_new+','+ivar
            data_str = data_str_new[1::]

        iparam = {}
        iparam = {'data': data_str}

        iparam[selection_params['parameter_name']] = selection_params['parameter'][i]

        if selection_params['parameter_name'] == 'woceline':
            # how the user assigns the start date should be improved
            iparam['section_start_date'] = selection_params['section_start_date'][i]

        # let's query data from the selected object
        api_output = avh.query(icollection, options=iparam, verbose='true',apikey=API_KEY, apiroot=get_route(icollection))

        # let's create api_output_formatted for each variable that was queried
        api_output_formatted_all = {}
        #api_output_formatted_all_var = []
        for ivar in data_str.split(','):
            if not ivar.isnumeric() and '~' not in ivar:
                #api_output_formatted_all_var.append(ivar)
                api_output_formatted_all[ivar] = format_api_output(api_output,selection_params,ivar,index_collection=i,API_KEY='')
                if 'region_tag' in selection_params.keys():
                    api_output_formatted_all[ivar]['region_tag']=selection_params['regions_tag'][i]
                else:
                    api_output_formatted_all[ivar]['region_tag']= ''
                    
                api_output_formatted_all[ivar]['startDate']     = min(api_output_formatted_all[ivar]['timestamp']).strftime("%Y-%m-%d")
                api_output_formatted_all[ivar]['endDate']       = max(api_output_formatted_all[ivar]['timestamp']).strftime("%Y-%m-%d")
                api_output_formatted_all[ivar]['varname_title'] = ivar[0].upper()+ivar[1::]
        api_output_formatted_list.append(api_output_formatted_all)
    return api_output_formatted_list


def api_output_formatted_list_1var_plot_lons_lats_map(api_output_formatted_list):
    # let's plot a map for the locations of point data (if the selected collections are not for point data, there will be no plot)
    for i_api_output_formatted_all in api_output_formatted_list:
        for ivar in i_api_output_formatted_all.keys():
            i_api_output_formatted =  i_api_output_formatted_all[ivar]
            #print(i_api_output_formatted.keys())
            if 'data' in i_api_output_formatted.keys():
                map_lons_lats(i_api_output_formatted['longitude'], i_api_output_formatted['latitude'],dx=20,dy=20)
                plt.title(i_api_output_formatted['collection']+' '+i_api_output_formatted['varname']+', '+i_api_output_formatted['region_tag']+'\n'+i_api_output_formatted['startDate'][0:10]+' to '+i_api_output_formatted['endDate'][0:10])

def api_output_formatted_list_1var_plot_profiles(api_output_formatted_list):
    # let's plot data for each of the point data (if the selected collections are not 
    # for point data, vertical profiles will be plotted for all the x,y points in the grid, at each time)
    for i_api_output_formatted_all in api_output_formatted_list:
        for ivar in i_api_output_formatted_all.keys():
            i_api_output_formatted =  i_api_output_formatted_all[ivar]
            #print(i_api_output_formatted.keys())
            if 'data_xarray' in i_api_output_formatted.keys():
                flag_xarray = np.logical_and('longitude' in list(i_api_output_formatted['data_xarray'].coords) and 'latitude' in list(i_api_output_formatted['data_xarray'].coords),'timestamp' in list(i_api_output_formatted['data_xarray'].coords))
            else:
                flag_xarray = False
            if 'data' in i_api_output_formatted.keys() or flag_xarray:
                plt.figure(figsize=(5,8))
                if 'data' in i_api_output_formatted.keys():
                    for i,idata in enumerate(i_api_output_formatted['data']):
                        plt.plot(idata,i_api_output_formatted['levels'][i],'k')
                elif flag_xarray:
                    # let's stack data_xarray so that one dimension is 'levels' and all the other dimensions are collapsed into a new dimension
                    xar_stacked = i_api_output_formatted['data_xarray'].stack(ind=("longitude","latitude","timestamp"))
                    plt.plot(xar_stacked['data'].values,i_api_output_formatted['data_xarray'].coords['levels'].values,'k')
                
                bfr_title = i_api_output_formatted['collection']+' '+i_api_output_formatted['varname']+', '+i_api_output_formatted['region_tag']
                if i_api_output_formatted['startDate'] and i_api_output_formatted['endDate']:
                    bfr_title = bfr_title + '\n'+i_api_output_formatted['startDate'][0:10]+' to '+i_api_output_formatted['endDate'][0:10]
                plt.title(bfr_title)
                plt.xticks(size=14);plt.yticks(size=14)
                plt.ylabel('Vertical level (m or dbar)',size=14)
                plt.xlabel(i_api_output_formatted['varname_title']+', '+ i_api_output_formatted['data_units'],size=14)
                plt.gca().invert_yaxis()

def api_output_formatted_list_1var_plot_map(api_output_formatted_list,ilev=0,itime=0):
    # plot the map for one timestep/level for each of the gridded products
    # ilev and itime are integers representing the index of intrest
    for i_api_output_formatted_all in api_output_formatted_list:
        for ivar in i_api_output_formatted_all.keys():
            i_api_output_formatted =  i_api_output_formatted_all[ivar]

            #print(i_api_output_formatted.keys())
            
            if 'data_xarray' in i_api_output_formatted.keys():
                if 'longitude' in list(i_api_output_formatted['data_xarray'].coords) and 'latitude' in list(i_api_output_formatted['data_xarray'].coords):
                    plt.figure()
                    i_api_output_formatted['data_xarray']['data'][:,:,ilev,itime].plot()
                    time_info = ''
                    if i_api_output_formatted['startDate'] and i_api_output_formatted['endDate']:
                        time_info = '\n'+i_api_output_formatted['startDate'][0:10]+' to '+i_api_output_formatted['endDate'][0:10]
                    plt.title(i_api_output_formatted['collection']+' '+i_api_output_formatted['varname']+','+i_api_output_formatted['region_tag']+time_info)

def api_output_formatted_list_1var_plot_horizontal_and_time_ave(api_output_formatted_list,colors):     
    # plot horizontal average using xarray objects
    # colors are a list with valid color names e.g. ['r' 'violet']
    
    figure_vars = []
    
    for i,i_api_output_formatted_all in enumerate(api_output_formatted_list):
        figure_vars.append(list(i_api_output_formatted_all.keys()))
        
    itest = len(figure_vars[0])
    for i in figure_vars:
        if len(i) != itest:
            print('Something is wrong, not all the items in api_output_formatted_list have the same number of variables')
            checkwhy
    
    for inum,inum_var in enumerate(figure_vars[0]):
        plt.figure(figsize=(6,8))
        leg = []
        for i,i_api_output_formatted_all in enumerate(api_output_formatted_list):
        
            i_api_output_formatted =  i_api_output_formatted_all[figure_vars[i][inum]]
    
            if 'data_xarray' in i_api_output_formatted.keys():
                if 'longitude' in list(i_api_output_formatted['data_xarray'].coords) and 'latitude' in list(i_api_output_formatted['data_xarray'].coords):
                    horiz_ave = xarray_regional_mean(i_api_output_formatted['data_xarray'], form='area').mean(dim='timestamp')
                elif 'index' in list(i_api_output_formatted['data_xarray'].coords):
                    horiz_ave = i_api_output_formatted['data_xarray'].mean(dim='index')
                #print(horiz_ave.sizes)
                horiz_ave['data'].plot(y='levels',color=colors[i],yincrease=False)
                plt.ylabel('Vertical level (m or dbar)',size=14)
                plt.xlabel(i_api_output_formatted['varname_title']+', '+ i_api_output_formatted['data_units'],size=14)
                plt.xticks(size=14);plt.yticks(size=14)
                time_info = ''
                if i_api_output_formatted['startDate'] and i_api_output_formatted['endDate']:
                    time_info = '\n'+i_api_output_formatted['startDate'][0:10]+' to '+i_api_output_formatted['endDate'][0:10]

                leg.append(i_api_output_formatted['collection']+' '+i_api_output_formatted['varname']+', '+i_api_output_formatted['region_tag']+time_info)
    plt.legend(leg)
    plt.show()