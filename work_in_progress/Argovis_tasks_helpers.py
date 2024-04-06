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
def show_variable_names_for_collections(collections_list,API_KEY):
    for icollection in collections_list:
        try:
            print(avh.query(icollection+'/vocabulary', options={'parameter': 'data'}, verbose='true',apikey=API_KEY, apiroot=get_route(icollection)) )
        except:
            print('No data parameter for vocabulary query')
        try:
            bfr = avh.query(icollection+'/vocabulary', verbose='true',apikey=API_KEY, apiroot=get_route(icollection))
            print(*bfr[0]['data'],sep=',')
        except:
            print('Needs data parameter for vocabulary query')
    
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
                grids_opt  = {
                            "id": api_output[0]['metadata'][0]
                            }
                grids_meta = avh.query('grids/meta', options=grids_opt, verbose='true',apikey=API_KEY, apiroot=get_route(selection_params['collections'][index_collection]))
                       
                api_output_formatted['data_units'] = grids_meta[0]['data_info'][2][grids_meta[0]['data_info'][0].index(varname)][0] 
                
                api_output_formatted['data_xarray'] = grids_to_xarray(api_output,grids_meta,varname)
                
    return api_output_formatted

# function to get a list of formatted api output
def get_api_output_formatted_list_1var(selection_params,API_KEY):

    api_output_formatted_list = []

    for icl,icollection in enumerate(selection_params['collections']):
        for i,ireg in enumerate(selection_params['regions']):
            for istart,iend in zip(selection_params['startDate'],selection_params['endDate']):

                print('>>>>>>>>> '+icollection+' '+selection_params['varnames'][icl]+', '+selection_params['regions_tag'][i]+' '+istart[0:10]+' to '+iend[0:10])

                iparam = {}
                iparam = {'data': selection_params['varnames'][icl]+selection_params['varnames_qc'][icl]+selection_params['data_extra'][icl]}
                if ireg:
                    iparam[selection_params['regions_type'][i]] = ireg
                if istart and 'glodap' not in icollection:
                    iparam['startDate'] = istart
                else:
                    istart = ''
                if iend and 'glodap' not in icollection:
                    iparam['endDate']   = iend
                else:
                    iend = ''
                api_output = avh.query(icollection, options=iparam, verbose='true',apikey=API_KEY, apiroot=get_route(icollection)) 
                api_output_formatted = format_api_output(api_output=api_output,selection_params=selection_params,varname=selection_params['varnames'][icl],index_collection=icl,API_KEY=API_KEY) # please note that we specify varname as there may be more than one variable requested in other cases (in this specific function, we are focusing on comparing the same variable across datasets)

                # include some more info from selection_params
                api_output_formatted['region']   =ireg
                api_output_formatted['startDate']=istart
                api_output_formatted['endDate']  =iend

                api_output_formatted['region_type']=selection_params['regions_type'][i]
                api_output_formatted['region_tag']=selection_params['regions_tag'][i]

                api_output_formatted['varname_title']=selection_params['varname_title']

                #print(api_output_formatted.keys())
                api_output_formatted_list.append(api_output_formatted)
    return api_output_formatted_list
    

def api_output_formatted_list_1var_plot_lons_lats_map(api_output_formatted_list):
    # let's plot a map for the locations of point data (if the selected collections are not for point data, there will be no plot)
    for i_api_output_formatted in api_output_formatted_list:
        #print(i_api_output_formatted.keys())
        if 'data' in i_api_output_formatted.keys():
            map_lons_lats(i_api_output_formatted['longitude'], i_api_output_formatted['latitude'],dx=20,dy=20)
            plt.title(i_api_output_formatted['collection']+' '+i_api_output_formatted['varname']+', '+i_api_output_formatted['region_tag']+'\n'+i_api_output_formatted['startDate'][0:10]+' to '+i_api_output_formatted['endDate'][0:10])

def api_output_formatted_list_1var_plot_profiles(api_output_formatted_list):
    # let's plot data for each of the point data (if the selected collections are not for point data, there will be no plot)
    for i_api_output_formatted in api_output_formatted_list:
        #print(i_api_output_formatted.keys())
        if 'data' in i_api_output_formatted.keys():
            plt.figure(figsize=(5,8))
            for i,idata in enumerate(i_api_output_formatted['data']):
                plt.plot(idata,i_api_output_formatted['levels'][i],'k')
            plt.title(i_api_output_formatted['collection']+' '+i_api_output_formatted['varname']+', '+i_api_output_formatted['region_tag']+'\n'+i_api_output_formatted['startDate'][0:10]+' to '+i_api_output_formatted['endDate'][0:10])
            plt.xticks(size=14);plt.yticks(size=14)
            plt.ylabel('Vertical level (m or dbar)',size=14)
            plt.xlabel(i_api_output_formatted['varname_title']+', '+ i_api_output_formatted['data_units'],size=14)
            plt.gca().invert_yaxis()

def api_output_formatted_list_1var_plot_map(api_output_formatted_list,ilev=0,itime=0):
    # plot the map for one timestep/level for each of the gridded products
    # ilev and itime are integers representing the index of intrest
    for i_api_output_formatted in api_output_formatted_list:
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
    plt.figure(figsize=(5,8))
    leg = []
    for i,i_api_output_formatted in enumerate(api_output_formatted_list):
        #print(i_api_output_formatted.keys())
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

# get profiles and perform horizontal ave, given a list of regions of interest, collections, and other info
# def get_profiles_in_regions_and_horiz_ave(collections,varname,varname_qc,varname_levels,interp_levels,regions_list_source,regions_list_source_type,regions_list_source_tags,startDate,endDate,API_KEY=''):
#     # input:
#     # collections: list of collections of interest, e.g. ['argo', 'grids/glodap'] or ['argo', 'grids/rg09']. Please note, these collections should include profiles (e.g. grids/ are stored as profiles, this is different for BSOSE)
#     # varname: list with the names of the variable of interest for each of the collections, e.g. ['doxy', 'oxygen']
#     # varname_qc: list with the qc values of the variable of interest for each of the collections, e.g. [',1', '']
#     # varname_levels: list with the names of the variable with vertical levels for each of the collections (for the grids route, this should be ''), e.g. ['pressure','']
#     # interp_levels: list with the levels to use for the vertical interpolation, e.g. list(range(10,2001))[0::20]
#     # regions_list_source: a list of regions to be used in the 'box' or 'polygon' searches (if empty [], no region is specified in the API query), e.g. for 'box' one possibility (including two regions) is [ [[-179.5,45.5],[-170.5,50.5]], [[-50,45],[-40,50]] ]
#     # regions_list_source_type: 'box' or 'polygon' depending on what regions_list is
#     # API_KEY: string with API_KEY or ''
#     #
#     # notes: the time range does not apply to glodap as for that product only the time mean is available
#     output_vars = ['regions_list_data_raw', 'regions_list_data_raw_xarray', 'regions_list_data_vert_interp', 'regions_list_data_horiz_ave', 'regions_list_data_horiz_ave_levels', 'regions_list_data_time', 'regions_list', 'regions_list_collections', 'regions_list_tags']
    
#     for ivar in output_vars:
#         globals()[ivar] = []
    
#     ############ Get the data for each of the regions of interest, for each collection
#     for iiireg,ireg in enumerate(regions_list_source):
#         for icol_ind,icollection in enumerate(collections):
#             # print and store the region and collection for each item in the output lists
#             print('>>>>>>> Region '+str(ireg)+' , '+icollection+' collection')
#             regions_list.append(ireg)
#             regions_list_collections.append(icollection)
#             regions_list_tags.append(regions_list_source_tags[iiireg])
#             ###### get profiles of interest using Argovis API (query based on qc if requested above for a collection)
#             iparam = {'data': varname[icol_ind]+varname_qc[icol_ind]}
            
#             # if provided, search in a region of interest
#             if ireg and (regions_list_source_type=='box' or regions_list_source_type=='polygon'):
#                 iparam[regions_list_source_type] = ireg
                
#             # the time range does not apply to glodap as only the time mean is available for glodap
#             if 'glodap' in icollection:
#                 print('For the glodap product, only the time mean is available')
#             else:
#                 iparam['startDate'] = startDate
#                 iparam['endDate']   = endDate
#             api_output = avh.query(icollection, options=iparam, verbose='true',apikey=API_KEY, apiroot=get_route(icollection)) 

#             ###### store data as is and interpolated
#             if api_output: # len(api_output)>0
#                 ## interpolate profiles (if not from a grid)   
#                 interpolated_profiles      = []
#                 timestamps                 = []
#                 lons                       = []
#                 lats                       = []
#                 for i in list(range(0,len(api_output)-1)):
#                     if len(api_output[i]['data'][0]) > 1:
#                         if 'grids' not in icollection:
#                             interpolated_profiles.append(interpolate_profiles(profile=api_output[i],levels_varname=varname_levels[icol_ind],levels_new=interp_levels))
#                         timestamps.append(dateutil.parser.isoparse(api_output[i]['timestamp'])) 
#                         lons.append(api_output[i]['geolocation']['coordinates'][0])
#                         lats.append(api_output[i]['geolocation']['coordinates'][1])
                
#                 # store what is needed for profile data
#                 if 'grids' not in icollection:
#                     # shape variable into something appropriate
#                     data = [x['data'] for x in interpolated_profiles]
#                     data = [[level[varname[icol_ind]] for level in x] for x in data]
#                     # store interpolated profiles
#                     regions_list_data_vert_interp.append(data)
#                     data = np.transpose(data)
#                     regions_list_data_horiz_ave.append(np.nanmean(data,1))
#                     regions_list_data_horiz_ave_levels.append(interp_levels)
#                     regions_list_data_raw_xarray.append([])
#                     regions_list_data_raw.append([x['data'] for x in api_output])
#                     ## quick plot of profiles in region
#                     map_lons_lats(lons,lats,dx=20,dy=20)

#                 ## store timestamp
#                 regions_list_data_time.append(timestamps)

#                 # store what is needed for grids
#                 if 'grids' in icollection:

#                     grids_opt  = {
#                                 "id": api_output[0]['metadata'][0]
#                                 }
#                     grids_meta = avh.query('grids/meta', options=grids_opt, verbose='true',apikey=API_KEY, apiroot=get_route(icollection))

#                     xar = grids_to_xarray(api_output,grids_meta)
#                     # store data
#                     regions_list_data_vert_interp.append([]) # no need to interpolate for grids
#                     regions_list_data_raw_xarray.append(xar)
#                     xar_h_ave = xarray_regional_mean(xar)['data'].mean(axis=1).values.flatten()
#                     regions_list_data_horiz_ave.append(xar_h_ave)
#                     regions_list_data_horiz_ave_levels.append(grids_meta[0]['levels'][0:len(xar_h_ave)])
#                     regions_list_data_raw.append([[x['data'][0],grids_meta[0]['levels'][0:len(x['data'][0])]] for x in api_output])

#                     # to see some info about the grid (including the units of the vertical level variable)
#                     grids_meta
                    
#     output = {}
#     for ivar in output_vars:
#         output[ivar] = globals()[ivar]
#     return output

# # function to plot the output of get_profiles_in_regions_and_horiz_ave: plot the horizontal average for (vertically) interpolated profiles and for the gridded product
# # input:
# # data_reg: output of the function get_preofiles_in_regions_and_horiz_ave
# # data_reg_cols     = ['k', 'b', 'r', 'm']
# # xlabel_tag        = 'Salinity, psu'
# def profiles_in_regions_and_horiz_ave_plot1d_horiz_ave(data_reg,data_reg_cols,xlabel_tag):
#     plt.figure(figsize=(10,8))
#     for i,idata in enumerate(data_reg['regions_list_data_horiz_ave']):
#         plt.plot(idata,data_reg['regions_list_data_horiz_ave_levels'][i],color=data_reg_cols[i],linewidth=3)
#     plt.gca().invert_yaxis()
#     plt.ylabel('Vertical level (m or dbar)',size=16)
#     plt.xlabel(xlabel_tag,size=16)
#     plt.legend([a_+', '+b_ for a_, b_ in zip(data_reg['regions_list_tags'],data_reg['regions_list_collections'])])

# # function to plot the output of get_profiles_in_regions_and_horiz_ave: let's look at all the raw profiles that were vertically interpolated (except for the gridded products) to then compute the horizontal average above
# # input:
# # data_reg: output of the function get_preofiles_in_regions_and_horiz_ave
# # data_reg_cols     = ['k', 'b', 'r', 'm']
# # xlabel_tag        = 'Salinity, psu'
# def profiles_in_regions_and_horiz_ave_plot1d_all(data_reg,data_reg_cols,xlabel_tag):
#     plt.figure(figsize=(15,8))
#     for i,idata in enumerate(data_reg['regions_list_data_raw']):
#         plt.subplot(1,len(data_reg['regions_list_data_horiz_ave']),i+1)    
#         for iidata in idata:
#             plt.plot(iidata[0],iidata[1],color=data_reg_cols[i])
#         plt.gca().invert_yaxis()
#         if i==0:
#             plt.ylabel('Vertical level (m or dbar)',size=16)
#         plt.title(data_reg['regions_list_tags'][i]+', '+data_reg['regions_list_collections'][i])
#         plt.xlabel(xlabel_tag,size=16)

# # function to plot the output of get_profiles_in_regions_and_horiz_ave: let's look at all the vertically interpolated profiles that were used to compute the horizontal average above
# # input:
# # data_reg: output of the function get_preofiles_in_regions_and_horiz_ave
# # data_reg_cols     = ['k', 'b', 'r', 'm']
# # xlabel_tag        = 'Salinity, psu'
# def profiles_in_regions_and_horiz_ave_plot1d_all_vert_interp(data_reg,data_reg_cols,xlabel_tag):
#     plt.figure(figsize=(15,8))
#     for i,idata in enumerate(data_reg['regions_list_data_vert_interp']):
#         if idata:
#             plt.subplot(1,len(data_reg['regions_list_data_horiz_ave']),i+1)    
#             for iidata in idata:
#                 plt.plot(iidata,data_reg['regions_list_data_horiz_ave_levels'][i],color=data_reg_cols[i])
#             plt.gca().invert_yaxis()
#             if i==0:
#                 plt.ylabel('Vertical level (m or dbar)',size=16)
#             plt.title(data_reg['regions_list_tags'][i]+', '+data_reg['regions_list_collections'][i]+' (interpolated)')
#             plt.xlabel(xlabel_tag,size=16)

# # function to plot the output of get_profiles_in_regions_and_horiz_ave: for each product, plot the raw profiles color coded by month/group of months of the year
# # input:
# # data_reg: output of the function get_preofiles_in_regions_and_horiz_ave
# # month_groups      = [[12, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11]]
# # month_groups_cols = ['dodgerblue', 'violet', 'orangered', 'gold']
# # month_groups_tags = ['DJF', 'MAM', 'JJA', 'SON']
# # xlabel_tag        = 'Salinity, psu'
# def profiles_in_regions_and_horiz_ave_plot1d_all_col_by_monthgroup(data_reg,month_groups,month_groups_cols,month_groups_tags,xlabel_tag):
#     fig = plt.figure(figsize=(15,8))
#     for i,idata in enumerate(data_reg['regions_list_data_raw']):
#         itime = data_reg['regions_list_data_time'][i]

#         plt.subplot(1,len(data_reg['regions_list_data_horiz_ave']),i+1)    
#         for ii,iidata in enumerate(idata):
#             for imm,mm in enumerate(month_groups):
#                 if ii<=len(itime)-1 and itime[ii].month in mm:
#                     plt.plot(iidata[0],iidata[1],color=month_groups_cols[imm])
#                     # glodap is a time mean
#                     if 'glodap' in data_reg['regions_list_collections'][i]:
#                         plt.plot(iidata[0],iidata[1],color='k')
#         plt.gca().invert_yaxis()
#         if i==0:
#             plt.ylabel('Vertical level (m or dbar)',size=16)
#         plt.title(data_reg['regions_list_tags'][i]+', '+data_reg['regions_list_collections'][i])
#         plt.xlabel(xlabel_tag,size=16)
#     for i,ileg in enumerate(month_groups_tags):
#         text = fig.text(0.5+((i-1)*(1/len(month_groups_cols)/3)), 0.02,ileg,horizontalalignment='center', wrap=True,color=month_groups_cols[i],size=16) 