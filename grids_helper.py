import requests, xarray, pandas, math, datetime, copy
import numpy as np
from datetime import datetime, timedelta
from argovisHelpers import helpers as avh
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# set up parameteres
def def_activity_param(str_activity,region_selected,section_selected_lat,depth_level_selected,API_KEY,API_PREFIX):
    if str_activity == 'activity_enso':
        dic = {'grid_name': 'temperature_rg', \
              'plot_title': 'Temperature, degC', \
              'region': def_reg_activity(str_region=region_selected), \
              'latitude_band_selected': def_zonal_section_activity(str_latitude=section_selected_lat), \
              'year_neutral_conditions': [2013], \
              'year_enso': [2007, 2009, 2010, 2015], \
              'month': 12, \
              'day_grid': 15, \
              'levels': def_levels_activity(str_levels=depth_level_selected), \
              'levels_section': '0,250', \
              'font_size_section': '26',\
              'font_size_map': '32',\
              'long_conversion_type': 'long20_380',\
              'cf_levels_maps':  np.arange(22,33,1), \
              'cf_levels_sections': np.arange(6,34,2), \
              'cf_levels_sections_line': np.arange(20,22,4), \
              'apikey': API_KEY, \
              'apiroot': API_PREFIX, \
              }
        # check the region name for the map
        if region_selected=='equatorial_atlantic' or region_selected=='equatorial_indian':
            print('HINT: this region will not help you describe the ENSO phases, try another one!')
        elif region_selected=='equatorial_pacific':
            print('Good job selecting the region for the map!')
        else:
            print('Before displaying the map, please select a valid name for the region')
            
        # check the level for the map
        if depth_level_selected=='near_surface':
            print('Good job selecting the depth level for the map!')
        elif depth_level_selected=='near_1500m':
            print('HINT: this depth level will not help you describe the ENSO phases, try another one!')
        else:
            print('Before displaying the map, please select a valid name for the depth level')
            
        # check the latitude for the section
        if section_selected_lat=='equatorial':
            print('Good job selecting the latitude for the section!')
        elif section_selected_lat=='subpolar':
            print('HINT: this latitude will not help you describe the ENSO phases, try another one!')
        elif section_selected_lat=='midlatitude':
            print('HINT: this latitude will not help you describe the ENSO phases, try another one!')
        else:
            print('Before displaying the section, please select a valid name for the latitude')
            
        return dic
    else:
        print('Please select a valid name for the activity')
        return

def def_reg_activity(str_region):
    # returns a list, e.g. [west long, east long, south_lat, north_lat]
    if str_region == 'equatorial_pacific':
        return [150.5,-119.5,-7,7]
    elif str_region == 'equatorial_indian':
        return [45,96,-7,7]
    elif str_region == 'equatorial_atlantic':
        return [-45.5,6.5,-7,7]
    else:
        return []
    
def def_levels_activity(str_levels):
    # returns a string, i.e. [0, 5]
    if str_levels == 'near_surface':
        return '0,5' 
    elif str_levels == 'near_1500m':
        return '1450,1550'
    else:
        return []
    
def def_zonal_section_activity(str_latitude):
    # returns a list, e.g. [south_lat, north_lat]
    if str_latitude == 'equatorial':
        return [-0.5,0.5]
    elif str_latitude == 'subpolar':
        return [54.5,56.5]
    elif str_latitude == 'midlatitude':
        return [44.5,46.5]
    else:
        return []

# query
def create_boxstr_for_query(longitude_west,longitude_east,latitude_south,latitude_north):
    boxstr = '[[' + str(longitude_west) + ',' + str(latitude_south) + \
        '], ['  + str(longitude_west) + ',' + str(latitude_north) + \
        '], [' + str(longitude_east) + ',' + str(latitude_north) + \
        '], [' + str(longitude_east) + ',' + str(latitude_south) + \
        '], [' + str(longitude_west) + ',' + str(latitude_south) + ']]'
    return boxstr

def query_grid_by_region_month_year(grid_name,region_str,levels, \
                                    long_conversion_type,\
                                    month_start,year_start,\
                                    month_end,year_end,\
                                    API_KEY,API_PREFIX):
    params = {
      "startDate": f'{year_start}'+'-'+f'{month_start:02}'+'-01T00:00:00Z',
      "endDate": f'{year_end}'+'-'+f'{month_end:02}'+'-31T00:00:00Z',
      "polygon": region_str,
      "data": grid_name,
      "presRange": levels,
      "compression": 'array',
    }

    data_raw = avh.query('grids/'+grid_name, options=params, apikey=API_KEY, apiroot=API_PREFIX)
    
    #     #### this below is not needed since we have presRange specified now (hence the pressure levels are returned in data_raw)

    #     metadata_params = {
    #             "id": data_raw[0]['metadata']
    #             }

    #     metadata = avh.query('grids/meta', options=metadata_params, apikey=API_KEY, apiroot=API_PREFIX)

    #     data = {**metadata,**data_raw}

    
    return xargrid(grid=data_raw, depths=data_raw[0]['levels'],long_conversion_type=long_conversion_type)

# process
def xargrid(grid, depths,long_conversion_type):
    # given the json response <grid> of a request to /grids/{gridName},
    # a list <depths> of the corresponding depths for these grid documents
    # return an xarray object with coordinates time, lat, lon, depth, and measurement value.
    
    lat = []
    lon = []
    time = []
    meas = []
    pressure = []
    for p in grid:
        for i, e in enumerate(p['data']):
            bfr_lon = p['geolocation']['coordinates'][0];
            if long_conversion_type=='long20_380':
                if bfr_lon < 20:
                    bfr_lon = bfr_lon + 360
            lon.append(bfr_lon)
            lat.append(p['geolocation']['coordinates'][1])
            time.append(avh.parsetime(p['timestamp']))
            meas.append(p['data'][i][0])
            pressure.append(depths[i])
            
    df = pandas.DataFrame({"latitude": lat, 
                           "longitude": lon, 
                           "time": time, 
                           "pressure": pressure, 
                           "data": meas}).set_index(["latitude","longitude","time","pressure"])
    return df.to_xarray()


def areaweighted_region_mean(dxr):
    # given an xarray dataset <grid> for a given depth and time,
    # calculate the mean of the gridded data variable, weighted by grid cell area
    weights = np.cos(np.deg2rad(dxr.latitude))
    weights.name = "weights"
    dxr_weighted = dxr.weighted(weights)
    
    return dxr_weighted.mean(("longitude", "latitude"))
    
# visualization
def run_activity_maps(activity,str_year):
    for year in activity[str_year]:
        data = query_grid_by_region_month_year(grid_name=activity['grid_name'],\
                                               region_str = create_boxstr_for_query(longitude_west=activity['region'][0],\
                                                                                 longitude_east=activity['region'][1], \
                                                                                 latitude_south=activity['region'][2],\
                                                                                 latitude_north=activity['region'][3]), \
                                            long_conversion_type=activity['long_conversion_type'],\
                                            levels=activity['levels'],\
                                            month_start=activity['month'],year_start=year,\
                                            month_end=activity['month'],year_end=year,\
                                            API_KEY=activity['apikey'],API_PREFIX=activity['apiroot'])

        plot_map_enso_activity(data=data["data"].mean(dim="time").mean(dim="pressure"),cf_levels=activity['cf_levels_maps'],\
                ylim_bottom=activity['region'][2]+1,ylim_top=activity['region'][3]-1,\
                 plot_title=str(year) + ' ' + activity['plot_title'],font_size=activity['font_size_map'])
        
def run_activity_sections(activity,str_year):
    for year in activity[str_year]:
        data = query_grid_by_region_month_year(grid_name=activity['grid_name'],\
                                               region_str = create_boxstr_for_query(longitude_west=activity['region'][0],\
                                                                                 longitude_east=activity['region'][1], \
                                                                                 latitude_south=activity['region'][2],\
                                                                                 latitude_north=activity['region'][3]), \
                                            long_conversion_type=activity['long_conversion_type'],\
                                            levels=activity['levels_section'],\
                                            month_start=activity['month'],year_start=year,\
                                            month_end=activity['month'],year_end=year,\
                                            API_KEY=activity['apikey'],API_PREFIX=activity['apiroot'])
        plot_section(data=data["data"].mean(dim="latitude").mean(dim="time"),xaxis=data["longitude"],yaxis=data["pressure"],\
                     cf_levels=activity['cf_levels_sections'],\
                     cf_levels_line=activity['cf_levels_sections_line'],\
                     ylim_bottom=200,plot_title=str(year) + ' ' + activity['plot_title'],font_size=activity['font_size_section'])
# print(data)
        
def plot_section(data,xaxis,yaxis,cf_levels,cf_levels_line,ylim_bottom,plot_title,font_size):
    # just to fix bug
    plt.plot()
    plt.rcParams['font.size'] = font_size
    plt.close()
    #
    data.plot.contourf(y="pressure",yincrease=False,aspect=2, size=6,levels=cf_levels)
    # data["data"].mean(dim="latitude").mean(dim="time")
    if len(cf_levels_line) != 0:
#         plt.contour(data["longitude"],data["pressure"],data["data"].mean(dim="latitude").mean(dim="time").transpose(),\
#                     levels=cf_levels_line)
        plt.contour(xaxis,yaxis,data.transpose(),levels=cf_levels_line)
    plt.ylim(bottom=ylim_bottom)
    plt.rcParams['font.size'] = font_size
    plt.title(plot_title)
    plt.xlabel('Longitude, Degrees East')
    plt.ylabel('Pressure, dbar')
    plt.show()
    return
    

def plot_map_enso_activity(data,cf_levels,ylim_bottom,ylim_top,plot_title,font_size):
    # just to fix bug
    plt.plot()
    plt.rcParams['font.size'] = font_size
    plt.close()
    #
    data.plot.contourf(levels=cf_levels,aspect=4, size=6)
    plt.ylim(bottom=ylim_bottom,top=ylim_top)
    plt.rcParams['font.size'] = font_size
    plt.title(plot_title)
    plt.xlabel('Longitude, Degrees East')
    plt.ylabel('Latitude, Degrees North')
    
    plt.gca().add_patch(patches.Rectangle((190, -5), 50, 10,\
                       fill = False,
                       color = "red",
                       linewidth = 10))

    plt.show()
    
def plot_line_pos_neg(data,data_time,data_time_delta_num,data_ylabel,data_title,font_size):
    plt.plot()
    plt.rcParams['font.size'] = font_size
    plt.close()
    
    plt.figure(figsize=(15, 5))
    for n in np.arange(-2,3,1):
        plt.axhline(y = n,color='k')
    for n in np.arange(0,len(data_time),data_time_delta_num):
        plt.axvline(x = data_time[n],color='k')
    #data_reg_ave_anom["data"].plot(size=10,aspect=3)
    plt.plot(data_time,data,linewidth=5,color='k')
    plt.xlim([data_time[0],data_time[-1]])

    plt.gca().fill_between(data_time,data,\
                           where=(data>0),color='red')
    plt.gca().fill_between(data_time,data,\
                           where=(data<0),color='blue')
    plt.ylabel(data_ylabel)
    plt.rcParams['font.size'] = font_size
    if data_title != '':
        plt.title(data_title)