import requests, xarray, pandas, math, datetime, copy
import numpy as np
from datetime import datetime, timedelta
from argovisHelpers import helpers as avh
import matplotlib.pyplot as plt

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
    
    metadata_params = {
        "id": data_raw[0]['metadata']
        }

    metadata = avh.query('grids/meta', options=metadata_params, apikey=API_KEY, apiroot=API_PREFIX)
    
    return xargrid(grid=data_raw, depths=metadata[0]['levels'],long_conversion_type=long_conversion_type)

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
    plt.show()
    return
    

def plot_map(data,cf_levels,ylim_bottom,ylim_top,plot_title,font_size):
    # just to fix bug
    plt.plot()
    plt.rcParams['font.size'] = font_size
    plt.close()
    #
    data.plot.contourf(levels=cf_levels,aspect=4, size=6)
    plt.ylim(bottom=ylim_bottom,top=ylim_top)
    plt.rcParams['font.size'] = font_size
    plt.title(plot_title)
    plt.show()