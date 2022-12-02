import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import time, json, pprint, copy, math, dateutil
from datetime import datetime, timedelta
import scipy.interpolate
import numpy as np
import pandas

from argovisHelpers import helpers as avh

def polygon_lon_lat(polygon_str):
    # polygon_str: string value of polygon search parameter, ie "[[lon0,lat0],[lon1,lat1],...,[lon0,lat0]]"
    # convert the polygon shape to lon and lat and save in a dictionary
    polygon_lon_lat_dict = {'lon': [float(i) for i in ((polygon_str.replace('[','')).replace(']','')).split(',')[0::2]], \
                    'lat': [float(i) for i in ((polygon_str.replace('[','')).replace(']','')).split(',')[1::2]]
                   }
    return polygon_lon_lat_dict   
    
def padlist(ls, length, token=None):
    x = copy.copy(ls)
    x = [float(i) for i in x]
    if len(ls) < length:
        tail = [token]*(length - len(ls))
        x += tail
    return x

def varrange(profiles, var):
    # given a list of profiles and a variable name,
    # return the globally min, max limits of that variable
    
    data = [p['data'] for p in profiles]
    data = [j for sub in data for j in sub]
    data = [level[var] for level in data if var in level]
    return [min(data), max(data)]

def interpolate(profile, levels):
    # given a <profile> and a list of desired pressure <levels>,
    # return a profile with profile.data levels at the desired pressure levels, with all available data interpolated to match
    # drop all QC and note `data_interpolated` in profile.data_warnings
    
    data_names = ['pressure']
    interpolated_data = [levels]
    for key in profile['data_keys']:
        if '_argoqc' not in key and '_woceqc' not in key and key!='pressure':
            finites = [(level['pressure'], level[key]) for level in profile['data'] if level['pressure'] is not None and level[key] is not None and not math.isnan(level['pressure']) and not math.isnan(level[key])]
            pressure = [x[0] for x in finites]
            data = [x[1] for x in finites]
            data_names.append(key)
            levels_in_range = [x for x in levels if x<pressure[-1]]
            interpolated_data.append(padlist(scipy.interpolate.pchip_interpolate(pressure, data, levels_in_range), len(levels), 0 ) )
    
    interpolated_levels = list(zip(*interpolated_data))
    data = [{data_names[i]:d[i] for i in range(len(data_names))} for d in interpolated_levels]
    interpolated_profile = copy.deepcopy(profile) # don't mutate the original
    interpolated_profile['data'] = data
    if 'data_warnings' in interpolated_profile:
        interpolated_profile['data_warnings'].append('data_interpolated')
    else:
        interpolated_profile['data_warnings'] = ['data_interpolated']
    return interpolated_profile

### plotting
def plot_xycol(x,y,d2col,x_tag,y_tag,d2col_tag,fontsz = 20):
    fig, ax = plt.subplots()
    
    plt.scatter(x,y,c=d2col)
    ax.invert_yaxis()
    ax.tick_params(axis='both', which='major', labelsize=fontsz)
    
    plt.xlabel(x_tag,fontsize=fontsz)
    plt.ylabel(y_tag,fontsize=fontsz)
    cb = plt.colorbar()
    cb.set_label(d2col_tag,fontsize=fontsz)
    cb.ax.tick_params(labelsize=fontsz) 
    # plt.show()
    
def simple_map(longitudes, latitudes, z=None, polygon=None, title='', fig=None, figIndex=None, marker=None, secondaries=None, bounding_box=None):
    fontsz = 20
    if fig:
        # ax = fig.add_subplot(figIndex[0], figIndex[1], figIndex[2], projection=ccrs.LambertConformal(cutoff=-60))
        ax = fig.add_subplot(figIndex[0], figIndex[1], figIndex[2], projection=ccrs.PlateCarree())
    else:
        fig = plt.figure(figsize=(20,10))
        # ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal(cutoff=-60))
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    if bounding_box:
        ax.set_extent(bounding_box)
    gl = ax.gridlines(draw_labels=True,color='black')
    if z:
        s = ax.scatter(longitudes, latitudes, c=z, transform=ccrs.PlateCarree(), s=100)
        plt.colorbar(s, pad=0.1)
    else:
        s = ax.scatter(longitudes, latitudes,transform=ccrs.PlateCarree(), s=100)
    
    if polygon:
        plt.plot(polygon_lon_lat(polygon)['lon'],polygon_lon_lat(polygon)['lat'],'-k',linewidth=10,transform=ccrs.PlateCarree()) 
    if marker:
        plt.plot(marker[0],marker[1],'Xr', transform=ccrs.PlateCarree(), markersize=20)
    if secondaries:
        for sec in secondaries:
            ax.scatter(sec['lon'], sec['lat'],transform=ccrs.PlateCarree(), color='red', marker='x', s=100)
    
    plt.rcParams['font.size'] = fontsz
    
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.LAND)
    plt.title(title, fontdict={'fontsize':fontsz})
    
def argo_heatmap(profiles, variable, interp_levels):
    interpolated_profiles = [interpolate(p,interp_levels) for p in profiles]
    timestamps = [dateutil.parser.isoparse(p['timestamp']) for p in interpolated_profiles]
    
    # pcolormesh needs arrays of bin boundaries; use the x,y values as the lower bound like a histogram, 
    # and set the final upper bound to make a bin the same width as previous bin
    level_bounds = copy.copy(interp_levels)
    level_bounds.append(interp_levels[-1] + interp_levels[-1] - interp_levels[-2])
    time_bounds = copy.copy(timestamps)
    time_bounds.append(timestamps[-1] + (timestamps[-1] - timestamps[-2]) )
    
    # shape variable into something appropriate
    data = [x['data'] for x in interpolated_profiles]
    data = [[level[variable] for level in x] for x in data]
    data = np.transpose(data)
    
    fig, ax = plt.subplots(figsize=(15, 10))
    pcm = ax.pcolormesh(time_bounds, level_bounds, data, vmin=varrange(profiles,variable)[0], vmax=varrange(profiles,variable)[1])
    ax.invert_yaxis()
    ax.set_xlabel('Timestamp')
    ax.set_ylabel('Pressure [mbar]')
    fig.colorbar(pcm, label=variable + ' [' + profiles[0]['units'][variable]+']')
    type(ax)
    
def hurrplot(bracket_points, tc, argo, var,colorBefore,colorAfter, maxpress=150, line=False):
    markers = ['.', 'v', 'P', 'X']
    markernames = ['dot', 'triangle', 'plus', 'x']
    fnt_sz = 40
    if var == 'temperature':
        var_units = ', deg C'
    elif var == 'salinity':
        var_units = ', psu'
    else:
        var_units = ''
    #fig, axs = plt.subplots(len(bracket_points), figsize=(10,10*len(bracket_points)))
    #n = 0
    for i in bracket_points:
        colo = argo[i]
        hurrtime = avh.parsetime(tc[i]['timestamp'])
        hurrspeed = tc[i]['data'][0]['wind']
        hurrpress = tc[i]['data'][0]['surface_pressure']
        hurrlon = tc[i]['geolocation']['coordinates'][0]
        hurrlat = tc[i]['geolocation']['coordinates'][1]
        annotation = 'Hurricane track longitude: ' + str(hurrlon) + '\nHurricane track latitude: ' + str(hurrlat) + \
        '\nHurricane track timestamp: ' + str(hurrtime) + '\nHurricane wind speed [kt]: ' + str(hurrspeed) + \
        '\nHurricane surface pressure [mb]: ' + str(hurrpress)
        np = 0
        for p in colo:
            print(p['data_keys'])
            ptime = avh.parsetime(p['timestamp'])
            c = colorBefore
            if ptime > hurrtime:
                c = colorAfter
            time2hurricane = str(hurrtime.replace(microsecond=0) - ptime.replace(microsecond=0))
            pressure = [level[p['data_keys'].index('pressure')] for level in p['data']]
            d = [level[p['data_keys'].index(var)] for level in p['data']]
            cutoff = next((i for i,v in enumerate(pressure) if v>maxpress))-1
            psub = pressure[0:cutoff]
            dsub = d[0:cutoff]
            
            if np == 0:
                fig = plt.figure(figsize=(10,10))
                ax  = fig.add_subplot(111)
            
            ax.scatter(dsub, psub, c=[c]*len(psub), marker=markers[np%len(markers)])
            if line:
                ax.plot(dsub, psub, c=c,linewidth=2)
            #ax.set(xlabel=var+var_units, ylabel='Pressure, dbar')
            plt.xlabel(var+var_units, fontsize=fnt_sz)
            plt.ylabel('pressure, dbar', fontsize=fnt_sz)
            ax.tick_params(axis='both', which='major', labelsize=fnt_sz)
            annotation += '\n\nArgo profile ' + p['_id'] + '\nmarker: ' + markernames[np%len(markers)]  + \
            '\nProfile longitude: ' + str(p['geolocation']['coordinates'][0]) + '\nProfile latitude: ' + \
            str(p['geolocation']['coordinates'][1]) + '\nProfile timestamp: ' + str(p['timestamp'])
            np+=1
        ax.invert_yaxis()
        ax.text(ax.get_xlim()[1] + 0.05*(ax.get_xlim()[1]-ax.get_xlim()[0]),0.9*maxpress, annotation, fontsize=fnt_sz)
        #n+=1  

def compare_plots(profile_list_1, var1, profile_list_2, var2):
    fig = plt.figure()
    fig, ax1 = plt.subplots(figsize=(15, 10))
    markers = ['o', '^', 's', '*', 'x']
    for i, p in enumerate(profile_list_1):
        pressures = [level['pressure'] for level in p['data']]
        var = [level[var1] for level in p['data']]
        ax1.scatter(var, pressures, s=10, c='b', marker=markers[i%len(markers)], label=p['_id'] + ' / ' + var1)
    for i, p in enumerate(profile_list_2):
        pressures = [level['pressure'] for level in p['data']]
        var = [level[var2] for level in p['data']]
        ax1.scatter(var, pressures, s=10, c='r', marker=markers[i%len(markers)], label=p['_id'] + ' / ' + var2)
        
    ax1.invert_yaxis()
    ax1.set_xlabel(var1 + ' [' + profile_list_1[0]['units'][var1] + '], ' + var2 + ' [' + profile_list_2[0]['units'][var2] + ']')
    ax1.set_ylabel('Pressure [dbar]')
    plt.legend(loc='upper left', frameon=False)
