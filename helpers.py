import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import random, numpy, pandas, matplotlib
from argovisHelpers import helpers as avh
import xarray as xr

def grids_to_xarray(grids,grids_meta):
    data_list        = []
    data_list_lev    = []
    data_list_lon    = []
    data_list_lat    = []
    data_list_tstamp = []
    for x in grids:
        for ix,x_at_lev in enumerate(x['data'][0]):
            data_list.append(x_at_lev)
            data_list_lev.append(grids_meta[0]['levels'][ix])
            data_list_lon.append(x['geolocation']['coordinates'][0])
            data_list_lat.append(x['geolocation']['coordinates'][1])
            data_list_tstamp.append(x['timestamp'])
    bfr_lon = numpy.array(data_list_lon)
    bfr_lon[bfr_lon<20] = bfr_lon[bfr_lon<20]+360

    data_list_lon = bfr_lon.tolist()

    data_dict = {'data': data_list,'latitude':data_list_lat, 'longitude': data_list_lon,'levels': data_list_lev} 

    data_df = pandas.DataFrame(data_dict)   
    df_rows = pandas.DataFrame(data_df).set_index(["latitude", "longitude","levels" ])
    return xr.Dataset.from_dataframe(df_rows)

def polygon_lon_lat(polygon_str):
    # polygon_str: string value of polygon search parameter, ie "[[lon0,lat0],[lon1,lat1],...,[lon0,lat0]]"
    # convert the polygon shape to lon and lat and save in a dictionary
    polygon_lon_lat_dict = {'lon': [float(i) for i in ((polygon_str.replace('[','')).replace(']','')).split(',')[0::2]], \
                    'lat': [float(i) for i in ((polygon_str.replace('[','')).replace(']','')).split(',')[1::2]]
                   }
    return polygon_lon_lat_dict

def simple_map(longitudes, latitudes, z=None, zlabel=None, markers=None, polygon=None, polygonweight=5, title='', fig=None, figIndex=None, spot=None, suppress_colorbar=False, font_size=20, secondaries=None):
    if fig:
        ax = fig.add_subplot(figIndex[0], figIndex[1], figIndex[2], projection=ccrs.LambertConformal())
    else:
        fig = plt.figure(figsize=(20,10))
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal())
    gl = ax.gridlines(draw_labels=True,color='black')
    if z:
        s = ax.scatter(longitudes, latitudes, c=z, transform=ccrs.PlateCarree())
        if not suppress_colorbar:
            cb = plt.colorbar(s, pad=0.1)
            if zlabel:
                cb.set_label(zlabel, rotation=270, labelpad=20, size=font_size)
    elif markers:
        lons = numpy.array(longitudes)
        lats = numpy.array(latitudes)
        marks = numpy.array(markers)
        ms = set(markers)
        for m in ms:
            s = ax.scatter(lons[marks==m],lats[marks==m], marker=m,transform=ccrs.PlateCarree())
    else:
        s = ax.scatter(longitudes, latitudes,transform=ccrs.PlateCarree())
    if polygon:
        plt.plot(polygon_lon_lat(str(polygon))['lon'],polygon_lon_lat(str(polygon))['lat'],'-k',transform=ccrs.PlateCarree(), linewidth=polygonweight, color='red')
    if secondaries:
        for sec in secondaries:
            ax.scatter(sec['lon'], sec['lat'],transform=ccrs.PlateCarree(), color=sec['color'], marker='x', s=150, linewidth=3)
    if spot:
        plt.plot(spot[0],spot[1],'Xr', transform=ccrs.PlateCarree(), markersize=20)
    plt.rcParams['font.size'] = font_size
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.LAND)
    plt.title(title, fontdict={'fontsize':font_size})
    
def random_color():
    return "#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])

def marker_map(items):
    # pass in a list of unique <items>, get back a dictionary mapping them onto different markers.
    options = ['.', 'o', 'v', '^', '<', '>', '1', '2', '3', '4', '8', 's', 'p', 'P', '*', 'h', 'H', '+', 'x', 'X', 'D', 'd', '$a$', '$b$', '$c$', '$d$', '$e$', '$f$', '$g$', '$h$']
    return {item: options[i%len(options)] for i, item in enumerate(items)}

def map_count_in_bins(df, startDate, endDate, outname='mapbins',
                      x_bound=[-180,181],y_bound=[-90,91],dx=1,dy=1,fpath='./',vmin_map=1,vmax_map=1000): 
    
    x_edges = numpy.arange(x_bound[0],x_bound[1],dx)
    y_edges = numpy.arange(y_bound[0],y_bound[1],dy)
    
    lons = df['longitudes']
    lats  = df['latitudes']
    datetimes  = df['timestamps']
    
    fig = plt.figure(figsize=(21, 7))
    
    h = numpy.histogram2d(lons, lats,[x_edges, y_edges])
    ax = plt.axes(projection=ccrs.PlateCarree())

    plt.pcolormesh(x_edges[0:-1]+dx/2, y_edges[0:-1]+dy/2, h[0].transpose(),
             transform=ccrs.PlateCarree(),norm=matplotlib.colors.LogNorm(),vmin=vmin_map,vmax=vmax_map)

    ax.stock_img()
    ax.coastlines()
    
    plt.colorbar()

    plt.savefig(fpath+outname+'_'+startDate[0:10]+'_'+endDate[0:10]+'.png')
    plt.show()
    
def mapping_df(api_returns):
    if 'geolocation' in api_returns[0]:
        longitudes = [x['geolocation']['coordinates'][0] for x in api_returns]
        latitudes  = [x['geolocation']['coordinates'][1] for x in api_returns]
        timestamps = [avh.parsetime(x['timestamp']) for x in api_returns]
    else:
        longitudes = [x[1] for x in api_returns]
        latitudes  = [x[2] for x in api_returns]
        timestamps = [avh.parsetime(x[3]) for x in api_returns]
        
    return pandas.DataFrame(zip(longitudes,latitudes,timestamps), columns=['longitudes', 'latitudes', 'timestamps'])

def level_df(api_returns, measurements, per_level_pressures=None, timesteps=None, index=None):
    # api_returns: list of data documents
    # measurements: list of valid measurement names from data_info[0] if data included; can also include 'latitude', 'longitude' 'timestamp', 'id' and 'months'
    # per_level_pressure: list of pressures, per what's found on a gridded product's metadata
    # timesteps: list of timesteps per what's found on a timeseries metadata
    # index: subset list of measurements to turn into dataframe index
    # returns a dataframe with one column for everything in measurements, and one row for every level in every document in api_returns
    
    ## make flat lists of all requested variables, corresponding by index
    flat_meas = {}
    
    if timesteps:
        ts = [avh.parsetime(t) for t in timesteps]
    
    for m in measurements:
        if m == 'months':
            if timesteps:
                mths = [x.month for x in ts]
                flat_meas['months'] = [mths*len(api_returns)]
                flat_meas['months'] = [item for sublist in flat_meas['months'] for item in sublist]
            else:
                flat_meas['months'] = [avh.parsetime(x['timestamp']).month for x in api_returns]
                flat_meas['months'] = [[flat_meas['months'][i]]*len(api_returns[i]['data'][0]) for i, x in enumerate(flat_meas['months'])]
            flat_meas['months'] = [item for sublist in flat_meas['months'] for item in sublist]
        elif m == 'longitude':
            flat_meas['longitude'] = [x['geolocation']['coordinates'][0] for x in api_returns]
            flat_meas['longitude'] = [[flat_meas['longitude'][i]]*len(api_returns[i]['data'][0]) for i, x in enumerate(flat_meas['longitude'])]
            flat_meas['longitude'] = [item for sublist in flat_meas['longitude'] for item in sublist]
        elif m == 'latitude':
            flat_meas['latitude'] = [x['geolocation']['coordinates'][1] for x in api_returns]
            flat_meas['latitude'] = [[flat_meas['latitude'][i]]*len(api_returns[i]['data'][0]) for i, x in enumerate(flat_meas['latitude'])]
            flat_meas['latitude'] = [item for sublist in flat_meas['latitude'] for item in sublist]
        elif m == 'timestamp' and 'timestamp' in api_returns[0]:
            flat_meas['timestamp'] = [avh.parsetime(x['timestamp']) for x in api_returns]
            flat_meas['timestamp'] = [[flat_meas['timestamp'][i]]*len(api_returns[i]['data'][0]) for i, x in enumerate(flat_meas['timestamp'])]
            flat_meas['timestamp'] = [item for sublist in flat_meas['timestamp'] for item in sublist]
        elif m == 'id':
            flat_meas['id'] = [x['_id'] for x in api_returns]
            flat_meas['id'] = [[flat_meas['id'][i]]*len(api_returns[i]['data'][0]) for i, x in enumerate(flat_meas['id'])]
            flat_meas['id'] = [item for sublist in flat_meas['id'] for item in sublist]
        else:
            flat_meas[m] = [x['data'][x['data_info'][0].index(m)] for x in api_returns]
            flat_meas[m] = [item for sublist in flat_meas[m] for item in sublist]
    if per_level_pressures:
        flat_meas['pressure'] = [per_level_pressures*len(api_returns)]
        flat_meas['pressure'] = [item for sublist in flat_meas['pressure'] for item in sublist]
    if timesteps:
        flat_meas['timestamp'] = [ts*len(api_returns)]
        flat_meas['timestamp'] = [item for sublist in flat_meas['timestamp'] for item in sublist]
        
    cols = flat_meas.keys()
    flat_meas = [flat_meas[key] for key in flat_meas.keys()]
    df = pandas.DataFrame(zip(*flat_meas), columns=cols)
    if index:
        df = df.set_index(index)
    return df

def regional_mean(dxr, form='area'):
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
    