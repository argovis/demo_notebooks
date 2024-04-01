import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import random, numpy, pandas, matplotlib
from argovisHelpers import helpers as avh
import xarray as xr

def get_route(collection_name):
    
    if collection_name == 'drifters':
        return 'https://argovisbeta01.colorado.edu/dapi/'
    else:
        return 'https://argovis-api.colorado.edu/'
    
