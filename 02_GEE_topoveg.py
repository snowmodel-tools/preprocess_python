#!/usr/bin/env python
# coding: utf-8

# In[1]:


import ee
import geemap
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import requests
from rasterio.warp import transform
import xarray as xr
# Initialize the Earth Engine module.
ee.Initialize()


# In[2]:


#########################################################################
############################ USER INPUTS ################################
#########################################################################
# PATHS
# path to temporary folder to store tif files from gee
TIFpath = 'GEE_Downloads/'

# DOMAIN
# choose the modeling domain
domain = 'MT'

# path to directory where you want your output .tif and .asc files
dataPath = '/nfs/attic/dfh/Aragon2/CSOdmn/'+domain+'/'
#########################################################################


# # Download DEM and LC data from GEE

# In[3]:


# Download CFSv2 met data function
def get_topoveg(domain, OUTpath):
    
    #path to CSO domains
    domains_resp = requests.get("https://raw.githubusercontent.com/snowmodel-tools/preprocess_python/master/CSO_domains.json")
    domains = domains_resp.json()

    '''
    // These are the min and max corners of your domain in Lat, Long
    // Western Wyoming
    // Input the minimum lat, lower left corner
    '''
    minLat = domains[domain]['Bbox']['latmin']
    #// Input the minimum long, lower left corner
    minLong = domains[domain]['Bbox']['lonmin']
    #// Input the max lat, upper right corner
    maxLat = domains[domain]['Bbox']['latmax']
    #// Input the max Long, upper right corner
    maxLong = domains[domain]['Bbox']['lonmax']


    # This resolution for the NLCD and DEM outputs for the SnowModel domain in meters
    sm_resolution = int(domains[domain]['cellsize'])

    '''// Resolution for the PRISM output. This shoud change by Latitude of the domain
    // because the PRISM product spatial resolution is 2.5 minutes, which equals 150 arc seconds.
    // You can use this arc-second calculator to estimate the correct value for the PRISM resolution by latitude
    // https://opendem.info/arc2meters.html
    // This is one arc-second in meters for 43 degrees N Latitude'''
    one_arcsecond = 22.57
    PRISM_resolution = one_arcsecond * 150

    '''// Define the final output projection using EPSG codes'''
    epsg_code = domains[domain]['mod_proj']

    my_domain = ee.Geometry.Rectangle(**{'coords':[minLong,minLat,maxLong,maxLat],'proj': 'EPSG:4326','geodesic':True,});
    
        #################### adjust for AK ##########################
    #need to specifiy different GEE products for above 60 deg lat
    #if maxLat>60:
    
    # download 30m SRTM data
    #// NOTE: This only covers through 60 degrees latitude. See visualization layers.
    SRTM30 = ee.Image('USGS/SRTMGL1_003')
    #image = SRTM30.clip(my_domain).unmask()
    filename = os.path.join(OUTpath, 'DEM_'+domain+'.tif')
    geemap.ee_export_image(SRTM30, filename=filename, scale=sm_resolution, region=my_domain, crs = epsg_code)
    
    # download NLCD data
    NLCD = ee.ImageCollection('USGS/NLCD');
    landcover = NLCD.select('landcover');
    #// Define the timeframe of NLCD images to select. Currently, its set to the previous 5 years.
    landcoverfiltered=landcover.filterDate('2015-01-01','2020-01-01');
    #// Create a single image out of the image collection using the most common land cover 
    #// designation from the previous 5 years.
    lcsingle=landcoverfiltered.mode();
    filename = os.path.join(OUTpath, 'NLCD2016_'+domain+'.tif')
    geemap.ee_export_image(lcsingle, filename=filename, scale=sm_resolution, region=my_domain, file_per_band=False, crs = epsg_code)


# In[4]:


# execute GEE function
get_topoveg(domain, dataPath)


# # DEM

# In[5]:


# dem
def DEM2SM(INfile, OUTpath):
    da = xr.open_rasterio(INfile)
    
    #ascii header 
    head = "ncols\t"+str(da.shape[2])+"\n"     "nrows\t"+str(da.shape[1])+"\n"     "xllcorner\t"+str(int(min(da.x.values)-da.res[0]/2))+"\n"     "yllcorner\t"+str(int(min(da.y.values)-da.res[0]/2))+"\n"     "cellsize\t"+str(int(da.res[0]))+"\n"     "NODATA_value\t-9999"    
    
    np.savetxt(OUTpath+'DEM_'+domain+'.asc', np.squeeze(da.values), fmt='%d', header = head,comments='')


# # Landcover 
# 
# ### NLCD LC Codes
# 
# |Code|NLCD2016 Landclass|Code|NLCD2016 Landclass|
# | :-: | :-: | :-: | :-: |
# |11 | open water |51 | dwarf shrub|
# |12| ice / snow |52 | shrub/scrub|
# |21 | developed; open space |71 | grassland/herbaceous|
# |22 | developed; low intensity|72 | hedge/herbaceous|
# |23 | developed; med intensity|73 | lichens|
# |24 | developed; high intensity|74 | moss|
# |31 | barren; rock, sand, clay|81 | pasture/hay|
# |41 | deciduous forest|82 | cultivated crops|
# |42 | evergreen forest|90 | woody wetlands|
# |43 | mixed shrub|95 | emergent herbaceous wetlands|
# 
# 
# ### Snowmodel LC Codes
# 
# |Code  |Landcover Class |Code  |Landcover Class |
# | --- | --- | --- | --- |
# |1     | coniferous forest |13    | subalpine meadow  |      
# |2     | deciduous forest |14    | tundra (non-tussock) |      
# |3     | mixed forest |15    | tundra (tussock) |           
# |4     | scattered short-conifer |16    | prostrate shrub tundra | 
# |5     | clearcut conifer |17    | arctic gram. wetland |       
# |6     | mesic upland shrub |18    | bare |       
# |7     | xeric upland shrub |19    | water/possibly frozen |       
# |8     | playa shrubland |20    | permanent snow/glacier |         
# |9     | shrub wetland/riparian |21    | residential/urban |   
# |10    | erect shrub tundra |22    | tall crops |       
# |11    | low shrub tundra |23    | short crops |        
# |12    | grassland rangeland  |24    | ocean |    

# In[6]:


# landcover data 
def LC2SM(INfile,OUTpath):
    da = xr.open_rasterio(INfile)
    data = np.squeeze(da.values)

    #ascii header 
    head = "ncols\t"+str(da.shape[2])+"\n"     "nrows\t"+str(da.shape[1])+"\n"     "xllcorner\t"+str(int(min(da.x.values)-da.res[0]/2))+"\n"     "yllcorner\t"+str(int(min(da.y.values)-da.res[0]/2))+"\n"     "cellsize\t"+str(int(da.res[0]))+"\n"     "NODATA_value\t-9999"
    
    #reassign lc from NLCD to SM classes
    DIR=np.empty([da.shape[1],da.shape[2]])
    DIR[data == 11 ]=19
    DIR[data == 12 ]=20
    DIR[data == 21 ]=21
    DIR[data == 22 ]=21
    DIR[data == 23 ]=21
    DIR[data == 24 ]=21
    DIR[data == 31 ]=18
    DIR[data == 41 ]=2
    DIR[data == 42 ]=1
    DIR[data == 43 ]=6
    DIR[data == 51 ]=6
    DIR[data == 52 ]=6
    DIR[data == 71 ]=12
    DIR[data == 72 ]=12
    DIR[data == 73 ]=12
    DIR[data == 74 ]=12
    DIR[data == 81 ]=23
    DIR[data == 82 ]=22
    DIR[data == 90 ]=9
    DIR[data == 95 ]=9
    DIR.astype(int)
    np.savetxt(OUTpath+'NLCD2016_'+domain+'.asc', DIR, fmt='%d', header = head,comments='')


# # Lat long grids

# In[7]:


def LTLN2SM(INfile,OUTpath):
    da = xr.open_rasterio(INfile)

    # Compute the lon/lat coordinates with rasterio.warp.transform
    ny, nx = len(da.y), len(da.x)
    x, y = np.meshgrid(da.x, da.y)

    # Rasterio works with 1D arrays
    lon, lat = transform(da.crs, {'init': 'EPSG:4326'},
                         x.flatten(), y.flatten())
    lon = np.asarray(lon).reshape((ny, nx))
    lat = np.asarray(lat).reshape((ny, nx))

    #ascii header 
    head = "ncols\t"+str(da.shape[2])+"\n"     "nrows\t"+str(da.shape[1])+"\n"     "xllcorner\t"+str(int(min(da.x.values)-da.res[0]/2))+"\n"     "yllcorner\t"+str(int(min(da.y.values)-da.res[0]/2))+"\n"     "cellsize\t"+str(int(da.res[0]))+"\n"     "NODATA_value\t-9999"
    np.savetxt(OUTpath+'grid_lat_'+domain+'.asc', lat, fmt='%2.5f', header = head,comments='')
    np.savetxt(OUTpath+'grid_lon_'+domain+'.asc', lon, fmt='%4.5f', header = head,comments='')


# # Execute functions

# In[8]:


# generate topo
INfile = dataPath+'DEM_'+domain+'.tif'
DEM2SM(INfile, dataPath)
#generate veg
INfile = dataPath+'NLCD2016_'+domain+'.tif'
LC2SM(INfile, dataPath)
#generate lat lon grids
LTLN2SM(INfile,dataPath)


# In[ ]:




