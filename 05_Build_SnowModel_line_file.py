#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import geopandas as gpd
import codecs
import requests
import rasterio as rio 


# In[2]:


#########################################################################
############################ USER INPUTS ################################
#########################################################################

# DOMAIN
# choose the modeling domain
domain = 'CO_S'

# PATHS
# path to domain data folder 
dataPath = '/nfs/attic/dfh/Aragon2/CSOdmn/'+domain+'/'
#########################################################################


# In[3]:


def build_snowmodel_line(domain,dataPath):

    # path to DEM ascii
    DEMpath = dataPath+'DEM_'+domain+'.asc'

    #path to SNOTEL gdf
    gdfpath = dataPath+'CSO_SNOTEL_sites_'+domain+'.geojson'

    #path to VEG .asc
    VEGpath = dataPath+'NLCD2016_'+domain+'.asc'

    #path to lat .asc
    LATpath = dataPath+'grid_lat_'+domain+'.asc'

    #path to lon .asc
    LONpath = dataPath+'grid_lon_'+domain+'.asc'

    # VEG outfile path
    outVEGpath = dataPath+'NLCD2016_'+domain+'_line.asc'

    # DEM outfile path
    outDEMpath = dataPath+'DEM_'+domain+'_line.asc'

    # Line outfile path
    outFpath = dataPath+'snowmodel_line_pts.dat'

    #lon outfile path
    outLONpath = dataPath+'grid_lon_'+domain+'_line.asc'

    #lat outfile path
    outLATpath = dataPath+'grid_lat_'+domain+'_line.asc'

    # station data
    stn_gdf = gpd.read_file(gdfpath)

    #path to CSO domain
    domains_resp = requests.get("https://raw.githubusercontent.com/snowmodel-tools/preprocess_python/master/CSO_domains.json")
    domains = domains_resp.json()
    # CSO projection
    mod_proj = domains[domain]['mod_proj']

    #get metadata from .asc 
    def getmeta(file_name):
        myvars = {}
        lines = open(file_name, 'r').readlines()
        for i in range(5):
            line=lines[i]
            name, var = line.partition("\t")[::2]
            myvars[name.strip()] = int(var)
        return myvars

    myvars = getmeta(DEMpath)

    #Build northing easting array 
    #Northing
    # locate center of lower left cell
    st = myvars['yllcorner']+(myvars['cellsize']/2)
    #build northing array
    north = np.arange(st,st+myvars['nrows']*myvars['cellsize'],myvars['cellsize'])

    #Easting
    # locate center of lower left cell
    st = myvars['xllcorner']+(myvars['cellsize']/2)
    #build easting array
    east = np.arange(st,st+myvars['ncols']*myvars['cellsize'],myvars['cellsize'])

    #fortran indexing starts at 1
    #staion 1 to N of N stations
    count = np.zeros(stn_gdf.shape[0])
    #easting of pixel corresponding to station 
    stn_est = np.zeros(stn_gdf.shape[0])
    #northing of pixel corresponding to station
    stn_nor = np.zeros(stn_gdf.shape[0])
    #index of pixel corresponding to station 
    est_idx = np.zeros(stn_gdf.shape[0])
    #index of pixel corresponding to station 
    nor_idx = np.zeros(stn_gdf.shape[0])

    for z in range(stn_gdf.shape[0]):
        count[z] = z + 1
        lons = abs(stn_gdf.easting[z]-east)
        loIDX = [i for i, value in enumerate(lons) if value == np.min(abs(stn_gdf.easting[z]-east))] 
        stn_est[z] = east[loIDX[0]]
        est_idx[z] = loIDX[0] + 1
        lats = abs(stn_gdf.northing[z]-north)
        laIDX = [i for i, value in enumerate(lats) if value == np.min(abs(stn_gdf.northing[z]-north))]
        stn_nor[z] = north[laIDX[0]]
        nor_idx[z] = laIDX[0] + 1

    #Print out .dat file 
    f= open(outFpath,"w+")
    for z in range(count.shape[0]):
        f.write('{:08.0f}\t'.format(count[z])+'{:10.0f}\t'.format(est_idx[z])+'{:10.0f}\t'.format(nor_idx[z])+
                '{:10.0f}\t'.format(stn_est[z])+'{:10.0f}\t\n'.format(stn_nor[z]))
    f.close() 

    ## Extract topo, veg, lat, and lon files for snowmodel_line
    new = stn_gdf.to_crs(mod_proj)

    with rio.open(DEMpath) as src:
        rows, cols = rio.transform.rowcol(src.transform, new.geometry.centroid.x, new.geometry.centroid.y)

    # function to extract .asc files at SNOTEL locations
    def field2line(filepath,outpath):
        with codecs.open(filepath, encoding='utf-8-sig') as f:
            data = np.loadtxt(f,skiprows=6)
        data_line=[]
        for i in range(stn_gdf.shape[0]):
            info = str(int(data[rows[i],cols[i]]))
            data_line.append(info)
        lines = open(filepath, 'r').readlines()
        head = 'ncols\t1\nnrows\t'+str(stn_gdf.shape[0])+'\n'+lines[2]+lines[3]+lines[4]+lines[5]
        data = ''
        for s in data_line:
            data += s +'\n'
        f= open(outpath,"w+")
        f.write(head+data)
        f.close()    

    #extract DEM
    field2line(DEMpath,outDEMpath)
    #extract VEG
    field2line(VEGpath,outVEGpath)
    #extract LAT
    field2line(LATpath,outLATpath)
    #extract LON
    field2line(LONpath,outLONpath)


# In[4]:


build_snowmodel_line(domain,dataPath)


# In[ ]:




