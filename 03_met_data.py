#!/usr/bin/env python
# coding: utf-8


import ee
import geemap
import numpy as np
import matplotlib.pyplot as plt
import os
#from paths import *
import requests
import pandas as pd
import xarray as xr
from os import listdir
from datetime import datetime, timedelta, date
import contextlib

# Initialize the Earth Engine module
ee.Initialize()


#########################################################################
############################ USER INPUTS ################################
#########################################################################
# DOMAIN
# choose the modeling domain
domain = 'CHENA'
print(domain)
#path to CSO domains
domains_resp = requests.get("https://raw.githubusercontent.com/snowmodel-tools/preprocess_python/master/CSO_domains.json")
domains = domains_resp.json()

# PATHS
# path to temporary folder to store tif files from gee
TIFpath = 'GEE_Downloads_tmp/'
# path to where you want your output met .dat fime
OUTpath = '/nfs/depot/cce_u1/hill/dfh/op_snowmodel/chena_snowmodel/met/mm_'+domain+'_2016-2018.dat'

# TIME
# choose if want to set 'manual' or 'auto' date 
date_flag = 'manual'
# If you choose 'manual' set your dates below  
# This will start on the 'begin' date at 0:00 and the last iteration will 
# be on the day before the 'end' date below.
st_dt = '2016-10-01'#domains[domain]['st']
ed_dt = '2018-10-01'#domains[domain]['ed']
#########################################################################


# Date setup function
def set_dates(st_dt,ed_dt,date_flag):
    if date_flag == 'auto':
        # ###automatically select date based on today's date 
        hoy = date.today()
        antes = timedelta(days = 2)
        #end date 3 days before today's date
        fecha = hoy - antes
        eddt = fecha.strftime("%Y-%m-%d") 
        #start date
        if fecha.month <10:
            styr = fecha.year - 1
        else:
            styr = fecha.year
        stdt = str(styr)+'-10-01'
    elif date_flag == 'manual':
        stdt = st_dt
        eddt = (datetime.strptime(ed_dt,'%Y-%m-%d') + timedelta(days = 1)).strftime('%Y-%m-%d')
    return stdt, eddt


# Download CFSv2 met data function
def get_cfsv2(domain, TIFpath, stdt, eddt):
    # in GEE the last iteration is on the day before the 'end' date below
    # we adjust this here since it is not intuative
    #eddt = (datetime.strptime(eddt, '%Y-%m-%d')+timedelta(days = 1)).strftime('%Y-%m-%d')
    
    #create directory with initiation date for ensemble if it doesn't exist
    get_ipython().system('mkdir -p $TIFpath')

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

    #/ These are the min and max corners of your reanalysis in Lat, Long (create a slightly larger box)
    #// Input the minimum lat, lower left corner
    minLatMET = (minLat - 0.25);
    #// print(minLat2);
    #// Input the minimum long, lower left corner
    minLongMET = (minLong - 0.5);
    #// Input the max lat, upper right corner
    maxLatMET = (maxLat + 0.25);
    #// Input the max Long, upper right corner
    maxLongMET = (maxLong + 0.5);

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

    #// Name the DEM output
    dem_name = 'DEM'
    #// Name the Land Cover output
    lc_name = 'NLCD2016'

    my_domain = ee.Geometry.Rectangle(**{'coords':[minLong,minLat,maxLong,maxLat],'proj': 'EPSG:4326','geodesic':True,});
    my_domain_met = ee.Geometry.Rectangle([minLongMET,minLatMET,maxLongMET,maxLatMET])
    
    # download reanalysis data
    cfsv2 = ee.ImageCollection('NOAA/CFSV2/FOR6H')         .filterBounds(my_domain_met)         .filter(ee.Filter.date(stdt,eddt))

    data = cfsv2.select('Temperature_height_above_ground',         'Geopotential_height_surface',         'u-component_of_wind_height_above_ground',         'v-component_of_wind_height_above_ground',         'Pressure_surface',         'Specific_humidity_height_above_ground',         'Precipitation_rate_surface_6_Hour_Average')
    with contextlib.redirect_stdout(None):
        geemap.ee_export_image_collection(data, out_dir=TIFpath,region=my_domain_met,scale=22200,crs=epsg_code)


# function to check for missing dates
def missing_slice_check(stdt, eddt, TIFpath):
    # create a 6-hourly timeseries with no missing values from the start to end date
    timesin = pd.date_range(start=stdt, end=eddt, freq='6H')[:-1]
    for time in timesin:
        nam = time.strftime('%Y%m%d%H')
        
    # compile list of all tif files downloaded from gee
    gee_times =[]

    for file in listdir(TIFpath):
        if file.endswith("tif"):
            datetmp = datetime.strptime(file[:-4], '%Y%m%d%H')
            gee_times.append(datetmp)
    gee_times = sorted(gee_times)

    # check for to see if all time slices downloaded from GEE
    if len(timesin) != len(gee_times):
        print('gee is missing timeslice(s):\n',timesin[~timesin.isin(gee_times)].values)
        
        # if 4 or more consecutive timeslices are missing - quit the function
        duration = []
        for i in range(len(gee_times)-1):
            time_delta = gee_times[i+1] - gee_times[i]
            duration.append(time_delta.total_seconds()/60/60)
        if max(duration) >= 24:    
            print('at least one full day of met data is missing - quitting function')

        # if there are less than 4 missing consecutive time slices 
        # repeat the met data from the last valid time slice 
        else:
            missing_idx = np.where(~timesin.isin(gee_times))
            missing_dt = timesin[missing_idx]
            for j in range(len(missing_dt)):
                if np.squeeze(missing_idx)[j] == 0:
                    print('choose earlier start date so missing time slices can be filled in')
                else:
                    pre_dt=TIFpath+timesin[np.squeeze(missing_idx)[j]-1].strftime('%Y%m%d%H')+'.tif'
                    mis_dt = TIFpath+timesin[np.squeeze(missing_idx)[j]].strftime('%Y%m%d%H')+'.tif' 
                    get_ipython().system('cp $pre_dt $mis_dt')
                    print('replaced', timesin[np.squeeze(missing_idx)[j]].strftime('%Y%m%d%H'),' with ', timesin[np.squeeze(missing_idx)[j]-1].strftime('%Y%m%d%H'))



# Format gee files for SnowModel function
def MET2SM(TIFpath, OUTpath, stdt, eddt):
    # create a 6-hourly timeseries with no missing values from the start to end date
    timesin = pd.date_range(start=stdt, end=eddt, freq='6H')[:-1]
    
    #load first tif to get dimensions
    ar = xr.open_rasterio(TIFpath+timesin[0].strftime('%Y%m%d%H')+'.tif')
    
    # empty arrays for each met variable
    T = np.empty((len(timesin),ar.shape[1],ar.shape[2]))
    Z = np.empty((len(timesin),ar.shape[1],ar.shape[2]))
    U = np.empty((len(timesin),ar.shape[1],ar.shape[2]))
    V = np.empty((len(timesin),ar.shape[1],ar.shape[2]))
    P = np.empty((len(timesin),ar.shape[1],ar.shape[2]))
    H = np.empty((len(timesin),ar.shape[1],ar.shape[2]))
    PR = np.empty((len(timesin),ar.shape[1],ar.shape[2]))

    # extract met data from tifs 
    for i in range(len(timesin)):

        #load tif file
        nam = TIFpath+timesin[i].strftime('%Y%m%d%H')+'.tif'
        ar = xr.open_rasterio(nam)
        T[i,:,:] = ar[0,:,:]
        Z[i,:,:] = ar[1,:,:]
        U[i,:,:] = ar[2,:,:]
        V[i,:,:] = ar[3,:,:]
        P[i,:,:] = ar[4,:,:]
        H[i,:,:] = ar[5,:,:]
        PR[i,:,:] = ar[6,:,:]
        

    #number of timesteps per dat 
    pointsperday = 4

    #compute number of grid points and time steps from size of 3d matrix
    t,y,x=PR.shape
    gridpts=x*y
    tsteps=t

    #create y m d h vectors
    year = timesin.year
    month = timesin.month
    day = timesin.day
    hour = timesin.hour

    #create ID numbers for the grid points
    ID=1000000+np.linspace(1,gridpts,gridpts)

    #create matrices of x and y values
    X, Y = np.meshgrid(ar.x.values, ar.y.values)
    X=X.flatten(order='F')
    Y=Y.flatten(order='F')

    #elevation is static (doesn't change with time)
    elev=Z[1,:,:].flatten(order='F')

    #find number of grid points with <0 elevation. Note: this is related to the
    #subroutine met_data_check in the preprocess_code.f. that subroutine seems
    #to suggest that negative elevations are ok (say, death valley). But, the
    #code itself checks for negative elevations and stops execution is any
    #negatives are found.
    I = np.where(elev>=0)
    validgridpts=np.shape(I)[1]

    #remove data at points with neg elevations
    ID=ID[I]
    X=X[I]
    Y=Y[I]
    elev=elev[I]

    #we are now ready to begin our main loop over the time steps.
    fid= open(OUTpath,"w+")

    for j in range(tsteps):
        #first we write the number of grid points
        fid.write('{0:6d}\n'.format(validgridpts))

        #prep data matrix for this time step. First, grab the jth time slice
        Prtmp=PR[j,:,:].flatten(order='F')
        Htmp=H[j,:,:].flatten(order='F')
        Ptmp=P[j,:,:].flatten(order='F')
        Ttmp=T[j,:,:].flatten(order='F')
        Utmp=U[j,:,:].flatten(order='F')
        Vtmp=V[j,:,:].flatten(order='F')

        #remove data at points with neg elevations
        Prtmp=Prtmp[I]
        Htmp=Htmp[I]
        Ptmp=Ptmp[I]
        Ttmp=Ttmp[I]
        Utmp=Utmp[I]
        Vtmp=Vtmp[I]


        #convert precip rate to precip DEPTH (mm) during time interval
        Prtmp=Prtmp*24*3600/pointsperday

        #convert specific hum. to RH from Clausius-Clapeyron. T is still in K
        RHtmp=0.263*Ptmp*Htmp*(np.exp(17.67*(Ttmp-273.16)/(Ttmp-29.65)))**(-1)

        #compute wind speed
        SPDtmp=np.sqrt(Utmp**2+Vtmp**2)

        #compute wind direction. 0-360, with 0 being true north! 90 east, etc.
        DIRtmp=np.degrees(np.arctan2(Utmp,Vtmp))
        K=np.where(DIRtmp>=180)
        J=np.where(DIRtmp<180)
        DIRtmp[K]=DIRtmp[K]-180
        DIRtmp[J]=DIRtmp[J]+180

        #put T in C
        Ttmp=Ttmp-273.16

        for z in range(len(Prtmp)): 

            fid.write('{:5.0f}\t'.format(int(year[j]))+'{:5.0f}\t'.format(int(month[j]))+
                      '{:3.0f}\t'.format(int(day[j]))+'{:6.3f}\t'.format(hour[j])+
                      '{:9.0f}\t'.format(int(ID[z]))+'{:12.1f}\t'.format(X[z])+
                      '{:12.1f}\t'.format(Y[z])+'{:8.1f}\t'.format(elev[z])+
                      '{:9.2f}\t'.format(Ttmp[z])+'{:9.2f}\t'.format(RHtmp[z])+
                      '{:9.2f}\t'.format(SPDtmp[z])+'{:9.2f}\t'.format(DIRtmp[z])+
                      '{:9.2f}\n'.format(Prtmp[z]))
    fid.close()


# # RUN THE THANG


# set time parameters
stdt, eddt = set_dates(st_dt,ed_dt,date_flag)


# download GEE data
get_cfsv2(domain, TIFpath, stdt, eddt)


# fill in missing time slices or throw error if missing >4 slices
missing_slice_check(stdt, eddt, TIFpath)


# build SnowModel met file
MET2SM(TIFpath, OUTpath, stdt, eddt)


# delete directory with tif files 
get_ipython().system('rm -rf $TIFpath')

