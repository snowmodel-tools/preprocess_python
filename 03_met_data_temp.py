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
domain = 'ST'
print(domain)
#path to CSO domains
domains_resp = requests.get("https://raw.githubusercontent.com/snowmodel-tools/preprocess_python/master/CSO_domains.json")
domains = domains_resp.json()

# PATHS
# path to temporary folder to store tif files from gee
TIFpath = 'GEE_Downloads_tmp_ID/'
# path to where you want your output met .dat fime
OUTpath = '/nfs/attic/dfh/Aragon2/CSOdmn/ST/mm_'+domain+'_2015-2018.dat'

# TIME
# choose if want to set 'manual' or 'auto' date 
date_flag = 'manual'
# If you choose 'manual' set your dates below  
# This will start on the 'begin' date at 0:00 and the last iteration will 
# be on the day before the 'end' date below.
st_dt = '2015-10-01'#domains[domain]['st']
ed_dt = '2018-08-31'#domains[domain]['ed']
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


# build SnowModel met file
MET2SM(TIFpath, OUTpath, stdt, eddt)


# delete directory with tif files 
get_ipython().system('rm -rf $TIFpath')

