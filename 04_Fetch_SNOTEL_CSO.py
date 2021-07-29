#!/usr/bin/env python
# coding: utf-8

# # Fetch and export SNOTEL sites and daily time series data

# In[1]:


from collections import OrderedDict
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely import geometry as sgeom
import ulmo
import json
import requests


# In[2]:


#########################################################################
############################ USER INPUTS ################################
#########################################################################

# DOMAIN
# choose the modeling domain
domain = 'WY'

# PATHS
dataPath = '/nfs/attic/dfh/Aragon2/CSOdmn/'+domain+'/'


# DATES
st_dt = '2011-09-01'
ed_dt = '2016-09-30'
#########################################################################


# In[3]:


# functions to get SNOTEL stations as geodataframe
def sites_asgdf(ulmo_getsites, stn_proj):
    """ Convert ulmo.cuahsi.wof.get_sites response into a point GeoDataframe
    """
    
    # Note: Found one SNOTEL site that was missing the location key
    sites_df = pd.DataFrame.from_records([
        OrderedDict(code=s['code'], 
        longitude=float(s['location']['longitude']), 
        latitude=float(s['location']['latitude']), 
        name=s['name'], 
        elevation_m=s['elevation_m'])
        for _,s in ulmo_getsites.items()
        if 'location' in s
    ])

    sites_gdf = gpd.GeoDataFrame(
        sites_df, 
        geometry=gpd.points_from_xy(sites_df['longitude'], sites_df['latitude']),
        crs=stn_proj
    )
    return sites_gdf

def get_snotel_stns(domain):
    
    #path to CSO domains
    domains_resp = requests.get("https://raw.githubusercontent.com/snowmodel-tools/preprocess_python/master/CSO_domains.json")
    domains = domains_resp.json()

    #Snotel bounding box
    Bbox = domains[domain]['Bbox']

    # Snotel projection
    stn_proj = domains[domain]['stn_proj']
    # model projection
    mod_proj = domains[domain]['mod_proj']

    # Convert the bounding box dictionary to a shapely Polygon geometry using sgeom.box
    box_sgeom = sgeom.box(Bbox['lonmin'], Bbox['latmin'], Bbox['lonmax'], Bbox['latmax'])
    box_gdf = gpd.GeoDataFrame(geometry=[box_sgeom], crs=stn_proj)
    
    # WaterML/WOF WSDL endpoint url 
    wsdlurl = "http://hydroportal.cuahsi.org/Snotel/cuahsi_1_1.asmx?WSDL"

    # get dictionary of snotel sites 
    sites = ulmo.cuahsi.wof.get_sites(wsdlurl,user_cache=True)

    #turn sites to geodataframe 
    snotel_gdf = sites_asgdf(sites,stn_proj)
    
    #clip snotel sites to domain bounding box
    gdf = gpd.sjoin(snotel_gdf, box_gdf, how="inner")
    gdf.drop(columns='index_right', inplace=True)
    gdf.reset_index(drop=True, inplace=True)

    #add columns with projected coordinates 
    CSO_proj = gdf.to_crs(mod_proj)
    gdf['easting'] = CSO_proj.geometry.x
    gdf['northing'] = CSO_proj.geometry.y
    
    return gdf


# In[5]:


def fetch(sitecode, variablecode, start_date, end_date):
    print(sitecode, variablecode, start_date, end_date)
    values_df = None
    wsdlurl = "http://hydroportal.cuahsi.org/Snotel/cuahsi_1_1.asmx?WSDL"
    try:
        #Request data from the server
        site_values = ulmo.cuahsi.wof.get_values(
            wsdlurl, 'SNOTEL:'+sitecode, variablecode, start=start_date, end=end_date
        )
        #Convert to a Pandas DataFrame   
        values_df = pd.DataFrame.from_dict(site_values['values'])
        #Parse the datetime values to Pandas Timestamp objects
        values_df['datetime'] = pd.to_datetime(values_df['datetime'])
        #Set the DataFrame index to the Timestamps
        values_df.set_index('datetime', inplace=True)
        #Convert values to float and replace -9999 nodata values with NaN
        values_df['value'] = pd.to_numeric(values_df['value']).replace(-9999, np.nan)
        #Remove any records flagged with lower quality
        values_df = values_df[values_df['quality_control_level_code'] == '1']
    except:
        print("Unable to fetch %s" % variablecode)
    
    return values_df


# In[6]:


# returns daily timeseries
# https://www.wcc.nrcs.usda.gov/web_service/AWDB_Web_Service_Reference.htm#commonlyUsedElementCodes
# 'WTEQ': swe [in]
# 'SNWD': snow depth [in]
# 'PRCP': precipitation increment [in]
# 'PREC': precipitation accumulation [in]
# 'TAVG': average air temp [F]
# 'TMIN': minimum air temp [F]
# 'TMAX': maximum air temp [F]
# 'TOBS': observered air temp [F]
def get_snotel_data(gdf,sd_dt, ed_dt,var,units='metric'):
    '''
    gdf - pandas geodataframe of SNOTEL sites
    st_dt - start date string 'yyyy-mm-dd'
    ed_dt - end date string 'yyyy-mm-dd'
    var - snotel variable of interest 
    units - 'metric' (default) or 'imperial'
    '''
    stn_data = pd.DataFrame(index=pd.date_range(start=st_dt, end=ed_dt))
    

    for sitecode in gdf.code:
        try:
            data = fetch(sitecode,'SNOTEL:'+var+'_D', start_date=st_dt, end_date=ed_dt)
            #check for nan values
            if len(data.value[np.isnan(data.value)]) > 0:
                #check if more than 10% of data is missing
                if len(data.value[np.isnan(data.value)])/len(data) > .1:
                    print('More than 10% of days missing')
                    gdf.drop(gdf.loc[gdf['code']==sitecode].index, inplace=True)
                    continue
            stn_data[sitecode] = data.value
        except:
            gdf.drop(gdf.loc[gdf['code']==sitecode].index, inplace=True)     
    
    gdf.reset_index(drop=True, inplace=True)
    if units == 'metric':
        if (var == 'WTEQ') |(var == 'SNWD') |(var == 'PRCP') |(var == 'PREC'):
            #convert SNOTEL units[in] to [m]
            for sitecode in gdf.code:
                stn_data[sitecode] = 0.0254 * stn_data[sitecode]
        elif (var == 'TAVG') |(var == 'TMIN') |(var == 'TMAX') |(var == 'TOBS'):
            #convert SNOTEL units[F] to [C]
            for sitecode in gdf.code:
                stn_data[sitecode] = (stn_data[sitecode] - 32) * 5/9
    return gdf, stn_data


# # Execute Functions

# In[6]:


#get geodataframe of all SNOTEL sites in the domain
snotel_gdf = get_snotel_stns(domain)

#get xy coordinates of stations in gdf 

#get SWE timeseries 
domain_gdf, swe = get_snotel_data(snotel_gdf,st_dt,ed_dt,'WTEQ')
#get snow depth timeseries 
domain_gdf, hs = get_snotel_data(snotel_gdf,st_dt,ed_dt,'SNWD')
#get precipitation timeseries 
domain_gdf, pr = get_snotel_data(snotel_gdf,st_dt,ed_dt,'PRCP')
#get av temp timeseries 
domain_gdf, tav = get_snotel_data(snotel_gdf,st_dt, ed_dt,'TAVG')
#get min temp timeseries 
domain_gdf, tmn = get_snotel_data(snotel_gdf,st_dt, ed_dt,'TMIN')
#get max temp timeseries 
domain_gdf, tmx = get_snotel_data(snotel_gdf,st_dt, ed_dt,'TMAX')


# # Save Data

# In[9]:


# save geojson
out = dataPath + 'CSO_SNOTEL_sites.geojson'
domain_gdf.to_file(out, driver='GeoJSON')

#save swe
out = dataPath + 'SNOTEL_data_SWEDmeters'+st_dt+'_'+ed_dt+'.csv'
swe.to_csv(out)

#save hs
out = dataPath + 'SNOTEL_data_HSmeters'+st_dt+'_'+ed_dt+'.csv'
hs.to_csv(out)

#save pr
out = dataPath + 'SNOTEL_data_PRmeters'+st_dt+'_'+ed_dt+'.csv'
pr.to_csv(out)

#save tav
out = dataPath + 'SNOTEL_data_TAVGcelsius'+st_dt+'_'+ed_dt+'.csv'
tav.to_csv(out)

#save tmn
out = dataPath + 'SNOTEL_data_TMINcelsius'+st_dt+'_'+ed_dt+'.csv'
tmn.to_csv(out)

#save tmx
out = dataPath + 'SNOTEL_data_TMAXcelsius'+st_dt+'_'+ed_dt+'.csv'
tmx.to_csv(out)

