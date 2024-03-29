#!/usr/bin/env python
# coding: utf-8

# In[1]:


# create a JSON with bounding parameters for each modeling domain
# get ncols, nrows, ll corner coordinates, and cell size from domain ascii files. 


# In[10]:


WY = {
    "name": 'WY',
    "Bbox": {'latmax' : 44.582480,'latmin' : 42.363116,'lonmax': -109.477849,'lonmin': -111.155208,},
    "st": "2011-09-01", 
    "ed": "2016-09-30",
    "stn_proj": 'epsg:4326',
    "mod_proj": 'epsg:32612',
    "ncols": '1382',
    "nrows": '2476',
    "xll": '487200',
    "yll": '4690100', 
    "cellsize": '100'
}

OR = {
    "name": 'OR',
    "Bbox": {'latmax' : 44.74,'latmin' : 43.645,'lonmax': -121.11,'lonmin': -122.11,},
    "st": "2011-09-01", 
    "ed": "2016-09-30",
    "stn_proj": 'epsg:4326',
    "mod_proj": 'epsg:32610', 
    "ncols": '821',
    "nrows": '1231',
    "xll": '570400',
    "yll": '4832800', 
    "cellsize": '100'
}

UT = {
    "name": 'UT',
    "Bbox": {'latmax' : 40.7119,'latmin' : 40.4407,'lonmax': -111.461,'lonmin': -111.7986,},
    "st": "2011-09-01", 
    "ed": "2016-09-30",
    "stn_proj": 'epsg:4326',
    "mod_proj": 'epsg:32612',
    "ncols": '960',
    "nrows": '1012',
    "xll": '432270',
    "yll": '4476750', 
    "cellsize": '30'
}

WA = {
    "name": 'WA',
    "Bbox": {'latmax' : 48.9,'latmin' : 47.8299,'lonmax': -119.769,'lonmin': -121.1,},
    "st": "2011-10-01", 
    "ed": "2016-09-30",
    "stn_proj": 'epsg:4326',
    "mod_proj": 'epsg:32610',
    "ncols": '2052',
    "nrows": '2446',
    "xll": '639200',
    "yll": '5299100', 
    "cellsize": '50'
}


CA = {
    "name": 'CA',
    "Bbox": {'latmax' : 39.7,'latmin' : 38.4,'lonmax': -119.66,'lonmin': -120.85,},
    "st": "2011-09-01", 
    "ed": "2016-09-30",
    "stn_proj": 'epsg:4326',
    "mod_proj": 'epsg:32611',
    "ncols": '2164',
    "nrows": '2962',
    "xll": '163750',
    "yll": '4253500', 
    "cellsize": '50'
}

HJA = {
    "name": 'HJA',
    "Bbox": {'latmax' : 44.28115,'latmin' : 44.19873,'lonmax': -122.10189,'lonmin': -122.26163,},
    "st": "2011-09-01", 
    "ed": "2016-09-30",
    "stn_proj": 'epsg:4326',
    "mod_proj": 'epsg:32610',
    "ncols": '429',
    "nrows": '310',    
    "xll": '558900',
    "yll": '4894200', 
    "cellsize": '30'
}

CO = {
    "name": 'CO',
    "Bbox": {'latmax' : 40.8,'latmin' : 37.1,'lonmax': -104.8,'lonmin': -108.5,},
    "st": "2011-09-01", 
    "ed": "2016-09-30",
    "stn_proj": 'epsg:4326',
    "mod_proj": 'epsg:32613',
    "ncols": '3289',
    "nrows": '4166',
    "xll": '188900',
    "yll": '4105900', 
    "cellsize": '100'
}

ID = {
    "name": 'ID',
    "Bbox": {'latmax' : 44.483224,'latmin' : 43.559833,'lonmax': -114.934861,'lonmin': -116.462681,},
    "st": "2011-09-01", 
    "ed": "2016-09-30",
    "stn_proj": 'epsg:4326',
    "mod_proj": 'epsg:32612',
    "ncols": '2567',
    "nrows": '2192',
    "xll": '58750',
    "yll": '4830500', 
    "cellsize": '50'
}

AK = {
    "name": 'AK',
    "Bbox": {'latmax' : 61.0179627,'latmin' : 60.67770498,'lonmax': -148.7578855597,'lonmin': -149.38460973,},
    "st": "2011-09-01", 
    "ed": "2016-09-30",
    "stn_proj": 'epsg:4326',
    "mod_proj": 'epsg:32607',
    "ncols": '',
    "nrows": '',
    "xll": '',
    "yll": '', 
    "cellsize": '30'
}

CO_N = {
    "name": 'CO_N',
    "Bbox": {'latmax' : 40.194217,'latmin' : 38.370328,'lonmax': -105.407231,'lonmin': -107.538579,},
    "st": "2011-09-01", 
    "ed": "2016-09-30",
    "stn_proj": 'epsg:4326',
    "mod_proj": 'epsg:32613',
    "ncols": '1872',
    "nrows": '2056',
    "xll": '278200',
    "yll": '4246900', 
    "cellsize": '100'
}

CO_S = {
    "name": 'CO_S',
    "Bbox": {'latmax' : 38.038433,'latmin' : 37.304722,'lonmax': -106.463801,'lonmin': -108.106257,},
    "st": "2011-09-01", 
    "ed": "2016-09-30",
    "stn_proj": 'epsg:4326',
    "mod_proj": 'epsg:32613',
    "ncols": '1470',
    "nrows": '851',
    "xll": '224600',
    "yll": '4129600', 
    "cellsize": '100'
}

NH = {
    "name": 'NH',
    "Bbox": {'latmax' : 44.328173,'latmin' : 43.901758,'lonmax': -71.009589,'lonmin': -71.754594,},
    "st": "2011-10-01", 
    "ed": "2016-09-30",
    "stn_proj": 'epsg:4326',
    "mod_proj": 'epsg:32619',
    "ncols": '2034',
    "nrows": '1637',
    "xll": '278760',
    "yll": '4862910', 
    "cellsize": '30'
}

WA_SQ = {
    "name": 'WA_SQ',
    "Bbox": {'latmax' : 47.550,'latmin' : 47.189,'lonmax': -121.183,'lonmin': -121.740,},
    "st": "2011-10-01",
    "ed": "2016-09-30",
    "stn_proj": 'epsg:4326',
    "mod_proj": 'epsg:32610',
    "ncols": '1429',
    "nrows": '1366',
    "xll": '594780',
    "yll": '5226930',
    "cellsize": '30'
}

CSO_domains = {
    "WY": WY,
    "OR": OR,
    "UT": UT,
    "WA": WA,
    "CA":CA,
    "HJA":HJA,
    "CO":CO,
    "CO_N":CO_N,
    "CO_S":CO_S,
    "AK":AK,
    "ID":ID,
    "NH":NH,
    "WA_SQ":WA_SQ
}


# In[11]:


CSO_domains.keys()


# In[12]:


import json
json = json.dumps(CSO_domains)
f = open('/nfs/attic/dfh/Aragon2/Notebooks/preprocess_python/CSO_domains.json',"w")
f.write(json)
f.close()


# In[3]:


# #save baseline .par file
# import json

# json = json.dumps(CSO_domains)
# f = open('/nfs/attic/dfh/Aragon2/CSOdata/CSO_domains.json',"w")
# f.write(json)
# f.close()

