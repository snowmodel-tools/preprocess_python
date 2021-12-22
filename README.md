# preprocess_python
This is a location for python code used to acquire and format data for a SnowModel run. 

Preprocess inputs:
* DEM data
* NLCD landcover data
* Meteorological data (hourly, 3-hourly, or 6-hourly):
  * elevation of met nodes [m]
  * surface air temperature [C]
  * wind speed [m/s]
  * precipitation [mm/timestep]
  * wind direction [deg]
    * NOTE: can use u-wind and v-wind to calculate wind speed and direction
  * RH [%]
    * NOTE: can use specific humidity and surface pressure to calculate RH
 * PRISM data
  * mean monthly temperature
  * mean monthly precipitaiton 

## The notebooks should be completed in the following order for each modeling domain:

####       1. DomainBounds_2json.ipynb
Notebook to create a JSON containing bounding parameters for all modeling domains. The output json file can be referenced by later scripts/notebooks using the requests package. Add new domains here. 

NOTE: This notebook needs to be run twice. On the first run, the user designates all domains, and includes the following details:
* Domain name
* Bounding box with latmax, latmin, lonmax, lonmin
* Start date
* End date
* Station projection (epsg:4326 in US)
* Model projection
 
After the CSO_domains.json is pushed to github, 02_GEE_topoveg.ipynb can be run. Ues the DEM or landcover ascii files for each domain to fill in the ncols, nrows, xll, yll values for each domain.

### 02_GEE_topoveg.ipynb

Notebook to pull dem and landcover data from GEE and prep ascii files for input into SnowModel. These ascii files will provide missing information (ncols, nrows, xll, yll) for the json. 


### 03_met_data.py 

Notebook to pull meteorological data from GEE and prep ascii files for input into SnowModel. 


### 04_Fetch_SNOTEL_CSO.ipynb

Notebook to get SNOTEL station data within a modeling domain to be used for the calibration.

### 05_Build_SnowModel_line_file.ipynb

Notebook to create a file to run SnowModel in line mode for the calibration. This notebook generates input files so that Snowmodel is only run at the cell(s) that correspond to station data. 

Files saved out:
* snowmodel_line_pts.dat
* .asc of line dem values
* .asc of line veg values
* .asc of line lat values
* .asc of line lon values

### 06_GET_MET.ipynb

Notebook to compare downscaled met outputs from MicroMet to station met data.

