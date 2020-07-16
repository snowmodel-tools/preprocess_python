# preprocess_python
This is a location for python code used to acquire and format data for a SnowModel run. 

####       DomainBounds_2json 
Notebook to create a JSON containing bounding parameters for all modeling domains.

## The notebooks should be completed in the following order for each modeling domain:

####       1. Fetch_SNOTEL_CSO.ipynb
Notebook to build json to save domain information. Add new domains here. 

####       2. Fetch_SNOTEL_CSO.ipynb
Notebook to get SNOTEL station data within a modeling domain.

####       3. Build_SnowModel_line_file.ipynb
Notebook to create file to run SnowModel in line mode for calibration.


## (optional notebooks)

####       Baseline_par2json.ipynb
Notebook to extract .par values and save them to a JSON. 

####       Explore_SNOTEL.ipynb
Notebook to examine SNOTEL stations within a modeling domain.

####       Extract_data_from_raster.ipynb
Notebook to extract raster values at SNOTEL locations.

####       Set_baseSM.ipynb
Notebook to set the SnowModel .par file to baseline values.