"""coding: utf-8

4/25-22/2021

TODO:
- Pass the config yml file name as a run-time argument
- Add option to provide output directory path?
- Abstract the main function a bit (move chunks to functions)
- Implement tenacity-based (formerly retrying) checking for siphon/remote errors, and retrying

# # Extract and transform CFSv2 output subsets from NCEI THREDDS server
# 
# - 2021, 3/29: Successfully ran 30 days, 2019-11-01 to 2019-11-30. The server-side step (cell 11, which runs `get_subset_as_xrds`) took ~ 16 min ("CPU times: user 31 s, sys: 4.01 s, total: 35 s. Wall time: 15min 45s"). So, that's roughly 0.5 min per model day.
# - 2021, 3/22: It's working again!?
# - 2021: 3/10,7, 2/17,12
# - 2020: 12/16, 11-16,5. https://github.com/emiliom/
# - Run with conda environment https://github.com/snowmodel-tools/postprocess_python/blob/master/snowmodelaws_env.yml
# 
# Define and execute a remote request for CFSv2 model output. The user specifies the geographic bounding box (in lat-lon coordinates), the output projection (as EPSG code) and cell resolution (in meters), and the date range.
# 
# The overall approach may be grouped into several overarching steps:
# 1. Define the query parameters
# 2. Loop over each day in the requested day range
# 3. Dynamically construct the url for a siphon `TDSCatalog` object corresponding to a specific date (a THREDDS directory catalog)
# 4. Within the selected date, loop through the 4 model time steps ("cycles"), dynamically constructing the `grib2` file name and file url.
# 4. Execute subset and data access requests, one grib2 file at a time and returning a collection (list) of xarray datasets; one file corresponds to one timestep.
# 5. Concatenate the individual xarray datasets and clean up the result, simplifying the structure and adding conformant projection information.
# 6. Reproject and resample to the desired UTM projection.
# 7. Export to netcdf file
# 
# - Currently the code uses this model product: [CFSv2 Operational Analysis - 6-Hourly Surface and Radiative Fluxes (FLX)](https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/climate-forecast-system-version2-cfsv2#CFSv2%20Operational%20Analysis).
# - For additional background on CFSv2 see https://github.com/snowmodel-tools/preprocess_python/issues/17
# - CSO GEE JavaScript code used: https://github.com/snowmodel-tools/preprocess_javascript/blob/master/define_SM_inputs.js

# ## TO-DOs
# 3/22/2021
# 
# - `snowmodelaws` env is on Py3.6 b/c of ulmo! Once I issue a new release, we'll be able to move to 3.8 or 3.9
# - Add error catching to `get_subset_as_xrds`, to log and retry datasets (files) that ran into a failure during processing. The failure is typically on the NCEI TDS server end. This is the latest error, that's likely to be common:
# 
# ```python
# ~/miniconda/envs/snowmodelaws/lib/python3.6/site-packages/siphon/http_util.py in get(self, path, params)
#     485                                      'Server Error ({1:d}: {2})'.format(resp.request.url,
#     486                                                                         resp.status_code,
# --> 487                                                                         text))
#     488         return resp
#     489 
# 
# HTTPError: Error accessing https://www.ncei.noaa.gov/thredds/ncss/model-cfs_v2_anl_6h_flxf/2019/201911/20191126/cdas1.t00z.sfluxgrbf06.grib2/dataset.xml
# Server Error (502: Proxy Error)
# ```
"""

from datetime import datetime, timedelta

import yaml
import dateutil
from dateutil.parser import parse as duparse
import pandas as pd
from siphon.catalog import TDSCatalog
import xarray as xr
from xarray.backends import NetCDF4DataStore
from pyproj import CRS
from rasterio.warp import Resampling
import rioxarray


# CFSv2 variables used
VARIABLE_NAMES = [
    'Temperature_height_above_ground',
    'Geopotential_height_surface',
    'u-component_of_wind_height_above_ground',
    'v-component_of_wind_height_above_ground',
    'Pressure_surface',
    'Specific_humidity_height_above_ground',
    'Precipitation_rate_surface_6_Hour_Average',
    'Downward_Long-Wave_Radp_Flux_surface',
    'Downward_Short-Wave_Radiation_Flux_surface_6_Hour_Average',
]

def get_subset_as_xrds(tds_ds, var_lst, bbox):
    """ Request a geographical and variable subset for the specified dataset (one time step).
    Returns an xarray dataset.
    """
    # Interface with the data through the NetCDF Subset Service (NCSS)
    ncss = tds_ds.subset()

    # Create an NCSS query with our desired specifications
    query = ncss.query()
    query.all_times()
    query.variables(*var_lst)
    query.lonlat_box(**bbox)
    query.accept('netcdf')  # or netcdf4?

    # Use the query to obtain NetCDF data
    data = ncss.get_data(query)
    data_ds = xr.open_dataset(NetCDF4DataStore(data))
    
    return data_ds


def concatds_cleanup(ds_lst, crs, grid_mapping_name='crs'):
    """ Concatenate and clean up list of xarray datasets.
    Write in and clean up the CRS / `grid_mapping` attributes. Shift longitude values from 0 to 360 to -180 to +180.
    In order to be able to use `rio.reproject`, and to generate a less complex xarray dataset. Specifically:
    - Insert the information (height) into attributes in the corresponding variables, 
      from the `height_above_ground` and `height_above_ground1` variables
    - Remove (squeeze) those dimensions from the variables currently using them
    - Remove the coordinate variables themselves
    """
    
    # Assumes the data variables returned use two identical time dimensions, 'time' and 'time1'
    # Concatenate on time1, drop the orphaned time dim, then rename time1 to time
    ds = xr.concat(ds_lst, dim='time1').drop_vars('GaussLatLon_Projection')
    ds = ds.drop_dims('time').rename({'time1':'time'})
    
    ds.rio.write_crs(crs, grid_mapping_name=grid_mapping_name, inplace=True);

    # longitude rewrap
    lon_attrs = ds.lon.attrs
    ds.coords['lon'] = (ds.coords['lon'] + 180) % 360 - 180
    ds = ds.sortby(ds.lon)
    ds.lon.attrs = lon_attrs

    # height_above_ground simplification
    height_dim_names = ['height_above_ground', 'height_above_ground1']
    for dim_name in height_dim_names:
        for data_var_name in ds.data_vars:
            data_var = ds[data_var_name]
            if dim_name in ds[data_var_name].dims:
                data_var.attrs['height_above_ground'] = "{:.1f}m".format(ds[dim_name].values[0])
                ds[data_var_name] = data_var.squeeze(dim_name, drop=True)

    # TODO: Update the variable coords attribute, for consistency
    
    ds = ds.squeeze(height_dim_names, drop=True)
    
    return ds


def main():
    # ## Data request and filtering
    # Includes filtering and selection scheme.

    tds_base_url = "https://www.ncei.noaa.gov/thredds/catalog/model-cfs_v2_anl_6h_flxf"

    # Read run-time configuration from yaml file
    with open("CFSv2_config_WY.yml", 'r') as f:
        config = yaml.safe_load(f)

    # Geo
    geobbox = dict(
        south=config['geo']['bbox']['south'],
        west=config['geo']['bbox']['west'],
        north=config['geo']['bbox']['north'],
        east=config['geo']['bbox']['east'],
    )

    to_crs = f"epsg:{config['geo']['to_epsg_crs']}"
    to_resolution = config['geo']['to_epsg_crs']

    # Output file name
    nc_export_fname = config['nc_export_fname']

    # Create a sligthly larger box for the data query, relative to geobbox
    buffered_geobbox = geobbox.copy()
    buffered_geobbox['south'] = geobbox['south'] - config['geo']['bbox']['buffer_north']
    buffered_geobbox['north'] = geobbox['north'] + config['geo']['bbox']['buffer_north']
    buffered_geobbox['west'] = geobbox['west'] - config['geo']['bbox']['buffer_east']
    buffered_geobbox['east'] = geobbox['east'] + config['geo']['bbox']['buffer_east']
    # ----------------------------------------

    # ### Construct and execute the query
    # 
    # Step through days in the date range, setting the `TDSCatalog` for each day, then stepping 
    # through the 4 model time steps ("model cycles"). The path to the `TDSCatalog` is constructed 
    # dynamically based on the date; the path to each `grib2` dataset is also constructed dynamically.

    # https://www.ncei.noaa.gov/thredds/catalog/model-cfs_v2_anl_6h_flxf/2018/201808/20180831/catalog.html?dataset=cfs_v2_anl_6h_flxf/2018/201808/20180831/cdas1.t00z.sfluxgrbf06.grib2
    # 
    # https://www.ncei.noaa.gov/thredds/ncss/model-cfs_v2_anl_6h_flxf/2018/201808/20180831/cdas1.t00z.sfluxgrbf06.grib2/dataset.xml

    data_ds_lst = []
    for day in pd.date_range(start=config['start_date'], end=config['end_date']):
        print(day)
        cfsv2_cat = TDSCatalog(f"{tds_base_url}/{day:%Y}/{day:%Y%m}/{day:%Y%m%d}/catalog.xml")
        # NOTE: THIS SCHEME WILL PULL ENTIRE DAYS' WORTH OF DATA; IT DOES NOT REFINE THE QUERY BY TIME
        # Step through time steps (UTC times)
        fts = 6
        for t in [0, 6, 12, 18]:
            grb2_fname = f"cdas1.t{t:02d}z.sfluxgrbf{fts:02d}.grib2"
            print(f"    {grb2_fname}")
            grb2_tds = cfsv2_cat.datasets[grb2_fname]
            xrds = get_subset_as_xrds(grb2_tds, VARIABLE_NAMES, buffered_geobbox)
            for v in VARIABLE_NAMES:
                # Most of the variables are assigned time1, so it should be more efficient 
                # to rename on 'time'
                if 'time' in xrds[v].dims:
                    xrds[v] = xrds[v].rename({'time': 'time1'})
            
            data_ds_lst.append(xrds)

    # The concatenation step assumes the datasets in `data_ds_lst` were added in 
    # ascending chronological order, which has been the case so far.
    data_ds = concatds_cleanup(data_ds_lst, crs='epsg:4326')

    # ### Explore the resulting dataset
    def print_cellsizes(dimcoord, i):
        return dimcoord.values[i+1] - dimcoord.values[i]

    # print_cellsizes(data_ds.lon, 0), print_cellsizes(data_ds.lat, 0), print_cellsizes(data_ds.lat, 1)
    # ----------------------------------------

    # ## Reproject and resample dataset
    x_name = 'easting'
    y_name = 'northing'

    # resampling=Resampling.bilinear
    data_reproj_ds = data_ds.rio.reproject(to_crs, resolution=to_resolution, resampling=Resampling.nearest)
    data_reproj_ds = (
        data_reproj_ds
        .rio.set_spatial_dims('x', 'y')
        .rename({"x": x_name, "y": y_name}) 
    )

    data_reproj_ds[x_name].attrs['long_name'] = 'Easting'
    data_reproj_ds[y_name].attrs['long_name'] = 'Northing'

    # NOTE: When the reprojected resolution was set to 20km (possibly coarser than the 
    # original resolution), some cells were set to the `_FillValue` but xarray didn't seem 
    # to properly apply the `_FillValue`.
    # ----------------------------------------

    # ## Export to netcdf
    data_reproj_ds.to_netcdf(nc_export_fname)


if __name__ == '__main__':
    main()
