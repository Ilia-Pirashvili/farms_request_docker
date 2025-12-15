# #!/bin/env python

import pandas as pd
import xarray as xr
import numpy as np
import os
import pyproj
import rasterio
from rasterio.transform import from_bounds
from rasterio.features import geometry_mask
from pyproj import CRS, Transformer
import shapely
from shapely.geometry import Polygon, Point
from shapely.ops import transform as shapely_transform
from datetime import date

from app.date_functions import sorted_dates, DateRange#, complement_in_date_range


def find_dates_in_NetCDF(NetCDF_path: str, variable_name: str, attr_name: str ="no_data_time") -> list[date]:
  """
  Returns the datetime.dates()'s which have already been addressed in the NetCDF file.
    By addressed, we mean the dates where either data is stored,
    or the dates stored in the attribute where no data is expected.
  ---
  Input:
    - NetCDF_path  : String. Path to the NetCDF file.
    - variable_name: String. The variable name for which we will look for the missing dates.
    - attr_name    : String. Name of the attribute storing dates where no data is expected.
  Output:
    - List of datetime.dates()'s.
  """
  
  if not os.path.exists(NetCDF_path):
    return []
  
  # Use chunks='auto' to improve reading speed and memory usage
  # Use decode_times=False to avoid automatic datetime conversion which we'll handle more efficiently
  with xr.open_dataset(NetCDF_path, engine="h5netcdf") as ds:
    data_dates = set()
    if variable_name in ds.data_vars and 'time' in ds[variable_name].coords:
      # Get the time values directly as numpy array first
      time_values = ds[variable_name].time.values
      time_values = [pd.to_datetime(date).date() for date in time_values]
      # Process dates in batches for better performance
      batch_size = 1000
      for i in range(0, len(time_values), batch_size):
          batch = time_values[i:i+batch_size]
          # Convert in batch which is more efficient
          converted_dates = pd.to_datetime(batch).date  # pyright: ignore
          data_dates.update(converted_dates)
    elif 'time' in ds.coords:
      time_values = ds.time.values
      time_values = [pd.to_datetime(date).date() for date in time_values]
      batch_size = 1000
      for i in range(0, len(time_values), batch_size):
        batch = time_values[i:i+batch_size]
        # Convert in batch which is more efficient
        converted_dates = pd.to_datetime(batch).date   # pyright: ignore
        data_dates.update(converted_dates)
    
    no_data_dates = set()
    if attr_name in ds.attrs and ds.attrs[attr_name]:
      no_data_dates = set(pd.to_datetime(ds.attrs[attr_name]).date)
    
    return list(data_dates.union(no_data_dates))

def force_save_xr_ds(dataset: xr.Dataset, NetCDF_path: str):
  """
  Safely (?) (or the opposite :)) saves a dataset as a NetCDF.
  This is needed since if the NetCDF file is already open, you can not directly save to it.
  ---
  Input:
    - dataset: xr.dataset.
    - NetCDF_file: string. Path to a NetCDF file
  Output:
    - Void
  """
  try:
    dataset.to_netcdf(NetCDF_path, engine="h5netcdf")
  except:
    temp_path = NetCDF_path + ".tmp"
    dataset.to_netcdf(temp_path, engine="h5netcdf")
    os.replace(temp_path, NetCDF_path)

def add_missing_dates_as_attr(NetCDF_path: str, missing_dates: list[date], attr_name="no_data_time"):
  """
  Adds a list of missing dates as an attribute to a NetCDF file.
  ---
  Input:
    - Net_CDF_file: String. Path to a NetCDF file
    - missing_dates: list of datetime.date()'s
    - attr_name: String. Name of attribute where the missing dates are stored. Should generally not be changed from the default...
  Output:
    - Void
  """
  str_dates = [d.strftime("%Y-%m-%d") for d in missing_dates]
  with xr.open_dataset(NetCDF_path, engine="h5netcdf") as ds:
    ds = ds.load()
    if attr_name in ds.attrs:
      old_dates = ds.attrs[attr_name]
    else:
      old_dates = []
    old_dates.extend(str_dates)
    old_dates = list(set(old_dates))
    ds.attrs[attr_name] = sorted_dates(old_dates)
    force_save_xr_ds(ds, NetCDF_path)

def restrict_NetCDF_to_polygon(
  NetCDF_file: str,
  variable_name: str,
  polygon: Polygon,
  date_range: DateRange = DateRange(),
  polygon_crs: int | str = 4326,
  restrict_to_bbox: bool = True
) -> xr.Dataset:
  """
  Restricts the contents (xr.Dataset) of a NetCDF to a shapely polygon.
    - The dataset is expected to come from a tiff (or be with coordinates x and y and a variable)
    - NOTE, it will not overwrite the NetCDF, just return the restricted dataset.
  ---
  Input:
    - NetCDF_file: string. Path to file
    - variable_name: string. The variable form which we derive the crs. Otherwise, all data is returned.
    - polygon: shapely polygon we restrict to.
    - date_range: Pair of string of the type ["2014-02-02", "2015-09-07"] or the DateRange class.
    - polygon_crs: integer or string. The crs of the polygon we restrict to. Assumed to be 4326 (long/lat).
    - restrict_to_bbox: boolean. If True, restrict to the bounding box of the polygon. 
                        If False, restrict to the exact polygon shape. Default is True.
  Output:
    - The restricted dataset.
  """
  with xr.open_dataset(NetCDF_file, engine="h5netcdf") as ds:
    if date_range:
      start_date = pd.to_datetime(date_range[0]).date()
      end_date   = pd.to_datetime(date_range[1]).date()
      ds         = ds.sel(time=slice(start_date, end_date))
    try:
      NetCDF_crs = ds[variable_name].crs
    except:
      NetCDF_crs = ds.crs
    project    = pyproj.Transformer.from_crs(polygon_crs, NetCDF_crs, always_xy=True).transform
    polygon    = shapely_transform(project, polygon)   #shapley_transform = shapely.transform
    minx, miny, maxx, maxy = polygon.bounds
    ds = ds.sel(y=slice(maxy, miny), x=slice(minx, maxx))
    
    # If we're only restricting to bbox, return the dataset now
    if restrict_to_bbox:
      return ds
    
    # Otherwise, further restrict to the actual polygon shape
    # Create a mask for points inside the polygon
    x_coords = ds.x.values
    y_coords = ds.y.values
    mask = np.zeros((len(y_coords), len(x_coords)), dtype=bool)
    
    # Create points for each coordinate pair and check if they're in the polygon
    for i, y in enumerate(y_coords):
      for j, x in enumerate(x_coords):
        point = Point(x, y)
        mask[i, j] = polygon.contains(point)
    
    # Apply the mask to the dataset
    # First convert mask to xarray DataArray with proper coordinates
    mask_da = xr.DataArray(mask, coords=[('y', y_coords), ('x', x_coords)])
    
    # Apply mask to all data variables
    for var in ds.data_vars:
      # Handle different dimensions - mask should be applied only to x,y dimensions
      if 'x' in ds[var].dims and 'y' in ds[var].dims:
        # Find x,y dimensions in the variable
#        dims = ds[var].dims
        # Apply mask - need to broadcast appropriately
        ds[var] = ds[var].where(mask_da)
    
  return ds

def sync_ds_to_NetCDF(ds: xr.Dataset, NetCDF_path: str) -> xr.Dataset:
  """
  Adds the new data to a NetCDF file (or creates a new one).
  Handles minor coordinate differences by realigning grids when necessary.
  If the new dataset is a subset of the existing data, returns the subset without modifying the file.
  NOTE: For time, it simply drops duplicates. It is not making check to ensure old data is preferred!
  ---
  Input:
    - ds: xr.Dataset. The new dataset information to be synced with the NetCDF file
    - NetCDF_path: string. The path to the NetCDF file to be updated.
  Output:
    - updated xr.Dataset
  """
  if not os.path.exists(NetCDF_path):
    if len(ds.x.values) > 0:
      updated_ds = ds
    else:
      raise ValueError(f"There seems to be no data in the dataset. Did not create {NetCDF_path}."
                           "If you used a combinable model, make sure that your polygon actually lies in the multipolygon.")
  else:
    with xr.open_dataset(NetCDF_path, engine="h5netcdf") as old_ds:
      # Check if new dataset is a subset of the existing one
      def is_subset(new_coord, old_coord):
        """Check if new_coord is a subset of old_coord with tolerance
        Returns: 0 = not subset, 1 = proper subset, 2 = same set"""
        if len(new_coord) > len(old_coord):
          return 0
        
        # Handle different coordinate types
        if new_coord.dtype.kind == 'M':  # datetime coordinates
          # For datetime, use exact matching
          new_vals = set(new_coord.values)
          old_vals = set(old_coord.values)
          if new_vals.issubset(old_vals):
            return 2 if new_vals == old_vals else 1
          return 0
        else:
          # For numeric coordinates, use tolerance
          # Check if it's essentially the same set (within 1-2 pixel tolerance)
          if abs(len(new_coord) - len(old_coord)) <= 2:
            # Check if coordinates are very close (same grid, minor differences)
            matches = 0
            for new_val in new_coord.values:
              if any(abs(new_val - old_val) < 1e-10 for old_val in old_coord.values):
                matches += 1
            
            # If most coordinates match and length difference is small, consider it same set
            if matches >= len(new_coord) and abs(len(new_coord) - len(old_coord)) <= 2:
              return 2
          
          # Check if it's a proper subset
          for new_val in new_coord.values:
            if not any(abs(new_val - old_val) < 1e-10 for old_val in old_coord.values):
              return 0
          
          # It's a subset - check if it's the same length (same set) or smaller (proper subset)
          return 2 if len(new_coord) == len(old_coord) else 1
      
      # Check if coordinates are subsets
      x_subset = is_subset(ds.x, old_ds.x)
      y_subset = is_subset(ds.y, old_ds.y)
      
      # Check if time is subset (if time dimension exists)
      time_subset = 2  # Default to "same set" if no time dimension
      if 'time' in ds.dims and 'time' in old_ds.dims:
        time_subset = is_subset(ds.time, old_ds.time)
      
      # If new dataset is a complete subset, return it without modifying the file
      if x_subset > 0 and y_subset > 0 and time_subset > 0:
        # Only print if it's a proper subset (not same set)
        if x_subset == 1 or y_subset == 1 or time_subset == 1:
          print(f"New dataset is a subset of existing data in {NetCDF_path}. Returning subset without modification.")
        return ds
      
      # Check if coordinates match exactly
      if ds.x.equals(old_ds.x) and ds.y.equals(old_ds.y):
        # Coordinates already match, proceed with normal concat
        updated_ds = xr.concat(
          [ds, old_ds],
          dim="time"
        ).sortby("time").drop_duplicates(dim="time")
      else:
        # Coordinates differ, check if they're close enough to align
        x_diff = abs(len(ds.x) - len(old_ds.x))
        y_diff = abs(len(ds.y) - len(old_ds.y))
        # If the difference is small (e.g., 1-2 grid cells), try to align
        if x_diff <= max(len(ds.x) // 10, len(old_ds.x) // 10) + 1 and y_diff <= max(len(ds.y) // 10, len(old_ds.y) // 10) + 1:
          # Create a common set of coordinates
          # For multi-band case, we need to ensure all variables use the same coordinates
          # Choose the larger grid to avoid data loss
          target_x = old_ds.x if len(old_ds.x) >= len(ds.x) else ds.x
          target_y = old_ds.y if len(old_ds.y) >= len(ds.y) else ds.y
          # Interpolate the new dataset to match the old coordinates
          aligned_ds = ds.interp(x=target_x, y=target_y, method="nearest")
          # Now concatenate with the existing dataset
          updated_ds = xr.concat(
            [aligned_ds, old_ds],
            dim="time"
          ).sortby("time").drop_duplicates(dim="time")
        else:
          # Coordinates differ too much, provide more detailed error
          raise ValueError(
            f"Coordinate mismatch between new and existing data in {NetCDF_path}:\n"
            f"New dataset:     x={len(ds.x)}, y={len(ds.y)}\n"
            f"Existing dataset: x={len(old_ds.x)}, y={len(old_ds.y)}\n"
            f"If you used a combinable model, make sure that your polygon actually lies in the multipolygon."
            f"Old data variables: {old_ds.data_vars} <--> New data variables: {ds.data_vars}"
          )
  
  # Only save if we're not returning a subset
  if 'updated_ds' in locals():
    force_save_xr_ds(updated_ds, NetCDF_path)
    return updated_ds
  else:
    # This handles the case where we created a new file
    force_save_xr_ds(ds, NetCDF_path)
    return ds
 
def spatialy_merge_xr_datasets(ds1: xr.Dataset, ds2: xr.Dataset) -> xr.Dataset:
  combined = xr.combine_by_coords([ds1, ds2], combine_attrs="override")
  assert isinstance(combined, xr.Dataset), "Unexpected return type from combine_by_coords"  # This line is only here cause the lsp seems borked..
  return combined

def source_missing_dates(NetCDF_target: str, NetCDF_source: str, attr_name="no_data_time"):

  with xr.open_dataset(NetCDF_source, engine="h5netcdf") as ds:
    if attr_name in ds.attrs:
      source_no_data_times = ds.attrs[attr_name]
    else:
      source_no_data_times = []

  with xr.open_dataset(NetCDF_target, engine="h5netcdf") as ds:
    if attr_name in ds.attrs:
      target_no_data_times = ds.attrs[attr_name]
    else:
      target_no_data_times = []
    target_no_data_times.extend(source_no_data_times)
    target_no_data_times = list(set(target_no_data_times))

    ds.attrs[attr_name] = sorted_dates(target_no_data_times)
    force_save_xr_ds(ds, NetCDF_target)

def spatial_mean(ds: xr.Dataset) -> xr.Dataset:
    """Compute mean across spatial dimensions (x, y)."""
    return ds.mean(dim=["x", "y"], skipna=True)

def temporal_mean(ds: xr.Dataset) -> xr.Dataset:
    """Compute mean across the time dimension."""
    return ds.mean(dim="time", skipna=True)

def NetCDF_to_tiffs(
    ds: xr.Dataset | str, 
    variable_name: str, 
    polygon: shapely.Polygon, 
    save_path: str
) -> list[str]:
    """
    Extracts timesteps from NetCDF dataset, clips to polygon, and saves as GeoTIFFs.
    
    Returns
    -------
    list[str]
        List of full paths to all tiff files (existing and newly created)
    """
    # Load dataset if path is provided
    if isinstance(ds, str):
        ds = xr.open_dataset(ds)
    
    # Create save directory if it doesn't exist
    os.makedirs(save_path, exist_ok=True)
    
    # Get all timestamps
    timestamps = ds.time.values
    
    # Build list of all tiff file paths and check which are missing
    all_tiff_paths = []
    missing_indices = []
    
    for i, timestamp in enumerate(timestamps):
        timestamp_str = str(timestamp).replace(':', '-')
        output_file = os.path.join(save_path, f"{variable_name}_{timestamp_str}.tiff")
        all_tiff_paths.append(output_file)
        if not os.path.exists(output_file):
            missing_indices.append(i)
    
    # If all files exist, return early
    if not missing_indices:
        return all_tiff_paths
    
    # Get CRS from dataset
    dataset_crs = ds.crs
    
    # Reproject polygon to dataset CRS if needed
    if dataset_crs != "EPSG:4326":
        source_crs = CRS.from_epsg(4326)
        target_crs = CRS.from_string(dataset_crs)
        transformer = Transformer.from_crs(source_crs, target_crs, always_xy=True)
        
        polygon_coords = list(polygon.exterior.coords)
        transformed_coords = [transformer.transform(x, y) for x, y in polygon_coords]
        polygon = shapely.Polygon(transformed_coords)
    
    # Get coordinate arrays
    x_coords = ds.x.values
    y_coords = ds.y.values
    
    # Calculate transform for rasterio
    x_res = (x_coords[-1] - x_coords[0]) / (len(x_coords) - 1)
    y_res = (y_coords[-1] - y_coords[0]) / (len(y_coords) - 1)
    transform = from_bounds(
        x_coords[0] - x_res/2, 
        min(y_coords) - abs(y_res)/2,
        x_coords[-1] + x_res/2, 
        max(y_coords) + abs(y_res)/2,
        len(x_coords), 
        len(y_coords)
    )
    
    # Create mask for polygon
    mask = geometry_mask(
        [polygon],
        out_shape=(len(y_coords), len(x_coords)),
        transform=transform,
        invert=True
    )
    
    # Process missing timestamps
    for i in missing_indices:
        # Get data for this timestamp
        data = ds.isel(time=i)[variable_name].values
        
        # Apply mask
        masked_data = np.where(mask, data, np.nan)
        
        # Save as GeoTIFF
        timestamp_str = str(timestamps[i]).replace(':', '-')
        output_file = all_tiff_paths[i]
        
        with rasterio.open(
            output_file,
            'w',
            driver='GTiff',
            height=masked_data.shape[0],
            width=masked_data.shape[1],
            count=1,
            dtype=masked_data.dtype,
            crs=dataset_crs,
            transform=transform,
            nodata=np.nan
        ) as dst:
            dst.write(masked_data, 1)
    
    return all_tiff_paths

def interpolate_ds(ds, interpolation_method="linear"):
  daily_time = pd.date_range(
    ds.time.min().values,
    ds.time.max().values,
    freq="D"
  )
  ds_daily = ds.interp(time=daily_time, method=interpolation_method)
  return ds_daily
