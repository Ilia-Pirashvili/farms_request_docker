# #!/bin/env python

import numpy as np
import rasterio
from rasterio.mask import mask
from rasterio.io   import MemoryFile
import pyproj
from shapely.ops import transform as shapely_transform
from shapely.geometry import mapping
import xarray as xr
from datetime import date
import io
import pandas as pd
import os

from typing import  Optional, overload

from app.NetCDF_functions import force_save_xr_ds, add_missing_dates_as_attr, spatialy_merge_xr_datasets

def restrict_tiff_to_polygon(tiffs, polygon, transformer=None):
  """
  Input:
    • tiff_file = A path to a tiff file or an io.Bytes() tiff file or a list thereofs:
      (Basically, anything that can be read in by >> rasterio.open(tiff_file) <<)
    • polygon   = shapely polygon
  Output:
    • A rasterio.io.MemoryFile or a listthereof (depends on input).
  """

  def mask_image(src, polygon):
    masked_image, masked_transform = mask(
      src,
      [mapping(polygon)],
      crop=True,
      nodata=src.nodata
    )
    profile = src.profile.copy()
    profile.update({
      "transform": masked_transform,
      "height"   : masked_image.shape[1],
      "width"    : masked_image.shape[2],
      "nodata"   : src.nodata
    })

    memfile = MemoryFile()
    with memfile.open(**profile) as mem_dst:
      mem_dst.write(masked_image)
    return memfile

  list_output = True
  if not isinstance(tiffs, list | tuple):
    tiffs = [tiffs]
    list_output = False


  outputs = []
  first_tiff_file = tiffs[0]
  with rasterio.open(first_tiff_file) as src:
    # Transform polygon to match the CRS of the raster if needed
    if src.crs.to_epsg() != 4326:
      if transformer is None:
        transformer = pyproj.Transformer.from_crs(
          4326, src.crs.to_string(), always_xy=True).transform
      else:
        pass
      polygon = shapely_transform(transformer, polygon)
    
    outputs.append(mask_image(src, polygon))

  if len(tiffs) > 1:
    for tiff_file in tiffs[1:]:
      with rasterio.open(tiff_file) as src:
        outputs.append(mask_image(src, polygon))

  if list_output:
    return outputs
  else:
    return outputs[0]
 
def tiffs_to_new_NetCDF(
  tiffs: io.BytesIO | str | list[io.BytesIO | str],
  date: date | str = "",
  NetCDF_path: str = "",
  variable_names: Optional[str | list[str] | list[list[str]]] = None,
  save: bool = False
) -> xr.Dataset:
  """
  Creates an xr.Dataset and optionally saves it as a NetCDF file from tiff files with optional timestamps.
  Support for multi-band tiff files.
  
  Note: If save == True, then NetCDF_path must be provided!
  Helper function to tiff_to_NetCDF.
  ---
  Input:
      - tiffs: io.BytesIO of a tiff or string or a list thereof.
      - date: datetime.date()
      - NetCDF_path: string. Path to a desired NetCDF file (assumed not to exist yet)
      - variable_names: 
          * string: Will be used as prefix for all bands in all tiffs
          * list[str]: One name per tiff file (for single-band tiffs)
          * list[list[str]]: Nested list with names for each band in each tiff
  Output:
      - xr.Dataset: The dataset created from the tiff files
  """
  date = pd.to_datetime(date)  # xr can not handle datetime.date()
  
  # Ensure tiffs is a list
  if not isinstance(tiffs, list):
    tiffs = [tiffs]
  
  # First pass to determine number of bands for each tiff
  bands_per_tiff = []
  for tiff in tiffs:
    with rasterio.open(tiff) as src:
      bands_per_tiff.append(src.count)
  
  # Process variable names based on input type
  if not variable_names:  # Check if its and empty list, an empty string or None, all at once.
    # Generate default names for each band in each tiff
    variable_names = []
    for i, num_bands in enumerate(bands_per_tiff):
      if num_bands == 1:
        variable_names.append([f"value_{i+1}"])  # pyright: ignore
      else:
        variable_names.append([f"value_{i+1}_band_{j+1}" for j in range(num_bands)])  # pyright: ignore
  
  elif isinstance(variable_names, str):
    # Use the string as prefix for all bands
    prefix = variable_names
    variable_names = []
    for i, num_bands in enumerate(bands_per_tiff):
      if num_bands == 1:
        variable_names.append([f"{prefix}_{i+1}"])  # pyright: ignore
      else:
        variable_names.append([f"{prefix}_{i+1}_band_{j+1}" for j in range(num_bands)])  # pyright: ignore
  
  elif isinstance(variable_names, list):
    if all(isinstance(item, str) for item in variable_names):
      # Flat list of strings - one per tiff file
      if len(variable_names) < len(tiffs):
        raise ValueError("Not enough variable names provided for all tiff files")
      
      # Convert to nested list
      old_names = variable_names.copy()
      variable_names = []
      for i, num_bands in enumerate(bands_per_tiff):
        if num_bands == 1:
          variable_names.append([old_names[i]])  # pyright: ignore
        else:
          variable_names.append([f"{old_names[i]}_band_{j+1}" for j in range(num_bands)])  # pyright: ignore
    
    elif all(isinstance(item, list) for item in variable_names):
      # Already a nested list, check if lengths match bands
      if len(variable_names) != len(tiffs):
        raise ValueError(f"Number of variable name lists ({len(variable_names)}) must match number of tiff files ({len(tiffs)})")
      
      for i, (names, num_bands) in enumerate(zip(variable_names, bands_per_tiff)):
        if len(names) != num_bands:
          raise ValueError(f"Tiff {i+1} has {num_bands} bands but {len(names)} variable names were provided")
    
    else:
      raise ValueError("variable_names must be a string, list of strings, or list of lists of strings")
  
  # Process the first tiff to get metadata and start building the dataset
  with rasterio.open(tiffs[0]) as src:
    transform = src.transform
    height = src.height
    width = src.width
    crs = src.crs
    x_coords = np.arange(width) * transform[0] + transform[2]
    y_coords = np.arange(height) * transform[4] + transform[5]
  
  # Create an empty dataset with coordinates
  ds = xr.Dataset(
    coords={
      "time": [date],
      "y": y_coords,
      "x": x_coords
    },
    attrs={
      "crs": str(crs)
    }
  )
  
  # Process all tiffs and their bands
  for tiff_idx, tiff in enumerate(tiffs):
    with rasterio.open(tiff) as src:
      # Check dimensions
      if src.height != height or src.width != width:
        raise ValueError(f"All tiffs must have the same dimensions. "
                           f"Expected {height}x{width}, got {src.height}x{src.width} for tiff {tiff_idx+1}")
      
      # TO_DO: ensure same crs, starting point and x_step, y_steps.
      
      # Process each band in the current tiff
      for band_idx in range(1, src.count + 1):
        data = src.read(band_idx)
        var_name = variable_names[tiff_idx][band_idx - 1]
        ds[var_name] = (("time", "y", "x"), data[np.newaxis, :, :])
  
  if save:
    dirname = os.path.dirname(NetCDF_path)
    if dirname:  # If we save in the same dir as the script, this will be empty.
      os.makedirs(dirname, exist_ok=True)
    ds.to_netcdf(NetCDF_path, unlimited_dims=["time"])
  
  return ds

def update_NetCDF_with_tiffs(
  tiffs: io.BytesIO | str | list[io.BytesIO | str],
  date: date | str,
  NetCDF: xr.Dataset | str,
  variable_names: Optional[str | list[str]] = None,
  save: bool = True,
  merge_mode: str = "temporal"
) -> xr.Dataset:
  """
  Updates an existing NetCDF with tiff files (and corresponding date)
  tiff_to_NetCDF should be used instead.
  ---
  Input:
    - tiffs: io.BytesIO of a tiff data or a list thereof
    - date: datetime.date()
    - NetCDF: If string, it should be the path to an existing NetCDF, else, it should be an xr.Dataset.
    - varibale_names: String or list of strings. Name of the variable to be updated.
    - save: bool, determines if saved or not. NOTE: if save == True, NetCDF must be string (path where its saved).
    - merge_mode: string in ["temporal", "spatial"]. If temporal, merge will be on "time".
      If spatial, merge will be on "x" and "y".
  Output:
    - The merged dataset.
  """
  known_merge_types = ["temporal", "spatial"]
  if merge_mode not in known_merge_types:
    raise ValueError(f"Selected merge_mode: {merge_mode} is not in {known_merge_types}")

  if not isinstance(tiffs, list):
    tiffs = [tiffs]

  new_ds = tiffs_to_new_NetCDF(
    tiffs=tiffs,
    date=date,
    variable_names=variable_names,
    save=False
  )

  if isinstance(NetCDF, xr.Dataset):
    old_ds = NetCDF
  elif isinstance(NetCDF, str):
    with xr.open_dataset(NetCDF, engine="h5netcdf") as ds:
      old_ds = ds.load()

  if merge_mode == "temporal":
    updated_ds = xr.concat(
      [old_ds, new_ds],
      dim = "time"
    ).sortby("time").drop_duplicates(dim="time")
  else:  # merge_mode == "spatial":
    updated_ds = spatialy_merge_xr_datasets(old_ds, new_ds)

  if save and isinstance(NetCDF, str):
    force_save_xr_ds(updated_ds, NetCDF)
  return updated_ds

""" Needed only for the lsp """
@overload
def tiffs_to_NetCDF(
  tiffs: io.BytesIO | str | list[io.BytesIO | str],
  date: date | str = "",
  NetCDF: xr.Dataset = xr.Dataset(),
  variable_names: Optional[str | list[str]] = None,
  merge_mode: str = "temporal",
  save: bool = False
) -> xr.Dataset: ...
@overload
def tiffs_to_NetCDF(
  tiffs: io.BytesIO | str | list[io.BytesIO | str],
  date: date | str = "",
  NetCDF: str = "",
  variable_names: Optional[str | list[str]] = None,
  merge_mode: str = "temporal",
  save: bool = False
) -> xr.Dataset: ...

def tiffs_to_NetCDF(
  tiffs: io.BytesIO | str | list[io.BytesIO | str],
  date: date | str = "",
  NetCDF = None,
  variable_names: Optional[str | list[str]] = None,
  merge_mode: str = "temporal",
  save: bool = False
) -> xr.Dataset:
  """
  Takes tiff files and, depending on the settings, either creates a new xr.Dataset or updates an existing NetCDF file with it.
  Update can be temporal (on "time") or spatial (on "x" and "y").
  Saving is optional.
  """

  if NetCDF:
    ds = update_NetCDF_with_tiffs(
      tiffs=tiffs,
      date=date,
      NetCDF=NetCDF,
      variable_names=variable_names,
      merge_mode=merge_mode, save=save
    )
  else:
    ds = tiffs_to_new_NetCDF(
      tiffs=tiffs,
      date=date,
      variable_names=variable_names,
      save=save
    )
  return ds

def sync_tiff_dc_to_NetCDF(dc: dict, NetCDF_path: str, variable_name: str | list[str]) -> xr.Dataset:
  if not dc["dates"] or not dc["tiffs"]:
    print(f"WARNING: No data to save for {variable_name} at {NetCDF_path}. Dates: {len(dc.get('dates', []))}, Tiffs: {len(dc.get('tiffs', []))}")
    return xr.Dataset()
  
  if os.path.exists(NetCDF_path):
    with xr.open_dataset(NetCDF_path, engine="h5netcdf") as initial_ds:
      ds = initial_ds.load()
  else:
    ds = xr.Dataset()
  
  # Track the original dataset size/state to detect if we actually added data
  original_time_size = len(ds.time) if 'time' in ds.dims else 0
  successful_additions = 0
  
  for _, date, tiff in zip(range(len(dc["dates"])), dc["dates"], dc["tiffs"]):
    try:
      # Attempt to add the tiff data
      ds = tiffs_to_NetCDF(tiffs=[tiff], date=date, NetCDF=ds, variable_names=[variable_name], merge_mode="temporal", save=False)  # pyright: ignore
      
      # Check if the dataset actually changed (got bigger)
      new_time_size = len(ds.time) if 'time' in ds.dims else 0
      if new_time_size > original_time_size + successful_additions:
        successful_additions += 1
      else:
        print(f"Warning: No data added for {date} (possible API issue)")
        # Optionally restore the previous state if needed
            
    except Exception as e:
      print(f"Failed to process tiff for date {date}: {e}")
      continue

  # Only save if we actually added valid data
  if successful_additions > 0:
    force_save_xr_ds(ds, NetCDF_path)
    print(f"Saved {successful_additions}/{len(dc['dates'])} valid records to {NetCDF_path}")
    
    if os.path.exists(NetCDF_path):
      add_missing_dates_as_attr(NetCDF_path, dc["no_data_dates"])
  else:
    print(f"No valid data to save for {NetCDF_path} - file unchanged")
  
  return ds
