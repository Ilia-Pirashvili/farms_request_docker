#!/bin/env pyt_newhon

from typing import Optional
import xarray as xr
import pyproj
import pandas as pd
from shapely.ops import transform
import copy

from app.NetCDF_functions import sync_ds_to_NetCDF

# https://custom-scripts.sentinel-hub.com/custom-scripts/sentinel/sentinel-2/ for documentation
composite_rule = {
  "SAVI":  lambda ds: (ds["B08"] - ds["B04"]) / (ds["B08"] + ds["B04"] + 0.5),
  "NDVI":  lambda ds: (ds["B08"] - ds["B04"]) / (ds["B08"] + ds["B04"]),
  "EVI":   lambda ds: 2.5*(ds["B08"] - ds["B04"]) / (ds["B08"] + 6*ds["B04"] - 7.5*ds["B02"] + 1),
  "EVI2":  lambda ds: 2.4*(ds["B08"] - ds["B04"]) / (ds["B08"] + ds["B04"] + 1),
  "GNDVI": lambda ds: (ds["B08"] - ds["B03"]) / (ds["B08"] + ds["B03"]),
  "MSI":   lambda ds: (ds["B08"] - ds["B04"]) / (ds["B08"] + ds["B04"]),
  "NDMI":  lambda ds: (ds["B08"] - ds["B11"]) / (ds["B08"] + ds["B11"]),  # = NDII
  "NBR":   lambda ds: (ds["B08"] - ds["B12"]) / (ds["B08"] + ds["B12"]),
  "PSRI":  lambda ds: (ds["B04"] - ds["B02"]) / ds["B06"],
  "ARI":   lambda ds: 1 / ds["B03"] - 1 / ds["B05"],
  "ARVI":  lambda ds: (ds["B8A"] - ds["B04"] - 0.106*(ds["B04"] - ds["B02"])) / (ds["B8A"] + ds["B04"] - 0.106*(ds["B04"] - ds["B02"])),
  "SIPI":  lambda ds: (ds["B08"] - ds["B02"]) / (ds["B08"] - ds["B04"]),
  "NDYI":  lambda ds: (ds["B03"] - ds["B02"]) / (ds["B03"] + ds["B02"])
}

import xarray as xr
import numpy as np
import pandas as pd

def AGS(ds: xr.Dataset, min_obs: int = 8) -> xr.DataArray:
   """
   Discount AGS.
   Rolling NDVI growth index from Sentinel-2 time series.
   For each day, uses ±30-day windows to compute NDVI(next) - NDVI(previous),
   only if each window has >= min_obs observations.

   Parameters
   ----------
   ds : xr.Dataset
     Dataset with dimensions (time, y, x) and bands 'B04' and 'B08'.
   min_obs : int
     Minimum number of observations required per 30-day window.

   Returns
   -------
   xr.DataArray
     Growth index as a function of time (same shape as ds minus band dimension).
   """

   """
   Still not done...
   """

   if "time" not in ds.dims:
     raise ValueError("Dataset must include a 'time' dimension.")

   # Compute NDVI for each time
   denom = ds["B08"] + ds["B04"]
   ndvi = xr.where(denom != 0, (ds["B08"] - ds["B04"]) / denom, np.nan)

   times = pd.to_datetime(ndvi.time.values)
   ndvi_vals = []

   for t in times:
     prev_mask = (times >= t - pd.Timedelta(days=30)) & (times < t)
     curr_mask = (times >= t) & (times < t + pd.Timedelta(days=30))
     next_mask = (times >= t + pd.Timedelta(days=30)) & (times < t + pd.Timedelta(days=60))

     prev = ndvi.sel(time=times[prev_mask])
     curr = ndvi.sel(time=times[curr_mask])
     nxt  = ndvi.sel(time=times[next_mask])

     if len(prev.time) >= min_obs and len(curr.time) >= min_obs and len(nxt.time) >= min_obs:
       prev_mean = prev.mean(dim="time", skipna=True)
       next_mean = nxt.mean(dim="time", skipna=True)
       growth = next_mean - prev_mean
     else:
       growth = xr.full_like(ndvi.isel(time=0), np.nan)

     growth = growth.assign_coords(time=t)
     ndvi_vals.append(growth)

   growth_index = xr.concat(ndvi_vals, dim="time")
   growth_index.name = "AGS"
   growth_index.attrs["description"] = (
     "Rolling NDVI-based growth index (next - previous 30-day means)"
   )

   return growth_index

def _make_AGS(
  band_source: str | xr.Dataset,
  save_path: Optional[str] = None
) -> xr.Dataset:
  if isinstance(band_source, str):
    with xr.open_dataset(band_source, engine="h5netcdf") as bands_data:
      ds = bands_data
  else:
    ds = band_source.copy()

  ds_new = AGS(ds)
  ds_new = ds_new.to_dataset(name="AGS")

  if ds.attrs.get("crs") is not None:
    ds_new.attrs["crs"] = ds.attrs.get("crs")

  if save_path is not None:
    sync_ds_to_NetCDF(ds=ds_new, NetCDF_path=save_path)

  return ds_new

def make_composite(
  composite_name,
  band_source: str | xr.Dataset,
  save_path: Optional[str] = None
) -> xr.Dataset:
  if isinstance(band_source, str):
    with xr.open_dataset(band_source, engine="h5netcdf") as bands_data:
      ds = bands_data
  else:
    ds = band_source.copy()

  ds_new = composite_rule[composite_name](ds)
  ds_new = ds_new.to_dataset(name=composite_name)

  if ds.attrs.get("crs") is not None:
    ds_new.attrs["crs"] = ds.attrs.get("crs")

  if save_path is not None:
    sync_ds_to_NetCDF(ds=ds_new, NetCDF_path=save_path)

  return ds_new

def make_SAVI(  # Soil Adjusted Vegetation Index
  bands_source: str | xr.Dataset,
  save_path: Optional[str] = None
) -> xr.Dataset:
  return make_composite("SAVI", band_source=bands_source, save_path=save_path)

def make_NDVI(  # Normalized difference vegetation index
  bands_source: str | xr.Dataset,
  save_path: Optional[str] = None
) -> xr.Dataset:
  return make_composite("NDVI", band_source=bands_source, save_path=save_path)

def make_EVI(  # Enhanced vegetation index
  bands_source: str | xr.Dataset,
  save_path: Optional[str] = None
) -> xr.Dataset:
  return make_composite("EVI", band_source=bands_source, save_path=save_path)

def make_EVI2(  # Enhanced vegetation index2
  bands_source: str | xr.Dataset,
  save_path: Optional[str] = None
) -> xr.Dataset:
  return make_composite("EVI2", band_source=bands_source, save_path=save_path)

def make_GNDVI(  # Green Normalized Difference Vegetation Index
  bands_source: str | xr.Dataset,
  save_path: Optional[str] = None
) -> xr.Dataset:
  return make_composite("GNDVI", band_source=bands_source, save_path=save_path)

def make_MSI(  # Moisture stress index
  bands_source: str | xr.Dataset,
  save_path: Optional[str] = None
) -> xr.Dataset:
  return make_composite("MSI", band_source=bands_source, save_path=save_path)

def make_NDMI(  # Normalised moisture index (water content of plant canopies)
  bands_source: str | xr.Dataset,
  save_path: Optional[str] = None
) -> xr.Dataset:
  return make_composite("NDMI", band_source=bands_source, save_path=save_path)

def make_NBR(  # Normalized Burn Ratio
  bands_source: str | xr.Dataset,
  save_path: Optional[str] = None
) -> xr.Dataset:
  return make_composite("NBR", band_source=bands_source, save_path=save_path)

def make_PSRI(  # Plant Senescence Reflectance Index
  bands_source: str | xr.Dataset,
  save_path: Optional[str] = None
) -> xr.Dataset:
  return make_composite("PSRI", band_source=bands_source, save_path=save_path)

def make_ARI(  # Anthocyanin Reflectance Index
               # provides valuable information about the physiological status of plants,
               # being considered an indicator of various types of plant stresses.
  bands_source: str | xr.Dataset,
  save_path: Optional[str] = None
) -> xr.Dataset:
  return make_composite("ARI", band_source=bands_source, save_path=save_path)

def make_ARVI(  # Atmospherically Resistant Vegetation Index
                # Helps assess the greenness and vigor of vegetation (crops) more reliably
                # in conditions with significant atmospheric interference (dust, haze).
                # Useful for crop health monitoring under less-than-ideal atmospheric conditions.
  bands_source: str | xr.Dataset,
  save_path: Optional[str] = None
) -> xr.Dataset:
  return make_composite("ARVI", band_source=bands_source, save_path=save_path)

def make_SIPI(  # Structure Insensitive Pigment Index
                # Early indicator of stress before greenness changes
  bands_source: str | xr.Dataset,
  save_path: Optional[str] = None
) -> xr.Dataset:
  return make_composite("SIPI", band_source=bands_source, save_path=save_path)

def make_NDYI(  # Normalized Difference Yellowness Index
                # Niche but useful: for detecting flowering stages, early senescence,
                # pigment changes (e.g., in canola or flowering trees).
  bands_source: str | xr.Dataset,
  save_path: Optional[str] = None
) -> xr.Dataset:
  return make_composite("NDYI", band_source=bands_source, save_path=save_path)

def make_AGS(  # Plant Senescence Reflectance Index
  bands_source: str | xr.Dataset,
  save_path: Optional[str] = None
) -> xr.Dataset:
  return _make_AGS(band_source=bands_source, save_path=save_path)

def restrict_copernicus_in_farms_dc_to_point(
  farms_dc: dict,
  farms_df: pd.DataFrame,
  point_col_name: str,
  double_sector_name: bool = True,
  point_crs = 4326,
  debug: bool = False
) -> dict:
  """
  Restricts the copernicus data in a farms_dc style to a given polygon.
  """
  # Create a deep copy to avoid modifying the original dictionary
  farms_dc_copy = copy.deepcopy(farms_dc)
  proj_src = pyproj.CRS(point_crs)
  for sector in farms_dc_copy.keys():
    # To keep things in a single directory, my sector names are often (SOC_f, 212312).
    if double_sector_name:
      point = farms_df.loc[farms_df["Parzelle"] == sector[1]][point_col_name].values[0]
    else:
      point = farms_df.loc[farms_df["Parzelle"] == sector][point_col_name].values[0]
    for copernicus_key, ds in farms_dc_copy[sector]["copernicus"].items():
      if ds is not None:
        proj_ds = pyproj.CRS(ds.crs)
        transformer_ds = pyproj.Transformer.from_crs(proj_src, proj_ds, always_xy=True).transform
        point_transform = transform(transformer_ds, point)
        # Find the containing grid pixel instead of nearest
        # Grid coordinates increase in x-direction, decrease in y-direction
        # Grid values are given at left-top corner
        x_coords = ds.x.values
        y_coords = ds.y.values
        
        # DEBUG: Calculate relative position
        x_min, x_max = x_coords.min(), x_coords.max()
        y_min, y_max = y_coords.min(), y_coords.max()
        
        # Calculate relative position (0-1 scale)
        x_relative = (point_transform.x - x_min) / (x_max - x_min) if x_max != x_min else 0.5
        y_relative = (point_transform.y - y_min) / (y_max - y_min) if y_max != y_min else 0.5
        
        # Find the grid pixel that contains the point
        # For x: find the largest grid_x that is <= point_x
        x_mask = x_coords <= point_transform.x
        if not x_mask.any():
          # Point is west of all grid pixels, use the westernmost
          containing_x_idx = 0
          x_status = "OUTSIDE_WEST"
        else:
          containing_x_idx = x_mask.nonzero()[0][-1]
          x_status = "INSIDE"
        
        # For y: find the largest grid_y that is >= point_y
        # (since y decreases and grid value is at top-left)
        y_mask = y_coords >= point_transform.y
        if not y_mask.any():
          # Point is south of all grid pixels, use the southernmost
          containing_y_idx = len(y_coords) - 1
          y_status = "OUTSIDE_SOUTH"
        else:
          containing_y_idx = y_mask.nonzero()[0][-1]
          y_status = "INSIDE"
        
        # Select the containing pixel while keeping coordinates as dimensions
        # Use slice to maintain x and y as dimensions with size 1
        point_ds = ds.isel(
          x=slice(containing_x_idx, containing_x_idx + 1),
          y=slice(containing_y_idx, containing_y_idx + 1))
        
        # Add mean values for each data variable across the entire grid
        for var_name in ds.data_vars:
          # Calculate mean across spatial dimensions (x, y) for each time step
          var_mean = ds[var_name].mean(dim=['x', 'y'])
          # Add the mean as a new variable with '_mean' suffix
          point_ds[f"{var_name}_mean"] = var_mean
        
#        # DEBUG: Check if result contains NaN and print debug info
#        has_nan = point_ds.isnull().any().values.item() if hasattr(point_ds.isnull().any(), 'values') else point_ds.isnull().any()
        
        if debug:
          sector_name = sector[1] if double_sector_name else sector
          print(f"Sector {sector_name}, Dataset {copernicus_key}:")
          print(f"  Point relative position: x={x_relative:.3f}, y={y_relative:.3f}")
          print(f"  Grid bounds: x=[{x_min:.2f}, {x_max:.2f}], y=[{y_min:.2f}, {y_max:.2f}]")
          print(f"  Point coords: x={point_transform.x:.2f}, y={point_transform.y:.2f}")
          print(f"  Grid indices: x={containing_x_idx}/{len(x_coords)}, y={containing_y_idx}/{len(y_coords)}")
          print(f"  Status: x={x_status}, y={y_status}")
#           print(f"  Result has NaN: {has_nan}")
#           if has_nan:
#               print(f"  *** NaN detected for sector {sector_name} ***")
          print()
        
        farms_dc_copy[sector]["copernicus"][copernicus_key] = point_ds
  return farms_dc_copy
