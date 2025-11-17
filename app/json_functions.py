#!/bin/env python

import os
import numpy as np
import xarray as xr
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon
from typing import Optional

from app.NetCDF_functions import spatial_mean, temporal_mean

def mapping_polygon_to_farms_df(
  mapping_polygon,
  farm_id: Optional[str]=None,
  sector_id=None
) -> gpd.GeoDataFrame:
  """
  Intakes a GeoJSON Polygon and returns a farms_df with a shapely polygon.
  This is used as a bridge between the web interface and my functions.
  """

  user_data_path = "../user_data"

  if farm_id is None:
    farm_id="generic_farm"

  farm_path = os.path.join(user_data_path, farm_id)
  os.makedirs(farm_path, exist_ok=True)
  if sector_id is None:
    sector_num = 0
    sector_id = f"Parzelle_{sector_num}"
    sector_path = os.path.join(farm_path, sector_id)
    while os.path.isfile(sector_path):
      sector_num += 1
      sector_id = f"Parzelle_{sector_num}"
      sector_path = os.path.join(farm_path, sector_id)

  try:
    polygon = Polygon(*mapping_polygon['coordinates'])
  except:
    try:
      polygon = Polygon(mapping_polygon['coordinates'])
    except:
      raise ValueError("Please give mapping_polygon in GeoJSON")

  farms_df = gpd.GeoDataFrame({
    "Farm_ID": [farm_id],
    "Parzelle": [sector_id],
    "geometry": [polygon]
  })

  farms_df = farms_df.set_geometry("geometry").set_crs("4326")  # pyright: ignore
  return farms_df

# This is a placeholder. Needs to be solved MUCH better
def base_dir_dict_to_local_dir(
  response: dict,
  base_dir: str="/data",
  local_dir: str="files"
) -> dict:
  for key in response:
    response[key]["satellite_img_path"] = response[key]["satellite_img_path"].replace(base_dir, local_dir)
  return response

# Convert dataset to dictionary with basic functionality, such as spatial_mean, temporial_mean etc.
def ds_to_dict_of_basic_functionality(
  ds: xr.Dataset,
  desired_functions: list[str]
) -> dict:
  output = dict()

  if "full" in desired_functions:
    output["full"] = ds

  if "spatial_mean" in desired_functions:
    output["spatial_mean"] = spatial_mean(ds)

  if "temporal_mean" in desired_functions:
    output["temporal_mean"] = temporal_mean(ds)

  return output

def safe_to_json(obj):
    if isinstance(obj, pd.DataFrame):
        return obj.replace({np.nan: None}).to_dict(orient="records")
    if isinstance(obj, xr.Dataset):
        result = {}
        for k, v in obj.data_vars.items():
            da_dict = v.to_dict()
            da_dict = da_dict.get('data', da_dict)
            da_dict = safe_to_json(v.to_dict(data=True))
            result[k] = safe_to_json(v)
        return result
    if isinstance(obj, xr.DataArray):
        data_dict = obj.to_dict()
        # Strip the outer name layer if it exists
        if 'data_vars' in data_dict:
            data_dict = data_dict['data_vars'].get(obj.name, data_dict)
        # Replace NaNs in data
        if 'data' in data_dict:
            data_dict['data'] = _replace_nan(data_dict['data'])
        return data_dict
    if isinstance(obj, dict):
        return {k: safe_to_json(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [safe_to_json(i) for i in obj]
    if isinstance(obj, float) and pd.isna(obj):
        return None
    return obj

def _replace_nan(x):
    if isinstance(x, list):
        return [_replace_nan(i) for i in x]
    return None if pd.isna(x) else x

