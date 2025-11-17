#!/bin/env python

import xarray as xr
import time
import io
import os
import pandas as pd
from shapely.geometry   import Polygon, CAP_STYLE, JOIN_STYLE, mapping
import openeo
from datetime import datetime
from collections import defaultdict

from typing import Optional

from app.date_functions    import DateRange, group_dates_in_date_ranges, complement_in_date_range
from app.NetCDF_functions  import find_dates_in_NetCDF, restrict_NetCDF_to_polygon, sync_ds_to_NetCDF
from app.dict_functions    import sync_dicts
from app.tiff_functions    import restrict_tiff_to_polygon, sync_tiff_dc_to_NetCDF

def _discard_duplicate_dates(dates, log_file, job):
  date_map = defaultdict(list)
  for date in dates:
    date_map[date.date()].append(date)

  cleaned = []
  total_discarded = 0
  
  for date_key, dates in date_map.items():
    unique_dates = list(set(dates))

    if len(unique_dates) == 1:
      # All datetimes for this date are identical
      cleaned.append(unique_dates[0])
    else:
      # Multiple different times for the same date
      # Sort all unique times and calculate the middle point
      sorted_dates = sorted(unique_dates)
      earliest = sorted_dates[0]
      latest = sorted_dates[-1]
      middle = earliest + (latest - earliest) / 2

      # Log the discarded dates
      with open(log_file, "a") as f:
        discarded_times = [dt.strftime("%Y-%m-%d @ %H:%M:%S") for dt in sorted_dates]
        chosen_time = middle.strftime("%Y-%m-%d @ %H:%M:%S")
        f.write(f"{job} had {len(sorted_dates)} dates for {date_key}: {', '.join(discarded_times)}. Chosen {chosen_time}\n")

      cleaned.append(middle)
      total_discarded += len(sorted_dates) - 1

  # Log summary
  with open(log_file, "a") as f:
    f.write(f"{job} summary: {len(dates)} input dates -> {len(cleaned)} cleaned dates ({total_discarded} discarded)\n")
  
  return cleaned

def copernicus_batchjob__single_date_range(polygon, date_range, copernicus_model_params, retries=1,
                                           log_file="copernicus_log", sector_name=None):
  """
  For the dates in the date_range, this function returns the tiff_binaries with the corresponding dates of the given polygon.
  Input:
    - polygon = shapely polygon
    - date_range = E.g., ["2024-04-03", "2025-02-02"]
    - copernicus_model_params = dictionary; Required keys are:
      - model_name
      - if model_name == BIOPAR:
        - biopar_type
    - retries: intiger > 0. If set to negative value, the code might run forever.
      Sanity determines how often one should retry the downloads in case there is an error.
  Output:
    - A dictionary with keys:
      - tiffs = io.BytesIO of the copernicus model tiff
      - dates = the corresponding date (e.g., dc[date][2] is the date of dc[tiff][2])
  """

  
  """
    To use an other account, delete (or move)
    ~/.local/share/openeo-python-client/refresh-tokens.json
    and re-authenticate.

    Optionally, make a change of this part of the script... I.e., if moving is enough..
  """
  
  #  ---- Handle authentication ----
  eoconn = openeo.connect("openeo.dataspace.copernicus.eu").authenticate_oidc()

  print(f"Copernicus <{copernicus_model_params['model_name']}> API called for {date_range}")

  model_params = copernicus_model_params.copy()
  model_params[model_params["geometry"]] = mapping(polygon)
  model_params["temporal_extent"] = date_range
  model_params["spatial_extent"] = mapping(polygon)
  model_params["date"] = date_range

  model_name = model_params["model_name"]
  for key in model_params.keys() - {
    "date",
    "polygon",
    "aoi",
    "biopar_type",
    "namespace",
    "temporal_extent",
    "spatial_extent",
    "bands",
    "athmospheric_conditions"
  }:
    del model_params[key]

  if model_name == "SENTINEL2_L2A":
    for key in ["polygon", "date"]:
      del model_params[key]

    model = eoconn.load_collection(model_name, **model_params)
  else:
    model = eoconn.datacube_from_process(model_name, **model_params)

  try:  # Try might break paralellisation:...
    job        = model.execute_batch()
    results    = job.get_results()
    assets     = results.get_assets()
    byte_datas = [asset.load_bytes() for asset in assets]
    tiffs      = [io.BytesIO(data) for data in byte_datas]
    day_dates  = [asset.name.split("_")[-1].split(".")[0][:10] for asset in assets]
    day_dates  = [pd.to_datetime(date).date() for date in day_dates]
    dates      = []
    try:
      metadata   = results.get_metadata()
      safe_links = [l["href"] for l in metadata["links"] if l["rel"] == "derived_from"]
      date_strs  = [l.split("/")[-1].split("_")[2] for l in safe_links]
      dates      = [datetime.strptime(d, "%Y%m%dT%H%M%S") for d in date_strs]
      dates      = _discard_duplicate_dates(dates=dates, log_file=log_file, job=str(job))
    except Exception as e:
      print(e)
      with open("/data/user_data/copernicus_logs", "a") as f:
        f.write(f"job_id: {str(job)}")
        try:
          f.write(f"metadata: {str(metadata)}")  # pyright: ignore
        except:
          f.write("FAILED TO WRITE METADATA")
  except Exception as e:
    print(e)
    if retries > 0:
      print(f"{retries - 1} left.")
      retries -= 1
      time.sleep(5)
      return copernicus_batchjob__single_date_range(polygon, date_range, copernicus_model_params, retries=retries, sector_name=sector_name)
    else:
      print(f"All retries exhausted for {sector_name}, returning empty results")
    tiffs = day_dates = dates = []

  print(f"Length of tiffs: {len(tiffs)}")
  print(f"Length of dates: {len(dates)}")
  print(f"Length of day_dates: {len(day_dates)}")

  no_data_dates = complement_in_date_range(day_dates, date_range)
  if len(dates) != len(tiffs):
    if len(day_dates) != len(tiffs):
      # Neither dates nor day_dates match tiffs count
      with open(log_file, "a") as f:
        f.write(f"job_id: {job}\n")  # pyright: ignore
        if sector_name is not None:
          f.write(f"Sector name: {sector_name}\n")
        f.write("Neither dates nor day_dates agree with the size of tiffs!\n")
        f.write("___\n")
        f.write(f"dates:\n{dates}\n")
        f.write("___\n")
        if len(dates) != len(day_dates):
          f.write(f"day_dates:\n{day_dates}\n")
          f.write("___\n")
        f.write("___\n")
      # Fallback to dummy dates
      dates = list(pd.date_range(start="2017-01-01", periods=len(tiffs), freq="D"))
    else:
      # day_dates matches tiffs, but dates doesn't - filter dates to match day_dates
      dates_fixed = [date for date in dates if date.date() in day_dates]
      
      if len(dates_fixed) != len(tiffs):
        # Still doesn't match after filtering
        with open(log_file, "a") as f:
          f.write(f"job_id: {job}\n")  # #pyright: ignore
          if sector_name is not None:
            f.write(f"Sector name: {sector_name}\n")
          f.write("Dates do not agree with the size of tiffs after filtering!\n")
          f.write("___\n")
          f.write(f"original dates ({len(dates)}):\n{dates}\n")
          f.write("___\n")
          f.write(f"filtered dates ({len(dates_fixed)}):\n{dates_fixed}\n")
          f.write("___\n")
          f.write(f"day_dates ({len(day_dates)}):\n{day_dates}\n")
          f.write("___\n")
        # Use day_dates as fallback
        dates = [datetime.combine(day_date, datetime.min.time()) for day_date in day_dates]
      else:
        # Successfully filtered - use the filtered dates
        dates = dates_fixed

  return {"tiffs": tiffs, "dates": dates, "no_data_dates": no_data_dates}

def copernicus_batchjob(polygon, date_ranges, copernicus_model_params, retries, sector_name=None):
  """
  For the dates in union of the date_ranges, this function returns the tiff_binaries with the corresponding dates of the given polygon.

  Input:
    - polygon = shapely polygon
    - date_ranges = either a date_range (["2024-04-03", "2025-02-02"]) or a list of such date_ranges.
    - copernicus_model_params. WIll be passed to copernicus_batchjob__single_date_range
  Output:
    - A dictionary with keys:
      - tiffs = io.BytesIO of the copernicus model tiff
      - dates = the corresponding date (e.g., dc[date][2] is the date of dc[tiff][2])
  """
  
  if not isinstance(date_ranges[0], (list, tuple)):
    date_ranges = [date_ranges]
  copernicus_dc = dict()
  for date_range in date_ranges:
    copernicus_dc = sync_dicts(copernicus_dc, copernicus_batchjob__single_date_range(polygon, date_range, copernicus_model_params, retries=retries, sector_name=sector_name))
  return copernicus_dc

def sync_copernicus(
  polygon: Polygon,
  polygon_NetCDF: str,
  date_range: DateRange,
  copernicus_model_params: dict,
  cont_multipoly: Polygon = Polygon(),
  cont_multipoly_NetCDF: str = "",
  retries = 1,
  read_only: bool = True,
  sector_name: Optional[str] = None
) -> xr.Dataset:
  """
  NOTE: DOES NOT (YET) CREATE DIRECTORIES... TO_DO
  Updates a NetCDF file with the appropriate information from Copernicus.
  The aim is for this to work with numerous (all?) copernicus models. Currently tested for BIOPAR_fAPAR and SAVI.
  In the latter case, it EXPECTS a containing multipolygon with the containing filename.
  ---
  Input:
    - polygon = a shapely polygon
    - date_range = ["2024-04-03", "2025-02-02"] << Something like this.
    - copernicus_model_params: a dictionary. I need the keys
      - "combinable"   , a boolean determining weather I get only the 
      - "variable_name", a string or list of strings which will be written in the NetCDF file as the name(s) of the downloaded variable(s).
      - "needs_scaling", a boolean
    - cont_multipoly = If the model prefers returning bboxes, we group some polygons as a multipolygon
      and extract the information as one multipolygon. This makes it much more efficient.
      We then first try to extract the desired info (polygon) from this multipolygon, and if its not available, we add to it.
    - cont_multipoly_NetCDF = Path to where we store the model data of the cont_multipoly.
  Output:
    - dataset of the polygon
  """

  variable_names = copernicus_model_params["variable_name"]
  is_multi_band = isinstance(variable_names, list)
  
  if not read_only:
    if not copernicus_model_params["combinable"]:
      if is_multi_band:
        variable_name = variable_names[0]
      else:
        variable_name = variable_names
        # For multi-band case, check dates for the first variable (assuming all bands share same dates)
      missing_date_ranges = group_dates_in_date_ranges(
        complement_in_date_range(
          find_dates_in_NetCDF(polygon_NetCDF, variable_name=variable_name), date_range
        ), date_range
      )
      
      if len(missing_date_ranges) > 0:
        """ Needed since BIOPAR: fAPAR returned a trimmed polygon. Ergo: Ask for more and then restrict. """
        if copernicus_model_params["needs_scaling"]:
          buffered_polygon = polygon.buffer(0.001, cap_style=CAP_STYLE.flat, join_style=JOIN_STYLE.mitre)
          copernicus_poly_dc = copernicus_batchjob(buffered_polygon, missing_date_ranges, copernicus_model_params, retries=retries, sector_name=sector_name)
          copernicus_poly_dc["tiffs"] = [restrict_tiff_to_polygon(tiff, polygon) for tiff in copernicus_poly_dc["tiffs"]]
        else:
          copernicus_poly_dc = copernicus_batchjob(polygon, missing_date_ranges, copernicus_model_params, retries=retries, sector_name=sector_name)
        
        if any(copernicus_poly_dc[key] for key in copernicus_poly_dc.keys()):
          sync_tiff_dc_to_NetCDF(copernicus_poly_dc, polygon_NetCDF, variable_names)
    else:
      """ We can assume that our model is combinable. """

      """ Update the containing multipolygon NetCDF """
      if is_multi_band:
        variable_name = variable_names[0]
      else:
        variable_name = variable_names

      dates_in_NetCDF = find_dates_in_NetCDF(cont_multipoly_NetCDF, variable_name=variable_name)
      complement_of_NetCDF_dates = complement_in_date_range(dates_in_NetCDF, date_range)
      missing_date_ranges = group_dates_in_date_ranges(complement_of_NetCDF_dates, date_range)
        
      if len(missing_date_ranges) > 0:
        copernicus_multipoly_dc = copernicus_batchjob(
          cont_multipoly,
          missing_date_ranges,
          copernicus_model_params,
          retries=retries,
          sector_name=sector_name
        )
        sync_tiff_dc_to_NetCDF(
          copernicus_multipoly_dc,
          cont_multipoly_NetCDF,
          variable_names
        )

        if not os.path.exists(cont_multipoly_NetCDF):
          print(f"ERROR: Failed to create {cont_multipoly_NetCDF}")
          print(f"Copernicus data: dates={len(copernicus_multipoly_dc.get('dates', []))}, tiffs={len(copernicus_multipoly_dc.get('tiffs', []))}")
          return xr.Dataset()  # or handle appropriately

      """ Extract the tiffs and dates as a combined dc, and sync that dc with the polygon nc file. """
      if is_multi_band:
        # For multi-band, we need to process each band separately and then merge the datasets
        polygon_ds = xr.Dataset()
        for var_name in variable_names:
          band_ds = restrict_NetCDF_to_polygon(
            cont_multipoly_NetCDF,
            var_name,
            polygon,
            date_range,
            restrict_to_bbox=False
          )
          polygon_ds = xr.merge([polygon_ds, band_ds]) if polygon_ds else band_ds
      else:
        polygon_ds = restrict_NetCDF_to_polygon(
          cont_multipoly_NetCDF,
          variable_names,
          polygon,
          date_range,
          restrict_to_bbox=False
        )  # Causes problems with ultrasmall boxes.

      sync_ds_to_NetCDF(polygon_ds, polygon_NetCDF)

      """ Return the newly updated dataset """
  if os.path.exists(polygon_NetCDF):
    with xr.open_dataset(polygon_NetCDF, engine="h5netcdf") as ds:
      polygon_ds = ds.copy()  # TO_DO (is copy needed here?)
    return polygon_ds
  else:
    return xr.Dataset()
