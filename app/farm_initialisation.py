#!/bin/env python

import os
import time
import pyproj
import rasterio
import pandas as pd
import geopandas as gpd
from datetime import date

from typing import Optional

from app.defaults                import copernicus_model_parameters, crop_to_bplut_table, bplut_table
from app.esri_satellite_images   import sync_esri
from app.download_copernius_data import sync_copernicus
from app.open_meteo              import sync_weather
from app.date_functions          import DateRange
from app.background_data         import get_soil_map_200k, get_background_data_from_tiffs
from app.generic_functions       import merge_dfs
from app.pandas_functions        import save_pandas
from app.process_copernicus_data import make_SAVI, make_NDVI, make_EVI, make_EVI2, make_GNDVI, make_MSI, make_NDMI, make_NBR, make_PSRI, make_ARI
from app.process_copernicus_data import make_ARVI, make_SIPI, make_NDYI, make_AGS
from app.json_requests           import CopernicusData, WeatherData, CropData, Photosynthesis
from app.json_functions          import mapping_polygon_to_farms_df, ds_to_dict_of_basic_functionality
from app.NetCDF_functions        import NetCDF_to_tiffs, interpolate_ds
from app.photosynthesis          import build_crop_timeseries, align_crop_ts_dates_to_ds_times, fetch_par_polygon, multiply_ds_by_daily_df

user_data_path = "/data/user_data"
generic_data_path = "/data/background_data"


def get(
  mapping_polygon,
  farm_id: Optional[str] = None,
  sector_id: Optional[str] = None,
  date_range: Optional[DateRange] = None,
  satellite_img_size: Optional[int] = 0,
  weather_data: Optional[WeatherData] = WeatherData(),
  copernicus_data: Optional[CopernicusData] = None,
  crop_data: Optional[list[CropData]] = None,
  photosynthesis: Optional[Photosynthesis] = None
) -> dict:
  """
  Helper function for the main sync_farms function.
  ---
  Returns JSON data where possible (e.g., when returning numeric data such as photosynthesis).
  Examples where its not possible are satellite snapshots (images).
  TO_DO: allow data output specialisation.
  """

  farms_df = mapping_polygon_to_farms_df(
    mapping_polygon=mapping_polygon,
    sector_id=sector_id,
    farm_id=farm_id
  )
#   weather_frequences = []
#   if weather_data is not None:
#     if weather_data.daily is not None:
#       weather_frequences.append("daily")
#     if weather_data.hourly is not None:
#       weather_frequences.append("hourly")

  sector_data = sync_farms(
    farms=farms_df,
    date_range=date_range,
    satellite_img_size=satellite_img_size,
    weather_data=weather_data,
    copernicus_data=copernicus_data,
    crop_data=crop_data,
    background_data=[],
    copernicus_model_parameters=copernicus_model_parameters,
    photosynthesis=photosynthesis,
    save_background=False,
    copernicus_retries=1,
    read_only=False
  )

  return sector_data

""" Main sync function """
def sync_farms(
  farms,
  date_range: Optional[DateRange | list[str]] = None,
  copernicus_data: Optional[CopernicusData] = None,
  weather_data: Optional[WeatherData] = None,
  satellite_img_size: Optional[int] = 0,
  background_data: list[str] | str = [
    "soil_map_200k",
    "groundwater",
   ## tiffs
      # --- main
    "air_capacity",
    "eff_rooting_depth",
    "exch_freq_water_agr_soil",
    "cappilary_rise",
    "plant_avail_water_summer",
    "binding_str_copper",
    "binding_str_isoprot",
    "soil_depth",
    "soil_quality",
    "water_capacity",
      # --- outliers
    "geomorphy",
    "percolation_from_soil"
   ## ---
  ],
  copernicus_model_parameters = copernicus_model_parameters,
  crop_data: Optional[list[CropData]] = None,
  photosynthesis: Optional[Photosynthesis] = None,
  save_background = True,
  start_point: int = 0,
  end_point: Optional[int] = None,
  copernicus_retries: int = 0,
  read_only: bool = False,
  sleep_between: int = 0
):

  if date_range is None:
    today = date.today().strftime("%Y-%m-%d")
    date_range = ["2017-01-01", today]

# try:
  """ Start and Endpoint -- used in dl_ph.py, but also for testing """
  farms = farms.iloc[start_point:end_point]

  """ Dates """
  if isinstance(date_range, list):
    date_range = DateRange(date_range)

  """ Define global variables """
  # Photosynthesis_dc
  photosynthesis_dc = dict()
  photosynthesis_tiff_path = None

  # Copernicus
  fapar_dc = None
  fapar_tiff_path = None
  fapar = copernicus_model_parameters["Fraction_Absorbed_Photosynthetic_Active_Radiation"]

  nitrogen_dc = None
  nitrogen_tiff_path = None
  nitrogen = copernicus_model_parameters["Canopy_Chlorophyll_Content"]

  canopy_water_dc = None
  canopy_water_tiff_path = None
  canopy_water = copernicus_model_parameters["Canopy_Water_Content"]

  vegetation_cover_dc = None
  vegetation_cover_tiff_path = None
  vegetation_cover = copernicus_model_parameters["Fraction_of_Vegetation_Coverage"]

  leaf_area_dc = None
  leaf_area_tiff_path = None
  leaf_area = copernicus_model_parameters["Leaf_Area_Index"]

  bands_dc = None
#  bands_tiff_path = None
  bands = copernicus_model_parameters["bands"]

  savi_dc  = None
  savi_tiff_path = None
  ndvi_dc  = None
  ndvi_tiff_path  = None
  evi_dc   = None
  evi_tiff_path   = None
  evi2_dc  = None
  evi2_tiff_path  = None
  gndvi_dc = None
  gndvi_tiff_path = None
  msi_dc   = None
  msi_tiff_path   = None
  ndmi_dc  = None
  ndmi_tiff_path  = None
  nbr_dc   = None
  nbr_tiff_path   = None
  psri_dc  = None
  psri_tiff_path  = None
  ari_dc   = None
  ari_tiff_path   = None
  arvi_dc  = None
  arvi_tiff_path  = None
  sipi_dc  = None
  sipi_tiff_path  = None
  ndyi_dc  = None
  ndyi_tiff_path  = None
  ags_dc   = None
  ags_tiff_path   = None

  variability_xr_ds = None
  variability_map = copernicus_model_parameters["Variability_map"]

  # So that the gdf active geometry column does not have to be called "geometry".
  geometry_col_name = farms.active_geometry_name

  # Background Data
  bg_tiffs_dict_main = {
    "air_capacity"            : f"{generic_data_path}/Air_Capacity/lkwe1000_250.tif",
    "eff_rooting_depth"       : f"{generic_data_path}/Effective_rooting_depth/we1000_250.tif",
    "exch_freq_water_agr_soil": f"{generic_data_path}/Exchange_frequency_of_water_in_agricultural_soils/AHACGL1000_250.tif",
    "cappilary_rise"          : f"{generic_data_path}/Mean_annual_rate_capillary_rise/ka1000_250.tif",
    "plant_avail_water_summer": f"{generic_data_path}/Plant_available_water_in_summer_halfyear/WVPFL1000_250.tif",
    "binding_str_copper"      : f"{generic_data_path}/Relative_binding_strength__copper/FSMCu3dm1000_250.tif",
    "binding_str_isoprot"     : f"{generic_data_path}/Relative_binding_strength_of_isoproturon_in_topsoil/FOBIP3dm1000_250.tif",
    "soil_depth"              : f"{generic_data_path}/Soil_depth/PhysGru1000_250.tif",
    "soil_quality"            : f"{generic_data_path}/Soil_quality_rating_for_cropland/sqr1000_250_v10.tif",
    "water_capacity"          : f"{generic_data_path}/Water_capcaity/NFKWe1000_250.tif",
  }
  bg_tiffs_dict_outliers = {
    "geomorphy"            : f"{generic_data_path}/Geomorphographic_map/gmk1000_250.tif",
    "percolation_from_soil": f"{generic_data_path}/Mean_Annual_Rate_Percolation_from_Soil/swr1000_250.tif"
  }
  with rasterio.open(bg_tiffs_dict_main["air_capacity"]) as src:
    bg_tiff_transformer_main = pyproj.Transformer.from_crs(
      4326, src.crs.to_string(), always_xy=True).transform

  """ Set flags """
  # Background data
  if isinstance(background_data, str):
    background_data = [background_data]
  bg_flags = background_data.copy()

  if "all" in background_data:
    bg_flags = [
      "soil_map_200k",
      "groundwater",
      "tiffs"
    ]

  if "tiffs" in background_data:
    bg_flags += list(bg_tiffs_dict_main.keys())
    bg_flags += list(bg_tiffs_dict_outliers.keys())

  bg_tiffs_dict_main = {
    key: val for key, val in bg_tiffs_dict_main.items() if key in bg_flags
  }
  bg_tiffs_dict_outliers = {
    key: val for key, val in bg_tiffs_dict_outliers.items() if key in bg_flags
  }

  # Weather
  weather_variables = {"daily": [], "hourly": []}
  if weather_data is not None:
    if weather_data.daily:
      weather_variables["daily"] = weather_data.daily
    if weather_data.hourly:
      weather_variables["hourly"] = weather_data.hourly

  """ Load big datas and pass them loaded to the iterrow """
  # Soil Map 200k. -- I think, date is NOV 2022.
  if "soil_map_200k" in bg_flags:
    soil_map_200k_file        = f"{generic_data_path}/Soil_map_200k/soil_map_200k.gpkg"
    soil_map_200k_geo_layer   = gpd.GeoDataFrame(gpd.read_file(soil_map_200k_file, layer='geo_layer'))
    soil_map_200k_info_layer  = gpd.GeoDataFrame(gpd.read_file(soil_map_200k_file, layer='info_layer'))
    soil_map_200k_bf_id_layer = gpd.GeoDataFrame(gpd.read_file(soil_map_200k_file, layer='bf_id_layer'))

  output = {}
  for idx, sector in farms.iterrows():
    """ Define local variables """
    polygon  = sector[geometry_col_name]
    save_dir = f"{user_data_path}/Farm_{sector.Farm_ID}/{sector.Parzelle}"
    background_save_path = f"{save_dir}/background_data.csv"

    """ Make nescessary directories """
    os.makedirs(save_dir, exist_ok=True)
  
    """ Get satellite image """
    if satellite_img_size and satellite_img_size > 0:
      sync_esri(polygon=polygon, save_dir=save_dir, max_size=satellite_img_size)
  
    """ Get sattelite land covers """
    # I will need fapar for photosynthesis:
    # Note that spatial_mean is not needed, I need temporal or full.
    if photosynthesis is not None:
      if copernicus_data is None:
        copernicus_data = CopernicusData()
      if copernicus_data.fapar is None:
        copernicus_data.fapar = ["temporal_mean", "full"]
      elif "temporal_mean" not in copernicus_data.fapar and "full" not in copernicus_data.fapar:
        copernicus_data.fapar.extend(["temporal_mean", "full"])

    if copernicus_data is not None:
      if copernicus_data.fapar:
        fapar_xr_ds = sync_copernicus(
          polygon=polygon,
          date_range=date_range,
          copernicus_model_params=fapar,
          polygon_NetCDF=f"{save_dir}/{fapar['variable_name']}.nc",
          retries=copernicus_retries,
          read_only=read_only,
          sector_name=sector.Parzelle
        )
        fapar_dc = ds_to_dict_of_basic_functionality(ds=fapar_xr_ds, desired_functions=copernicus_data.fapar)
        fapar_tiff_path = NetCDF_to_tiffs(
          ds=fapar_xr_ds,
          variable_name=fapar['variable_name'],
          polygon=polygon,
          save_path=save_dir
        )
      if copernicus_data.nitrogen:
        nitrogen_xr_ds = sync_copernicus(
          polygon=polygon,
          date_range=date_range,
          copernicus_model_params=nitrogen,
          polygon_NetCDF=f"{save_dir}/{nitrogen['variable_name']}.nc",
          retries=copernicus_retries,
          read_only=read_only,
          sector_name=sector.Parzelle
        )
        nitrogen_dc = ds_to_dict_of_basic_functionality(ds=nitrogen_xr_ds, desired_functions=copernicus_data.nitrogen)
        nitrogen_tiff_path = NetCDF_to_tiffs(
          ds=nitrogen_xr_ds,
          variable_name=nitrogen['variable_name'],
          polygon=polygon,
          save_path=save_dir
        )
      if copernicus_data.vegetation_cover:
        vegetation_cover_xr_ds = sync_copernicus(
          polygon=polygon,
          date_range=date_range,
          copernicus_model_params=vegetation_cover,
          polygon_NetCDF=f"{save_dir}/{vegetation_cover['variable_name']}.nc",
          retries=copernicus_retries,
          read_only=read_only,
          sector_name=sector.Parzelle
        )
        vegetation_cover_dc = ds_to_dict_of_basic_functionality(ds=vegetation_cover_xr_ds, desired_functions=copernicus_data.vegetation_cover)
        vegetation_cover_tiff_path = NetCDF_to_tiffs(
          ds=vegetation_cover_xr_ds,
          variable_name=vegetation_cover['variable_name'],
          polygon=polygon,
          save_path=save_dir
        )
      if copernicus_data.leaf_area:
        leaf_area_xr_ds = sync_copernicus(
          polygon=polygon,
          date_range=date_range,
          copernicus_model_params=leaf_area,
          polygon_NetCDF=f"{save_dir}/{leaf_area['variable_name']}.nc",
          retries=copernicus_retries,
          read_only=read_only,
          sector_name=sector.Parzelle
        )
        leaf_area_dc = ds_to_dict_of_basic_functionality(ds=leaf_area_xr_ds, desired_functions=copernicus_data.leaf_area)
        leaf_area_tiff_path = NetCDF_to_tiffs(
          ds=leaf_area_xr_ds,
          variable_name=leaf_area['variable_name'],
          polygon=polygon,
          save_path=save_dir
        )
      if copernicus_data.canopy_water:
        canopy_water_xr_ds = sync_copernicus(
          polygon=polygon,
          date_range=date_range,
          copernicus_model_params=canopy_water,
          polygon_NetCDF=f"{save_dir}/{canopy_water['variable_name']}.nc",
          retries=copernicus_retries,
          read_only=read_only,
          sector_name=sector.Parzelle
        )
        canopy_water_dc = ds_to_dict_of_basic_functionality(ds=canopy_water_xr_ds, desired_functions=copernicus_data.canopy_water)
        canopy_water_tiff_path = NetCDF_to_tiffs(
          ds=canopy_water_xr_ds,
          variable_name=canopy_water['variable_name'],
          polygon=polygon,
          save_path=save_dir
        )
      if copernicus_data.bands:
        bands_xr_ds = sync_copernicus(
          polygon=polygon,
          date_range=date_range,
          copernicus_model_params=bands,
          polygon_NetCDF=f"{save_dir}/{'_'.join(bands['variable_name'])}.nc",
          cont_multipoly=polygon,                                                    # EFFECTIVELY disables efficient download.
          cont_multipoly_NetCDF=f"{save_dir}/{'_'.join(bands['variable_name'])}.nc",
          retries=copernicus_retries,
          read_only=read_only,
          sector_name=sector.Parzelle
        )
        bands_dc = ds_to_dict_of_basic_functionality(ds=bands_xr_ds, desired_functions=copernicus_data.bands)
#        bands_tiff_path = NetCDF_to_tiffs(
#          ds=bands_xr_ds,
#          variable_name=bands['variable_name'],
#          polygon=polygon,
#          save_path=save_dir
#        )
      if copernicus_data.composites:
        #  All composites are derived from bands. While not all bands are needed, just get them all since dl-time is only marginally increased like this.
        if not os.path.exists(f"{save_dir}/{'_'.join(bands['variable_name'])}.nc"):
          if copernicus_data.bands:
            bands_xr_ds = sync_copernicus(
              polygon=polygon,
              date_range=date_range,
              copernicus_model_params=bands,
              polygon_NetCDF=f"{save_dir}/{'_'.join(bands['variable_name'])}.nc",
              cont_multipoly=polygon,                                                    # EFFECTIVELY disables efficient download.
              cont_multipoly_NetCDF=f"{save_dir}/{'_'.join(bands['variable_name'])}.nc",
              retries=copernicus_retries,
              read_only=read_only,
              sector_name=sector.Parzelle
            )
          
        if copernicus_data.composites.savi:
          savi_xr_ds = make_SAVI(
            bands_source=f"{save_dir}/{'_'.join(bands['variable_name'])}.nc",
            save_path=f"{save_dir}/SAVI.nc"
          )
          savi_dc = ds_to_dict_of_basic_functionality(ds=savi_xr_ds, desired_functions=copernicus_data.composites.savi)
          savi_tiff_path = NetCDF_to_tiffs(
            ds=savi_xr_ds,
            variable_name="SAVI",
            polygon=polygon,
            save_path=save_dir
          )

        if copernicus_data.composites.ndvi:
          ndvi_xr_ds = make_NDVI(
            bands_source=f"{save_dir}/{'_'.join(bands['variable_name'])}.nc",
            save_path=f"{save_dir}/ndvi.nc"
          )
          ndvi_dc = ds_to_dict_of_basic_functionality(ds=ndvi_xr_ds, desired_functions=copernicus_data.composites.ndvi)
          ndvi_tiff_path = NetCDF_to_tiffs(
            ds=ndvi_xr_ds,
            variable_name="NDVI",
            polygon=polygon,
            save_path=save_dir
          )
        if copernicus_data.composites.evi:
          evi_xr_ds = make_EVI(
            bands_source=f"{save_dir}/{'_'.join(bands['variable_name'])}.nc",
            save_path=f"{save_dir}/evi.nc"
          )
          evi_dc = ds_to_dict_of_basic_functionality(ds=evi_xr_ds, desired_functions=copernicus_data.composites.evi)
          evi_tiff_path = NetCDF_to_tiffs(
            ds=evi_xr_ds,
            variable_name="EVI",
            polygon=polygon,
            save_path=save_dir
          )
        if copernicus_data.composites.evi2:
          evi2_xr_ds = make_EVI2(
            bands_source=f"{save_dir}/{'_'.join(bands['variable_name'])}.nc",
            save_path=f"{save_dir}/evi2.nc"
          )
          evi2_dc = ds_to_dict_of_basic_functionality(ds=evi2_xr_ds, desired_functions=copernicus_data.composites.evi2)
          evi2_tiff_path = NetCDF_to_tiffs(
            ds=evi2_xr_ds,
            variable_name="EVI2",
            polygon=polygon,
            save_path=save_dir
          )
        if copernicus_data.composites.gndvi:
          gndvi_xr_ds = make_GNDVI(
            bands_source=f"{save_dir}/{'_'.join(bands['variable_name'])}.nc",
            save_path=f"{save_dir}/gndvi.nc"
          )
          gndvi_dc = ds_to_dict_of_basic_functionality(ds=gndvi_xr_ds, desired_functions=copernicus_data.composites.gndvi)
          gndvi_tiff_path = NetCDF_to_tiffs(
            ds=gndvi_xr_ds,
            variable_name="GNDVI",
            polygon=polygon,
            save_path=save_dir
          )
        if copernicus_data.composites.msi:
          msi_xr_ds = make_MSI(
            bands_source=f"{save_dir}/{'_'.join(bands['variable_name'])}.nc",
            save_path=f"{save_dir}/msi.nc"
          )
          msi_dc = ds_to_dict_of_basic_functionality(ds=msi_xr_ds, desired_functions=copernicus_data.composites.msi)
          msi_tiff_path = NetCDF_to_tiffs(
            ds=msi_xr_ds,
            variable_name="MSI",
            polygon=polygon,
            save_path=save_dir
          )
        if copernicus_data.composites.ndmi:
          ndmi_xr_ds = make_NDMI(
            bands_source=f"{save_dir}/{'_'.join(bands['variable_name'])}.nc",
            save_path=f"{save_dir}/ndmi.nc"
          )
          ndmi_dc = ds_to_dict_of_basic_functionality(ds=ndmi_xr_ds, desired_functions=copernicus_data.composites.ndmi)
          ndmi_tiff_path = NetCDF_to_tiffs(
            ds=ndmi_xr_ds,
            variable_name="NDMI",
            polygon=polygon,
            save_path=save_dir
          )
        if copernicus_data.composites.nbr:
          nbr_xr_ds = make_NBR(
            bands_source=f"{save_dir}/{'_'.join(bands['variable_name'])}.nc",
            save_path=f"{save_dir}/nbr.nc"
          )
          nbr_dc = ds_to_dict_of_basic_functionality(ds=nbr_xr_ds, desired_functions=copernicus_data.composites.nbr)
          nbr_tiff_path = NetCDF_to_tiffs(
            ds=nbr_xr_ds,
            variable_name="NBR",
            polygon=polygon,
            save_path=save_dir
          )
        if copernicus_data.composites.psri:
          psri_xr_ds = make_PSRI(
            bands_source=f"{save_dir}/{'_'.join(bands['variable_name'])}.nc",
            save_path=f"{save_dir}/psri.nc"
          )
          psri_dc = ds_to_dict_of_basic_functionality(ds=psri_xr_ds, desired_functions=copernicus_data.composites.psri)
          psri_tiff_path = NetCDF_to_tiffs(
            ds=psri_xr_ds,
            variable_name="PSRI",
            polygon=polygon,
            save_path=save_dir
          )

        if copernicus_data.composites.ari:
          ari_xr_ds = make_ARI(
            bands_source=f"{save_dir}/{'_'.join(bands['variable_name'])}.nc",
            save_path=f"{save_dir}/ari.nc"
          )
          ari_dc = ds_to_dict_of_basic_functionality(ds=ari_xr_ds, desired_functions=copernicus_data.composites.ari)
          ari_tiff_path = NetCDF_to_tiffs(
            ds=ari_xr_ds,
            variable_name="ARI",
            polygon=polygon,
            save_path=save_dir
          )
        if copernicus_data.composites.arvi:
          arvi_xr_ds = make_ARVI(
            bands_source=f"{save_dir}/{'_'.join(bands['variable_name'])}.nc",
            save_path=f"{save_dir}/arvi.nc"
          )
          arvi_dc = ds_to_dict_of_basic_functionality(ds=arvi_xr_ds, desired_functions=copernicus_data.composites.arvi)
          arvi_tiff_path = NetCDF_to_tiffs(
            ds=arvi_xr_ds,
            variable_name="ARVI",
            polygon=polygon,
            save_path=save_dir
          )

        if copernicus_data.composites.sipi:
          sipi_xr_ds = make_SIPI(
            bands_source=f"{save_dir}/{'_'.join(bands['variable_name'])}.nc",
            save_path=f"{save_dir}/sipi.nc"
          )
          sipi_dc = ds_to_dict_of_basic_functionality(ds=sipi_xr_ds, desired_functions=copernicus_data.composites.sipi)
          sipi_tiff_path = NetCDF_to_tiffs(
            ds=sipi_xr_ds,
            variable_name="SIPI",
            polygon=polygon,
            save_path=save_dir
          )
        if copernicus_data.composites.ndyi:
          ndyi_xr_ds = make_NDYI(
            bands_source=f"{save_dir}/{'_'.join(bands['variable_name'])}.nc",
            save_path=f"{save_dir}/ndyi.nc"
          )
          ndyi_dc = ds_to_dict_of_basic_functionality(ds=ndyi_xr_ds, desired_functions=copernicus_data.composites.ndyi)
          ndyi_tiff_path = NetCDF_to_tiffs(
            ds=ndyi_xr_ds,
            variable_name="NDYI",
            polygon=polygon,
            save_path=save_dir
          )
        if copernicus_data.composites.ags:
          ags_xr_ds = make_AGS(
            bands_source=f"{save_dir}/{'_'.join(bands['variable_name'])}.nc",
            save_path=f"{save_dir}/ags.nc"
          )
          ags_dc = ds_to_dict_of_basic_functionality(ds=ags_xr_ds, desired_functions=copernicus_data.composites.ags)
          ags_tiff_path = NetCDF_to_tiffs(
            ds=ags_xr_ds,
            variable_name="AGS",
            polygon=polygon,
            save_path=save_dir
          )
##     if (isinstance(copernicus_data, list) and "variability_map" in copernicus_data) or (isinstance(copernicus_data, str) and copernicus_data == "all"):
##       variability_xr_ds = sync_copernicus(
##         polygon=polygon,
##         date_range=date_range,
##         copernicus_model_params=variability_map,
##         polygon_NetCDF=f"{save_dir}/{variability_map['variable_name']}.nc",
##         retries=copernicus_retries,
##         read_only=read_only,
##         sector_name=sector.Parzelle
##       )
##     """
##     Here, if farms are close together, they can effectively download only once.
##     This needs to be set-up, however.
##     Since in my use case this is not relevant yet, I am effectively dissabling it.
##     """

    """ Photosynthesis """
    with open("/data/log.txt", "w") as f:
      f.write(f"Start of log\n")

    if photosynthesis is not None:
      crop_names_ts = build_crop_timeseries(
        crop_data=crop_data,
        default_crop_name=crop_to_bplut_table["missing_data"],
        min_days=crop_to_bplut_table["inactive_days_before_fallback_to_missing_data"]
      )

      # Define crop_data_ts from crop_names_ts using the LEA?? data
      crop_data_ts = crop_names_ts.copy()
      for key in bplut_table["__keys"].keys():

        crop_data_ts[key] = crop_data_ts["crop_name"].apply(lambda x: bplut_table[crop_to_bplut_table[x]][key])

      # Interpolate everything to daily data, if interpolation is not None, else leave as it is
      fapar_spatial_mean_ds = None
      fapar_full_ds = None
      times = None
      if photosynthesis.interpolation is not None:
        if photosynthesis.mean is not None:
          if "spatial_mean" in photosynthesis.mean:
            fapar_spatial_mean_ds = interpolate_ds(ds=fapar_dc["spatial_mean"], interpolation_method=photosynthesis.interpolation)  # pyright: ignore
            times = fapar_spatial_mean_ds.time.values
          if "full" in photosynthesis.mean:
            fapar_full_ds = interpolate_ds(ds=fapar_dc["full"], interpolation_method=photosynthesis.interpolation)  # pyright: ignore
            times = fapar_full_ds.time.values

      else:
        if photosynthesis.mean is not None:
          if "spatial_mean" in photosynthesis.mean:
            fapar_spatial_mean_ds = fapar_dc["spatial_mean"].copy()  # pyright: ignore
            times = fapar_spatial_mean_ds.time.values
          if "full" in photosynthesis.mean:
            fapar_full_ds = fapar_dc["full"].copy()  # pyright: ignore
            times = fapar_full_ds.time.values


      # Align crop_names_ts dates to the dates in fapar_dc (note that for all keys in fapar_dc, the dates are the same).
      assert times is not None  # This should never be None! If its none, we are not sending anything.

      crop_data_ts = align_crop_ts_dates_to_ds_times(df=crop_data_ts, times=times)
      # Get PAR data (which is daily) and restrict it to the dates in fapar_dc

      with open("/data/log.txt", "a") as f:
        f.write(f"fapar_full_ds__1:\n{fapar_full_ds}\n\n")

      par_df = fetch_par_polygon(
        geom=polygon,
        start_date=pd.to_datetime(times[0]).strftime("%Y-%m-%d"),
        end_date=pd.to_datetime(times[-1]).strftime("%Y-%m-%d"),
      )
      par_df = align_crop_ts_dates_to_ds_times(df=par_df, times=times)

      # For temporial_mean and full, define photosythesis_dc[mean key thingie] as fapar_dc[mean key thingie] * par * crop_data_ts[LUA_max?]
      if photosynthesis.mean is not None:
        if "spatial_mean" in photosynthesis.mean:
          fapar_spatial_mean_ds = multiply_ds_by_daily_df(ds=fapar_spatial_mean_ds, df=par_df, column_name="PAR")
          photosynthesis_dc["spatial_mean"] = multiply_ds_by_daily_df(ds=fapar_spatial_mean_ds, df=crop_data_ts, column_name="LUEmax")
        if "full" in photosynthesis.mean:
          fapar_full_ds = multiply_ds_by_daily_df(ds=fapar_full_ds, df=par_df, column_name="PAR")
          with open("/data/log.txt", "a") as f:
            f.write(f"starting\n")
            f.write(f"fapar_full_ds__2:\n{fapar_full_ds}\n\n")
          photosynthesis_xr_ds = multiply_ds_by_daily_df(ds=fapar_full_ds, df=crop_data_ts, column_name="LUEmax")
          photosynthesis_xr_ds = photosynthesis_xr_ds.rename({"BIOPAR_fAPAR": "GPP"})
          photosynthesis_dc["full"] = photosynthesis_xr_ds
          with open("/data/log.txt", "a") as f:
            f.write(f"ending\n")

          # Now, make it to GeoTiff

          with open("/data/log.txt", "a") as f:
            f.write(f"starting\n")
            f.write(f"photosynthesis_xr_ds:\n{photosynthesis_xr_ds}\n\n")
            f.write(f"photosynthesis_dc['full']:\n{photosynthesis_dc['full']}\n\n")
          photosynthesis_tiff_path = NetCDF_to_tiffs(
            ds=photosynthesis_dc["full"],
            variable_name="GPP",
            polygon=polygon,
            save_path=save_dir
          )
          with open("/data/log.txt", "a") as f:
            f.write(f"ending\n")

    """ Get weather data """
    if weather_variables is not None and any(weather_variables.values()):
      weather_hdf = f"{save_dir}/weather.h5"
      weather_dc = sync_weather(polygon=polygon, weather_hdf=weather_hdf, date_range=date_range, weather_variables=weather_variables)
    else:
      weather_dc = None

    """ Get background data """
    if os.path.exists(background_save_path):
      bg_df = pd.read_csv(background_save_path)
    else:
      bg_df = None
      
      # Soilmap 200k
      if "soil_map_200k" in bg_flags:
        soil_map_200k = get_soil_map_200k(
          polygon=polygon,
          soil_map_200k_geo_layer=soil_map_200k_geo_layer,      # pyright: ignore
          soil_map_200k_info_layer=soil_map_200k_info_layer,    # pyright: ignore
          soil_map_200k_bf_id_layer=soil_map_200k_bf_id_layer,  # pyright: ignore
        )
        soil_map_200k["Parzelle"] = sector.Parzelle
        bg_df = merge_dfs(df_old=bg_df, df_new=soil_map_200k)

      # Tiffs
      # --- Main tiff files (which have a common raster size)
      tiff_df = None
      if len(bg_tiffs_dict_main.keys()) > 0:
        tiff_df = get_background_data_from_tiffs(
          polygon = polygon,
          tiffs = list(bg_tiffs_dict_main.values()),
          variable_names = list(bg_tiffs_dict_main.keys()),
          transformer=bg_tiff_transformer_main
        ).to_dataframe().reset_index().drop(["time", "x", "y"], axis=1)
        tiff_df["Parzelle"] = sector.Parzelle

      # --- Remaining tiff files (which have a varried raster sizes)
      new_tiff_df = None
      if len(bg_tiffs_dict_outliers.keys()) > 0:
        for key, val in bg_tiffs_dict_outliers.items():
          new_tiff_df = get_background_data_from_tiffs(
            polygon=polygon,
            tiffs=val,
            variable_names=key
          ).to_dataframe().reset_index().drop(["time", "x", "y"], axis=1)
          new_tiff_df["Parzelle"] = sector.Parzelle


      tiff_df = merge_dfs(df_old=tiff_df, df_new=new_tiff_df)

      if tiff_df is not None:
        tiff_df[tiff_df <= -9999] = 0
        tiff_df = tiff_df.mean()
        tiff_df = pd.DataFrame(tiff_df).transpose()
        bg_df = merge_dfs(df_old=bg_df, df_new=tiff_df)

      # Grundwasser
      if "groundwater" in bg_flags:
        pass

      # Optionally save background data to file for fast retrival.
      if save_background:
        if bg_df is not None:
          save_pandas(bg_df, background_save_path, "csv")

    """ Combine derived data into a single dictionary, which will then be returned """
    output[(sector.Farm_ID, sector.Parzelle)] = dict()
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"] = dict()
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"]["photosynthesis"] = photosynthesis_dc
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"]["photosynthesis"]["tiff"] = photosynthesis_tiff_path
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"][fapar["variable_name"]] = fapar_dc
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"][fapar["variable_name"]]["tiff"] = fapar_tiff_path
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"][nitrogen["variable_name"]] = nitrogen_dc
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"][nitrogen["variable_name"]]["tiff"] = nitrogen_tiff_path
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"][canopy_water["variable_name"]] = canopy_water_dc
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"][canopy_water["variable_name"]]["tiff"] = canopy_water_tiff_path
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"][vegetation_cover["variable_name"]] = vegetation_cover_dc
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"][vegetation_cover["variable_name"]]["tiff"] = vegetation_cover_tiff_path
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"][leaf_area["variable_name"]] = leaf_area_dc
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"][leaf_area["variable_name"]]["tiff"] = leaf_area_tiff_path
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"]["savi"] = savi_dc
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"]["savi"]["tiff"] = savi_tiff_path
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"]["NDVI"] = ndvi_dc
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"]["NDVI"]["tiff"] = ndvi_tiff_path
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"]["EVI"] = evi_dc   
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"]["EVI"]["tiff"] = evi_tiff_path
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"]["EVI2"] = evi2_dc 
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"]["EVI2"]["tiff"] = evi2_tiff_path
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"]["GNDVI"] = gndvi_dc
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"]["GNDVI"]["tiff"] = gndvi_tiff_path
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"]["MSI"] = msi_dc
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"]["MSI"]["tiff"] = msi_tiff_path
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"]["NDMI"] = ndmi_dc 
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"]["NDMI"]["tiff"] = ndmi_tiff_path
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"]["NBR"] = nbr_dc  
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"]["NBR"]["tiff"] = nbr_tiff_path
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"]["PSRI"] = psri_dc 
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"]["PSRI"]["tiff"] = psri_tiff_path
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"]["ARI"] = ari_dc  
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"]["ARI"]["tiff"] = ari_tiff_path
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"]["ARVI"] = arvi_dc 
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"]["ARVI"]["tiff"] = arvi_tiff_path
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"]["SIPI"] = sipi_dc 
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"]["SIPI"]["tiff"] = sipi_tiff_path
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"]["NDYI"] = ndyi_dc 
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"]["NDYI"]["tiff"] = ndyi_tiff_path
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"]["AGS"] = ags_dc  
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"]["AGS"]["tiff"] = ags_tiff_path
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"][variability_map["variable_name"]] = variability_xr_ds
    output[(sector.Farm_ID, sector.Parzelle)]["copernicus"]["_".join(bands["variable_name"])] = bands_dc
    output[(sector.Farm_ID, sector.Parzelle)]["weather"] = weather_dc
    output[(sector.Farm_ID, sector.Parzelle)]["background"] = bg_df
    output[(sector.Farm_ID, sector.Parzelle)]["satellite_img_path"] = os.path.join(save_dir, "satellite_image.png")

    if idx % 200 == -1 % 200:
      print(f"{idx + 1} sectors finished.")
    """ If run through a script, we can still get some feedback -- this will send it via notify-send """
#     if verbose_script:
#       subprocess.run(["notify-send", f"{sector.Parzelle} ended @ {datetime.now().strftime("%H-%M-%S")}"])

    for i in range(sleep_between):
      print(f"Next call in {sleep_between - i} seconds")
      time.sleep(1)
  return output

def farms_dc_to_dataframes(
  farms_dc: dict,
  double_sector_names: bool = True
) -> tuple[Optional[pd.DataFrame], Optional[pd.DataFrame], Optional[pd.DataFrame]]:
  """
  Intakes a standard farms datafrme and for each key,
    all the farm sectors are combined into one dataframe.
  """

  weather_df = None
  bg_df = None
  copernicus_df = None

  for key, farm_dc in farms_dc.items():
    sector_name = key[1] if double_sector_names else key

    # Weather
    if farm_dc["weather"] is not None:
      sector_weather_df = farm_dc["weather"]["daily"]
      sector_weather_df["Parzelle"] = sector_name
      sector_weather_df["date"] = pd.to_datetime(sector_weather_df.date).dt.tz_localize(None)
    else:
      sector_weather_df = None

    # Background data
    if farm_dc["background"] is not None:
      sector_bg_df = farm_dc["background"]
      sector_bg_df["Parzelle"] = sector_name
    else:
      sector_bg_df = None

    # Copernicus
    sector_copernicus_df = None
    for _, sector_copernicus_ds in farm_dc["copernicus"].items():
      if sector_copernicus_ds is not None:
        sector_copernicus_model_df = sector_copernicus_ds.to_dataframe().reset_index()
        if "time" in sector_copernicus_model_df.columns:
          sector_copernicus_model_df = sector_copernicus_model_df.rename(
            {"time": "date"}, axis=1)

        sector_copernicus_model_df["Parzelle"] = sector_name

        if sector_copernicus_df is None:
          sector_copernicus_df = sector_copernicus_model_df
        else:
          sector_copernicus_df = sector_copernicus_df.merge(sector_copernicus_model_df,
                                            on=["Parzelle", "date", "x", "y"])

    if weather_df is None:
      weather_df = sector_weather_df
    elif sector_weather_df is not None:
      weather_df = pd.concat([weather_df, sector_weather_df])

    if bg_df is None:
      bg_df = sector_bg_df
    elif sector_bg_df is not None:
      bg_df = pd.concat([bg_df, sector_bg_df])

    if copernicus_df is None:
      copernicus_df = sector_copernicus_df
    elif sector_copernicus_df is not None:
      copernicus_df = pd.concat([copernicus_df, sector_copernicus_df])

  return copernicus_df, weather_df, bg_df
