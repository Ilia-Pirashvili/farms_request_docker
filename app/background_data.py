# #!/bin/env python

from shapely.geometry import Polygon, MultiPolygon
import geopandas as gpd
import io
import xarray as xr

from typing import Optional

from app.geometry_functions import intersect_polygon_with_multipolygon_gdf
from app.tiff_functions import tiffs_to_NetCDF, restrict_tiff_to_polygon
# 
def get_soil_map_200k(
  polygon: Polygon | MultiPolygon,
  soil_map_200k_geo_layer: gpd.GeoDataFrame,
  soil_map_200k_info_layer: gpd.GeoDataFrame,
  soil_map_200k_bf_id_layer: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
  """
  Return the soil_map_200k info for a given polygon as a GeoPandas GeoDataFrame
  """
  # Date of data is Nov 2022 (I think)
  soil_map_200k = intersect_polygon_with_multipolygon_gdf(polygon, soil_map_200k_geo_layer)
  soil_map_200k = soil_map_200k.merge(soil_map_200k_info_layer, how="left")
  soil_map_200k = soil_map_200k.merge(soil_map_200k_bf_id_layer, how="left")
  soil_map_200k = soil_map_200k.drop(["BF_ID", "TKLE_NR"], axis=1)
  return gpd.GeoDataFrame(soil_map_200k)

def get_background_data_from_tiffs(
  polygon: Polygon | MultiPolygon,
  tiffs: io.BytesIO | str | list[io.BytesIO | str],
  variable_names: Optional[str | list[str]],
  transformer=None,
  date: str = ""
) -> xr.Dataset:
  """
  Intakes a list of tiffs and polygon and returns an xarray dataset
  with the information of the tiffs in that polygon.
  The xarray dataset comes with multiple data variables that can be inputed.
  Can also be given a date-stamp.
  """

  if not isinstance(tiffs, list):
    tiffs = [tiffs]

  tiffs = restrict_tiff_to_polygon(tiffs=tiffs, polygon=polygon, transformer=transformer)

  ds = tiffs_to_NetCDF(
    tiffs=tiffs,
    variable_names=variable_names,
    date=date
  )

  return ds
