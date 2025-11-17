#!/bin/env python
 
from shapely.geometry import Polygon, MultiPolygon, CAP_STYLE, JOIN_STYLE
from pyproj import Geod
import math
import geopandas as gpd
 
def guess_crs(polygon: Polygon | MultiPolygon) -> int:
  """
  Very superficial way of guessing the crs of a polygon. Basically, it differentiates only
  between 4326 (geodesics etc) and 25832 (rasters) by saying that small (between -180 and 180)
  values means crs=4326 and large (between 200000 and 8000000) means crs = 25832.
  ---
  Input:
    - polygon: A shapely polygon or Multipolygon.
  Output:
    - The crs in integers.
  """

  """
  POLYGON ((4345751.346853818 3264282.3684222954, 4345745.090147552 3264216.239056682, 4346049.83882011 3264202.417751182, 4346053.895011213 3264272.987420751, 4345751.346853818
            3264282.3684222954)) has crs 3035
  """
  x, y = polygon.centroid.coords[0]

  if (-180 < x < 180) and (-180 < y < 180):
    return 4326
  if (200000 < x < 800000) and (2000000 < y < 8000000):
    return 25832
  raise AssertionError("failed to dedue crs")


def metric_dist_in_crs(dist: float, crs: int, lat: float | None = None) -> tuple[float, float] | float :
  """
  Takes a distane in meters (e.g., 10 meters) and returns the appropriate number that
  represents the same distance in the given crs. Currently only works in long/lat
  (4326).
  ---
  Input:
    - dist: The distance in meters (10 for 10 meters).
    - crs: The crs we want our distance to be transformed to. E.g., 4326
  Output:
    - Appropriate representing the new distinace in inputed crs.
    - E.g., if crs=4326, we return 2 values, representing (lat, lon), as they are
      different values!
    - If, as in grids, the same value represents the same lenght in x and y, we return
      only 1 value.
  """

  if crs == 4326:
    geod = Geod(ellps="WGS84")
    R = geod.a  # Semi-major axis (approx. 6378137 meters)
    delta_lat = (dist / R) * (180 / math.pi)
    
    if not lat:
      print("". join(("NOTE: Crs=4326 (geodesic geometry) is not isometric! Distances",
                      "are not preserved and depend on the latitude (x-axis).",
                      "It was assumed to be 50, which correspondes to central Europe",
                      "It is likely unprecise, however. Try to input the real latitude.")))
      lat = 50
    
    delta_lon = (dist / (R * math.cos(math.radians(lat)))) * (180 / math.pi)
    return delta_lat, delta_lon

  raise ValueError("The inputed crs is not (yet) recognised.")

def buffer_polygon(polygon: Polygon | MultiPolygon, buffer_dist: float = 10,
                   preserve_shape: bool = True, crs: int | None = None ) -> Polygon | MultiPolygon:

  if not crs:
    crs = guess_crs(polygon)
  lat = polygon.centroid.coords[0][0] if crs == 4326 else None
  transformed_buffer_dist = metric_dist_in_crs(dist=buffer_dist, crs=crs, lat=lat)
  if isinstance(transformed_buffer_dist, tuple):
    transformed_buffer_dist = max(transformed_buffer_dist)

  if preserve_shape:
    buffered_polygon = polygon.buffer(transformed_buffer_dist, cap_style=CAP_STYLE.flat, join_style=JOIN_STYLE.mitre)
  else:
    buffered_polygon = polygon.buffer(transformed_buffer_dist)
  return buffered_polygon

def intersect_polygon_with_multipolygon_gdf(
  polygon: Polygon | MultiPolygon,
  multipolygon_gdf: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
  polygon_partition_w_multipolygon_gdf_data = []
  intersecting_polygons = multipolygon_gdf[multipolygon_gdf.intersects(polygon)]

  for _, row in intersecting_polygons.iterrows():
    intersection = row.geometry.intersection(polygon)
    if not intersection.is_empty:
      intersected_row = row.copy()
      intersected_row["geometry"] = intersection
      polygon_partition_w_multipolygon_gdf_data.append(intersected_row)

  return gpd.GeoDataFrame(polygon_partition_w_multipolygon_gdf_data, crs=multipolygon_gdf.crs)
