#!/bin/env python

from shapely.geometry import Polygon, MultiPolygon
from io import BytesIO
from PIL import Image
import requests
import os

from app.geometry_functions import buffer_polygon
from app.image_functions import restrict_image
 
def get_esri_satellite_image(polygon: Polygon | MultiPolygon, border_size: int = 0, max_size: int = 512) -> Image.Image | None:
  """
  Fetchs and returns the satellite image of the bbox of a polygon.
  ---
  Input:
    - polygon: Shapely polygon or Multipolygon
    - max_size: How big the picture should be. This dirreclty correlates (rather equals) the resolution.
      NOTE, bigger does not increase the size (area), but rather in what resolution the picture should be downloaded.
  """

  """ TO_DO_?? max_size can be a factor or dpi and polygonsize? """

  polygon = buffer_polygon(polygon, border_size)

  min_lon, min_lat, max_lon, max_lat = polygon.bounds
  lon_diff = max_lon - min_lon
  lat_diff = max_lat - min_lat

  """ make sure we keep the aspect ratio """
  if lon_diff >= lat_diff:
    x_size = max_size
    y_size = (lat_diff * max_size) // lon_diff
  else:
    x_size = (lon_diff * max_size) // lat_diff
    y_size = max_size

  url = "".join((
    "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/export?",
    f"bbox={min_lon},{min_lat},{max_lon},{max_lat}",
    f"&bboxSR=4326&imageSR=4326",
    f"&size={x_size},{y_size}",
    "&format=png&f=image"))
  response = requests.get(url)

  if response.status_code == 200:
    return Image.open(BytesIO(response.content))

def sync_esri(polygon: Polygon, save_dir: str, border_size: int = 10, max_size: int = 512) -> Image.Image | None:
  """
  Download and save the various images required for the image analysis. These are:
    - The buffered image (default 10m, since that is the resolution of the photosynthesis)
    - The exact farm (this is also RETURNED)
  ---
  Input:
    - polygon: The shapely polygon defining the plot.
    - save_dir: Directory where the files will be saved.
    - border_size: The buffer size in meters.
    - max_size: The max pixel-size of the image. Directly corresponds to the resolution of the image.
  Output:
    - Returns the exact farm, not the buffered one.
  """
  
  os.makedirs(save_dir, exist_ok=True)
  """ This is the main file: We download this and the rest is extracted """
  if not os.path.exists(f"{save_dir}/bbox_with_border.png"):
    satellite_bbox_image = get_esri_satellite_image(polygon=polygon, border_size=border_size, max_size=max_size)
    if not satellite_bbox_image:
      raise ValueError("Unable to download satellite image. Define ways to deal with it...")
    satellite_bbox_image.save(f"{save_dir}/bbox_with_border.png", format="PNG")
  else:
    satellite_bbox_image = Image.open(f"{save_dir}/bbox_with_border.png")

  if border_size > 0:
    if not os.path.exists(f"{save_dir}/satellite_image.png"):
      satellite_image = restrict_image(satellite_bbox_image, polygon)
      satellite_image.save(f"{save_dir}/satellite_image.png", format="PNG")
    else:
      satellite_image = Image.open(f"{save_dir}/satellite_image.png")
  else:
    satellite_image = satellite_bbox_image

  return satellite_image
