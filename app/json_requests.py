#!/bin/env python

from pydantic import BaseModel
from typing import Optional

class PolygonModel(BaseModel):
  type: str
  coordinates: list

class PolygonID(BaseModel):
  farm_id: Optional[str] = None
  sector_id: Optional[str] = None

class WeatherData(BaseModel):
  daily: Optional[list[str]] = None
  hourly: Optional[list[str]] = None

class CopernicusComposites(BaseModel):
  ndvi: Optional[list[str]] = None             # All will be have mean, median, full; not bool.
  savi: Optional[list[str]] = None             # ==
  evi: Optional[list[str]] = None              # ==
  evi2: Optional[list[str]] = None             # ==
  gndvi: Optional[list[str]] = None            # ==
  msi: Optional[list[str]] = None              # ==
  ndmi: Optional[list[str]] = None             # ==
  nbr: Optional[list[str]] = None              # ==
  psri: Optional[list[str]] = None             # ==
  ari: Optional[list[str]] = None              # ==
  arvi: Optional[list[str]] = None             # ==
  sipi: Optional[list[str]] = None             # ==
  ndyi: Optional[list[str]] = None             # ==
  ags: Optional[list[str]] = None              # ==

class CopernicusData(BaseModel):
  photosynthesis: Optional[list[str]] = None   # All will be have mean, median, full; not bool.
  nitrogen: Optional[list[str]] = None         # ==
  canopy_water: Optional[list[str]] = None     # ==
  vegetation_cover: Optional[list[str]] = None # ==
  leaf_area: Optional[list[str]] = None        # ==
  bands: Optional[list[str]] = None            # Shold also allow to specify which bands, and
                                               # weather atmospheric conditions should be
                                               # included.
  composites: Optional[CopernicusComposites] = None  # List of composites, such as NDVI etc.

  retries: Optional[int] = None
  read_only: Optional[bool] = None

class SatelliteImage(BaseModel):
  max_size: int = 0

class BackgroundData(BaseModel):
  tiff_data: Optional[list] = None
  save_background_data: Optional[bool] = None

class RequestOptions(BaseModel):
  cache: Optional[bool] = None

class DataRequest(BaseModel):
  date_range: Optional[list[str]] = None
  satellite_image: SatelliteImage = SatelliteImage()
  weather_data: Optional[WeatherData]  = None
  copernicus_data: Optional[CopernicusData] = None
  background_data: Optional[BackgroundData] = None

class JSONRequest(BaseModel):
  geometry: PolygonModel
  id: PolygonID
  requests: DataRequest
  options: Optional[RequestOptions] = None
