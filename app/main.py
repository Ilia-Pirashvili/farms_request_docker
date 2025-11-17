from fastapi import FastAPI

import os
from typing import Optional

from app.date_functions      import DateRange
from app.json_functions      import base_dir_dict_to_local_dir, safe_to_json
from app.json_requests       import JSONRequest
from app.farm_initialisation import get

app = FastAPI()
BASE_DIR = "/data"

@app.get("/ls")
def ls(path: Optional[str] = None):
  path = os.path.join(BASE_DIR, path) if path else BASE_DIR
  if os.path.isdir(path):
    return {"files": os.listdir(path)}
  elif os.path.isfile(path):
    return "File exists but is a file, not a directory."
  else:
    return "Dude, get a working path"

@app.post("/request")
def request(req: JSONRequest):
  try:
    sane_response = get(
      mapping_polygon=req.geometry.model_dump(),
      farm_id=req.id.farm_id,
      sector_id=req.id.sector_id,
      date_range=DateRange(req.requests.date_range),
      satellite_img_size=req.requests.satellite_image.max_size,
      weather_data=req.requests.weather_data,
      copernicus_data=req.requests.copernicus_data
    )
    response = sane_response
    response = base_dir_dict_to_local_dir(response)
    if req.id.farm_id is None or req.id.sector_id is None:
      response = response[list(response.keys())[0]]
    else:
      response = response[((req.id.farm_id, req.id.sector_id))]
    return safe_to_json(response)
  except Exception as e:
    return {"error": str(e)}
