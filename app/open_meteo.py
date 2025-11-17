import pandas as pd
import openmeteo_requests
import requests_cache
from retry_requests import retry
from datetime import datetime, timedelta
import time

from app.dict_functions import find_common_dates_in_pandas_dict, sync_dicts
from app.date_functions import group_dates_in_date_ranges, complement_in_date_range
from app.HDF_functions import read_in_hdf_as_dict, write_dict_as_hdf

def get_weather(polygon, date_range, variables_list):
  def response_to_time_series(response, variables_list):
    def process_data(response_ts, variables_list):
      data = {
        "date": pd.date_range(
          start     = pd.to_datetime(response_ts.Time(),    unit="s", utc=True),
          end       = pd.to_datetime(response_ts.TimeEnd(), unit="s", utc=True),
          freq      = pd.Timedelta(seconds=response_ts.Interval()),
          inclusive = "left"
        )}
      for idx, variable in enumerate(variables_list):
        data[variable] = response_ts.Variables(idx).ValuesAsNumpy()
      return pd.DataFrame(data)

    time_series = dict()
    if len(variables_list["hourly"]) > 0:
      time_series["hourly"] = process_data(response.Hourly(), variables_list["hourly"])
    if len(variables_list["daily"]) > 0:
      time_series["daily"] = process_data(response.Daily(), variables_list["daily"])
    return time_series

  cache_session = requests_cache.CachedSession(".cache", expires_after=-1)
  retry_session = retry(cache_session, retries=5, backoff_factor=0.2)
  openmeteo     = openmeteo_requests.Client(session=retry_session)

  longitude, latitude = polygon.centroid.coords[0]
  params        = {
    "longitude": longitude,
    "latitude" : latitude,
    "hourly"   : variables_list["hourly"],
    "daily"    : variables_list["daily"]
  }

  # Openmeteo has different syntax for recent and old dates.
  # ## Archive works with dates ( > 5 days ) and forecast with dates ( < 30 days ).
  def split_date_range(date_range, threshold_days):
    start_date = pd.to_datetime(date_range[0])
    end_date   = pd.to_datetime(date_range[1])
    today      = datetime.today()

    part_1_end = min(today - timedelta(days=threshold_days + 1), end_date)
    part_2_start = max(today - timedelta(days=threshold_days), start_date)
  
    part_1 = [start_date.strftime("%Y-%m-%d"), part_1_end.strftime("%Y-%m-%d")] if part_1_end.date() >= start_date.date() else None
    part_2 = [part_2_start.strftime("%Y-%m-%d"), end_date.strftime("%Y-%m-%d")] if part_2_start.date() <= end_date.date() else None
    return part_1, part_2

  archive_dates, forecast_dates = split_date_range(date_range, threshold_days=5)
  archive_dc, forecast_dc = dict(), dict()
  if archive_dates:
    print("Call for archived weather API")
    url                  = "https://archive-api.open-meteo.com/v1/archive"
    params["start_date"] = archive_dates[0]
    params["end_date"]   = archive_dates[1]
    archive_responses    = openmeteo.weather_api(url, params=params)
    archive_dc           = response_to_time_series(archive_responses[0], variables_list)
  if forecast_dates:
    print("Call for current-forcast weather API")
    url                  = "https://api.open-meteo.com/v1/forecast"
    params["start_date"] = forecast_dates[0]
    params["end_date"]   = forecast_dates[1]
    forecast_responses   = openmeteo.weather_api(url, params=params)
    forecast_dc          = response_to_time_series(forecast_responses[0], variables_list)

  return sync_dicts(archive_dc, forecast_dc)

def sync_weather(polygon, weather_hdf, date_range, weather_variables):
  """
  Input:
    - polygon: shapely polygon
    - weather_hdf: string. Path to the weather_hdf file
    - date_range: A pair of strings of the type ["2014-02-02", "2015-09-07"]
    - weather_variables: dictionary with keys "daily" and "hourly" and lists as corresponding values.
      The entries of the list are name such as "temperature_2m", to be derived from the open_meteo api.
  Output:
    - The populated dictionariy with keys "daily" and "hourly", with cooresponding values being panda Dataframes.
  """
  weather_dc = read_in_hdf_as_dict(weather_hdf)
  missing_date_ranges = group_dates_in_date_ranges(
    complement_in_date_range(
      find_common_dates_in_pandas_dict(dc=weather_dc, key_variables=weather_variables), date_range
    ), date_range
  )
  for date_range in missing_date_ranges:
    weather_dc = sync_dicts(weather_dc, get_weather(polygon, date_range, weather_variables))
    time.sleep(1)
  weather_dc = {key: pd.DataFrame(weather_dc[key]).drop_duplicates(subset="date") for key in weather_dc.keys()}  # pd.DataFarme() is only there to quiet a false error from my lsp.
  """ Need to sort, since otherwise its unpredictable (and for things like NN's, its import once we drop the column names.) """
  weather_dc = {key: df[sorted([col for col in df.columns])]for key, df in weather_dc.items()}
  write_dict_as_hdf(weather_dc, weather_hdf)
  return weather_dc
