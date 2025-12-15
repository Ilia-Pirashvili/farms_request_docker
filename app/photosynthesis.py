import requests
import pandas as pd
import numpy as np
import xarray as xr
from collections import defaultdict
from shapely.geometry import Point, Polygon

def fetch_par_polygon(geom, start_date="2017-01-01", end_date=None):
  if end_date is None:
    end_date = start_date

  start_date = pd.to_datetime(start_date).strftime("%Y%m%d")
  end_date = pd.to_datetime(end_date).strftime("%Y%m%d")

  if isinstance(geom, Polygon):
    lon = geom.centroid.x
    lat = geom.centroid.y
  elif isinstance(geom, tuple | list):
    lon = geom[0]
    lat = geom[1]
  elif isinstance(geom, Point):
    lon = geom.x
    lat = geom.y

  url = (
    "https://power.larc.nasa.gov/api/temporal/daily/point"
    f"?parameters=ALLSKY_SFC_PAR_TOT&community=AG"
    f"&start={start_date}&end={end_date}"
    f"&latitude={lat}&longitude={lon}&format=JSON"  # pyright: ignore
  )

  r = requests.get(url, timeout=60)
  r.raise_for_status()
  data = r.json()

  values = data["properties"]["parameter"]["ALLSKY_SFC_PAR_TOT"]
  dates = sorted(values.keys())

  par_values = np.array([values[d] for d in dates], dtype="float32")
  parsed_dates = pd.to_datetime(dates, format="%Y%m%d")

  df = pd.DataFrame({"date": parsed_dates, "PAR": par_values})
  return df

def build_crop_timeseries(crop_data, default_crop_name="Grass", min_days=1):
  
  events = defaultdict(list)  # date -> list of ("plant"/"harvest", crop)
  for entry in crop_data:
    crop = entry.crop_name
    plant = pd.to_datetime(entry.plant_time)
    harvest = pd.to_datetime(entry.harvest_time)

    # integrity check: harvest must not be before plant
    if harvest < plant:
      raise ValueError(f"harvest_time before plant_time for crop {crop}")

    events[plant].append(("plant", crop))
    if (harvest - plant).days >= min_days:
      events[harvest].append(("harvest", default_crop_name))

  # decide single event per date: prefer any plant on that date over harvest
  records = []
  for date in sorted(events):
    day_events = events[date]
    plant_events = [c for t, c in day_events if t == "plant"]
    if plant_events:
      # if multiple plants same day, take the last one (preserves input order roughly)
      chosen = plant_events[-1]
    else:
      chosen = [c for t, c in day_events if t == "harvest"][-1]
    records.append((date, chosen))

  # remove consecutive duplicates (compress continuous same-crop intervals)
  compressed = []
  prev = None
  for date, crop in sorted(records):
    if crop != prev:
      compressed.append((date, crop))
      prev = crop

  df = pd.DataFrame(compressed, columns=["date", "crop_name"]).reset_index(drop=True)  # pyright: ignore
  return df

def align_crop_ts_dates_to_ds_times(df, times):
    """
    Align a static crop time-series dataframe to target times.

    For each target date in `times`, pick the last known row in `df` 
    whose date is <= target date. Values between df dates are considered constant.
    
    Parameters
    ----------
    df : pd.DataFrame
        Must have a 'date' column and other columns with values.
    times : array-like of datetime or str
        Target timestamps to align to (can include time info; only day is used).

    Returns
    -------
    pd.DataFrame
        Aligned dataframe, same length as `times`, with values filled as per last known row.
    """

    # Ensure 'date' column is datetime and df is sorted
    df = df.copy()
    df["date"] = pd.to_datetime(df["date"])
    df = df.sort_values("date")

    # Normalize target dates (drop time)
    target_dates = pd.to_datetime(times).normalize()
    target_df = pd.DataFrame({"date": target_dates})

    # Use merge_asof to pick last known df row for each target date
    out = pd.merge_asof(
        target_df,
        df,
        on="date",
        direction="backward"
    )

    return out

def multiply_ds_by_daily_df(ds, df, column_name):
    """
    Assumes ds is an xarray Dataset with coordinate 'time'.
    Multiplies all values in ds for each day by the corresponding
    daily value in df[column_name].
    """

    # Normalize dates to day resolution
    ds_days = pd.to_datetime(ds.time.values).normalize()
    df_days = pd.to_datetime(df["date"]).dt.normalize()

    # Build a pandas Series indexed by day
    daily_factor = (
        df.assign(_day=df_days)
          .set_index("_day")[column_name]
          .sort_index()
    )

    # Align factors to ds.time
    factors = xr.DataArray(
        daily_factor.loc[ds_days].values,
        coords={"time": ds.time},
        dims=("time",)
    )

    # Multiply all variables in the dataset
    out_ds = ds * factors

    # Xarray drops attributes by default (why would it not ...) We need to add it back in.
    out_ds.attrs = ds.attrs.copy()

    return out_ds
