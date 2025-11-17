#!/bin/env python

from datetime import datetime, date
import pandas as pd
from typing import Any, Optional
from pydantic_core import core_schema
from pydantic import GetCoreSchemaHandler


class DateRange:
  """
  Initialised by a pair of strings of the type ["2014-02-02", "2015-09-07"].
  ---
  These represent the start and end dates, respectively.
  As expected, start must be <= end.
  
  If no argument is provided, the instance represents a "Nil" state,
  which evaluates to False in a boolean context.
  """

  def __init__(self, dates: Optional[list[str]] = []) -> None:
    if dates:
      try:
        if len(dates) != 2:
          raise ValueError("DateRange must contain exactly two dates (start, end).")
        self.start = dates[0]
        self.end   = dates[1]
        if datetime.strptime(self.start, "%Y-%m-%d") > datetime.strptime(self.end, "%Y-%m-%d"):
          raise ValueError(f"The end date {self.end} may not be before the start date {self.start}.")
      except ValueError:
        raise ValueError("False input. Please input a valid pair of strings in the format >> %Y-%m-%d << representing the start and end dates.")
    else:
      self.start = None
      self.end   = None

  def __repr__(self) -> str:
    return f"DateRange(start={self.start}, end={self.end})"

  def __getitem__(self, idx) -> str:
    if self.start is None or self.end is None:
      raise IndexError("Cannot index a Nil DateRange")
    return self.start if idx == 0 else self.end

  def __bool__(self) -> bool:
    return self.start is not None

  # --- Pydantic v2 support ---
  @classmethod
  def __get_pydantic_core_schema__(cls, _source_type: Any, handler: GetCoreSchemaHandler):  # pyright: ignore
      # schema: list[str] of length 2
      return core_schema.no_info_after_validator_function(
          cls._validate,
          core_schema.list_schema(core_schema.str_schema(), min_length=2, max_length=2),
          serialization=core_schema.plain_serializer_function_ser_schema(
              lambda v: [v.start, v.end] if v else [],
              return_schema=core_schema.list_schema(core_schema.str_schema())
          )
      )

  @classmethod
  def _validate(cls, v: list[str]) -> "DateRange":
      return cls(v)

def complement_in_date_range(
  dates: list[date] | list[datetime],
  date_range: DateRange
) -> list[date]:
  """
  Returns all the dates in date_range, which are not in dates.
  ---
  Input:
    - dates: List of datetime.date()'s.
    - date_range: Pair of string of the type ["2014-02-02", "2015-09-07"] or the DateRange class.
  Output:
    - List of datetime.date()'s.
  """
  dates = [pd.to_datetime(date).date() for date in dates]
  dates_set = set(dates)
  
  start = pd.to_datetime(date_range[0]).date()
  end = pd.to_datetime(date_range[1]).date()
  
  # Generate all dates in range at once
  all_dates = pd.date_range(start, end, freq="D").date  # pyright: ignore
  
  complement_dates = [d for d in all_dates if d not in dates_set]
  
  return complement_dates
 
def group_dates_in_date_ranges(
  dates: list[date],
  date_range: DateRange
) -> list[DateRange]:
  """
  Groups a list of datestime.dates()'s inside the given date_range into a list of DateRanges.
  ---
  Input:
    - dates: A list (or tuple) or datetime.date()'s.
    - date_range: Pair of string of the type ["2014-02-02", "2015-09-07"] or the DateRange class.
  Output:
    - A list of pairs of strings of the type [["2014-02-02", "2015-05-04"], ["2015-07-02", "2015-09-07"]].
      Their union should cover exactly the date_range, but not include the dates.
  """
  start_date = pd.to_datetime(date_range[0]).date()
  end_date   = pd.to_datetime(date_range[1]).date()

  dates = sorted([date for date in dates if start_date <= date <= end_date])
    
  if not dates:
    return []
  
  result = []
  range_start = dates[0]
  for i in range(1, len(dates)):
    if (dates[i] - dates[i - 1]).days > 1:
      result.append([range_start.strftime("%Y-%m-%d"), dates[i - 1].strftime("%Y-%m-%d")])
      range_start = dates[i]
  result.append([range_start.strftime("%Y-%m-%d"), dates[-1].strftime("%Y-%m-%d")])
  return result

def sorted_dates(
  dates: list[str] | list[date],
  preserve_type: bool =True
) -> list[str] | list[date]:
  """
  Returns the sorted list of inputed datetime.dates(). If preserve_type == True, it will try to preserve type
  ---
  Input:
    - dates: list of dates in various formats. Only pd.datetime(date) has to be able to deal with it.
    - preserve_type: bool. If True, will return string or datetime.date() as desired.
  Output:
    - datetime.date() if preserve_type=False.
      Else, will try to return in the same type as the input.
  """
  datetime_dates = [pd.to_datetime(date).date() for date in dates]
  datetime_dates = sorted(datetime_dates)
  if preserve_type:
    if isinstance(dates[0], str):
      return [date.strftime("%Y-%m-%d") for date in datetime_dates]
  return datetime_dates
