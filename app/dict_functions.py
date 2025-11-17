#!/bin/env python
 
import pandas as pd


def sync_dicts(dict1, dict2):
  """
  Input:  2 dictionaries. Currently only handling lists, tuples and pandas dataframes.
          New ones can be added at will.
          For the same key, the dictionaries should have the same type of data!
  Output: One dictionaries with the union of keys. On the intersection, we do further sensible merges:
          Such as merging pandas dataframes and extending lists.
  NOTE:   tuples are returned as lists!
  """

  # If nothing is to do...
  if not dict2:
    return dict1
  if not dict1:
    return dict2

  combined_keys = set(dict1.keys()).union(dict2.keys())
  combined_dict = dict()

  for key in combined_keys:
    if key in dict1 and key in dict2:
      # Various merges, depending on what the values are. Can be added new ones as desired.

      ## Pandas Dataframe and Series -- Will be concatenated! Can be made more complex to also handle mergers, in the future.
      if (isinstance(dict1[key], (pd.DataFrame, pd.Series))) and (isinstance(dict2[key], (pd.DataFrame, pd.Series))):
        if "date" in dict1[key]:
          merged_value = pd.concat([dict1[key], dict2[key]], ignore_index=True).sort_values(by="date").reset_index(drop=True)
        else:
          merged_value = pd.concat([dict1[key], dict2[key]], ignore_index=True).reset_index(drop=True)
      ## Lists and tuples
      elif (isinstance(dict1[key], (list, tuple))) and (isinstance(dict2[key], (list, tuple))):
        ## In case it was a tuple. NOTE, tuples will be returned as lists...
        merged_value = list(dict1[key])
        merged_value.extend(list(dict2[key]))
      else:
        raise ValueError(f"Something is funky with dict1: {dict1} or dict2: {dict2}")

      combined_dict[key] = merged_value
    elif key in dict1:
      combined_dict[key] = dict1[key]
    elif key in dict2:
      combined_dict[key] = dict2[key]

  return combined_dict

def find_common_dates_in_pandas_dict(dc, key_variables=None, date_column_name="date"):
  """
  Finds the dates where all required columns are present across specified dataframes.
  
  Input:
  -----------
  dc : dict
      Dictionary of pandas dataframes
  key_variables : dict or None
      Dictionary mapping dataframe keys to lists of required column names.
      Example: {'daily': ['rain_sum', 'sunshine_duration'], 'hourly': ['snow_depth']}
      If None, all keys will be considered with all their columns.
  date_column_name : str
      Name of the date column in the dataframes (default: "date")
  
  Output:
  --------
  list
      Sorted list of datetime.date objects where ALL required columns have data
  
  Notes:
  ------
  - Returns dates that have complete data for all specified key_variables
  - If any required column is missing data for a date, that date is excluded
  - Useful for identifying which dates need additional data downloads
  """
  if key_variables is None:
    # If no key_variables specified, use all keys with all their columns
    examined_keys = dc.keys()
    key_variables = {key: [col for col in dc[key].columns if col != date_column_name] 
                    for key in examined_keys}
  else:
    # Only examine keys that exist in both key_variables and dc
    examined_keys = set(key_variables.keys()).intersection(set(dc.keys()))
  
  if not examined_keys:
    return []
  
  # Collect dates for each key where ALL required columns have non-null data
  date_sets = []
  
  for key in examined_keys:
    df = dc[key]
    required_columns = key_variables[key]
    
    # Check which required columns actually exist in the dataframe
    missing_columns = set(required_columns) - set(df.columns)
    if missing_columns:
      # If any required column is completely missing, no dates are valid for this key
      return []
    
    # Filter to rows where all required columns have non-null values
    mask = df[required_columns].notna().all(axis=1)
    valid_df = df[mask]
    
    # Extract dates and convert to date objects
    dates = pd.to_datetime(valid_df[date_column_name]).dt.date.unique()
    date_sets.append(set(dates))
  
  # Find intersection of all date sets (dates present in all keys with all required columns)
  if date_sets:
    common_dates = set.intersection(*date_sets)
    return sorted(common_dates)
  else:
    return []
