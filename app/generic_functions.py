#!/bin/env python

import pandas as pd

def flatten_list(list_of_lists):
  """
  Input:  List of lists
  Output: List flattended by one level.
  """

  return [element for list in list_of_lists for element in list]

def roof(a: float | int) -> int:
  return int(a) if int(a) == a else int(a) + 1

def merge_dfs(df_old, df_new, how="merge"):
  if df_old is not None:
    if df_new is not None:
      if how == "merge":
        df_old = df_old.merge(df_new)
      elif how == "concat":
        df_old = pd.concat([df_old, df_new], ignore_index=True)
  else:
    df_old = df_new
  return df_old
