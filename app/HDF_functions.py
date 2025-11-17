#!/bin/env python

import os
import h5py
import pandas as pd

def read_in_hdf_as_dict(hdf_file_path):
  if not os.path.exists(hdf_file_path):
    return dict()
  with h5py.File(hdf_file_path, "r") as hdf_file:
    keys = list(hdf_file.keys())
    data_dict = {key: pd.read_hdf(hdf_file_path, key=key) for key in keys}
  return data_dict

def write_dict_as_hdf(pandas_dict, hdf_path):
  with pd.HDFStore(hdf_path, mode="w") as hdf_store:
    for key, df in pandas_dict.items():
      hdf_store.put(key, df, format="table")
