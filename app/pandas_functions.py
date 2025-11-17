#!/bin/env python

import pandas as pd
import geopandas as gpd
from sklearn.preprocessing import OneHotEncoder
from sklearn.compose import ColumnTransformer


def numerify_objects(
  df: pd.DataFrame,
  columns_to_encode: list[str],
) -> pd.DataFrame:
  """
  Input:
    - df: pandas dataframe
    - columns_to_encode:
    - mandatory_columns:
    - empty_values:
  Return:
    Pandas Dataframe with only numeric values and vector_basis-type columns
  """

  ct = ColumnTransformer(
    transformers=[("ohe", OneHotEncoder(sparse_output=False), columns_to_encode)],
    remainder="passthrough"
  )

  # Get encoded data (as numpy)
  array_encoded    = ct.fit_transform(df)

  # Get column names
  passthrough_cols = [col for col in df.columns if col not in columns_to_encode]
  encoded_cols     = ct.named_transformers_["ohe"].get_feature_names_out(columns_to_encode)
  all_cols         = list(encoded_cols) + passthrough_cols

  # Make dataframe with encoded data and column names
  df_encoded = pd.DataFrame(array_encoded, columns=all_cols)  # pyright: ignore

  return df_encoded

def save_pandas(
  df: pd.DataFrame | gpd.GeoDataFrame,
  save_path: str,
  save_format: str
):
  if save_format == "parquet":
    df.set_geometry("geometry")
    df.to_parquet(f"{save_path}.parquet")
  if save_format == "csv":
    if not save_format.split(".")[-1] == "csv":
      save_path = f"{save_path}.csv"
    df.to_csv(save_path, index=False)
