#!/bin/env python
import copy

################################################################################
# Copernicus
################################################################################

copernicus_base_namespace = "https://openeo.dataspace.copernicus.eu/openeo/processes/u:3e24e251-2e9a-438f-90a9-d4500e576574"
copernicus_model_parameters = {
  "bands": {
    "model_name"   : "SENTINEL2_L2A",
    "geometry"     : "polygon",
    "hidden_days"  : 0,
    "combinable"   : True,  # ?
    "needs_scaling": False,
    "variable_name": ["B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B11", "B12", "SCL", "AOT", "WVP"],  # Input band names here.
    "bands"        : ["B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B11", "B12", "SCL", "AOT", "WVP"]
    # https://custom-scripts.sentinel-hub.com/custom-scripts/sentinel-2/bands/
  },
  "BIOPAR": {
    "model_name"   : "BIOPAR",
    "geometry"     : "polygon",
    "hidden_days"  : 0, ### IF MONTHLY_CHECK, THIS NEEEEDS! to be 0!!!
    "combinable"   : False,
    "needs_scaling": True,
    "namespace"    : f"{copernicus_base_namespace}/BIOPAR"
  },
  "Variability_map": {
    "model_name"   : "variabilitymap",
    "geometry"     : "polygon",
    "hidden_days"  : 0, ### IF MONTHLY_CHECK, THIS NEEEEDS! to be 0!!!
    "combinable"   : False,
    "needs_scaling": True,
    "variable_name": "variability_map",
    "namespace"    : "https://raw.githubusercontent.com/ESA-APEx/apex_algorithms/refs/heads/main/algorithm_catalog/vito/variabilitymap/openeo_udp/variabilitymap.json"
  }
}

biopar_data = {
  "Fraction_Absorbed_Photosynthetic_Active_Radiation": "fAPAR",
  "Canopy_Chlorophyll_Content": "CCC",
  "Canopy_Water_Content": "CWC",
  "Fraction_of_Vegetation_Coverage": "fCOVER",
  "Leaf_Area_Index": "LAI"
}

for key, value in biopar_data.items():
  copernicus_model_parameters[key] = copy.deepcopy(copernicus_model_parameters["BIOPAR"])
  copernicus_model_parameters[key]["biopar_type"] = value
  copernicus_model_parameters[key]["variable_name"] = f"BIOPAR_{value}"
