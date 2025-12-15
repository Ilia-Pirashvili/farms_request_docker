#!/bin/env python
import copy
import numpy as np

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

################################################################################
# LUE_max
################################################################################

bplut_table = {
  # Source https://modis-land.gsfc.nasa.gov/pdf/MOD17C61UsersGuideV10Feb2021.pdf
  "__keys": {
    "LUEmax": np.nan,
    "Tmin_min": np.nan,
    "Tmin_max": np.nan,
    "VPD_min": np.nan,
    "VPD_max": np.nan,
    "SLA": np.nan,
    "Q10": np.nan,
    "froot_leaf_ratio": np.nan,
    "livewood_leaf_ratio": np.nan,
    "leaf_mr_base": np.nan,
    "froot_mr_base": np.nan,
    "livewood_mr_base": np.nan,
  },
  "ENF": {  # Evergreen Needleleaf Forest
    "LUEmax": 0.000962,  # KgC/m2/d/MJ
    "Tmin_min": -8.00,  # C
    "Tmin_max": 8.31,  # C
    "VPD_min": 650.0,  # Pa
    "VPD_max": 4600.0,  # Pa
    "SLA": 14.1,  # LAI/KgC
    "Q10": 2.0,
    "froot_leaf_ratio": 1.2,
    "livewood_leaf_ratio": 0.182,
    "leaf_mr_base": 0.00604,
    "froot_mr_base": 0.00519,
    "livewood_mr_base": 0.00397
  },
  "EBF": {  # Evergreen Broadleaf Forest
    "LUEmax": 0.001268,
    "Tmin_min": -8.00,
    "Tmin_max": 9.09,
    "VPD_min": 800.0,
    "VPD_max": 3100.0,
    "SLA": 25.9,
    "Q10": 2.0,
    "froot_leaf_ratio": 1.1,
    "livewood_leaf_ratio": 0.162,
    "leaf_mr_base": 0.00604,
    "froot_mr_base": 0.00519,
    "livewood_mr_base": 0.00397
  },
  "DNF": {  # Deciduous Needleleaf Forest
    "LUEmax": 0.001086,
    "Tmin_min": -8.00,
    "Tmin_max": 10.44,
    "VPD_min": 650.0,
    "VPD_max": 2300.0,
    "SLA": 15.5,
    "Q10": 2.0,
    "froot_leaf_ratio": 1.7,
    "livewood_leaf_ratio": 0.165,
    "leaf_mr_base": 0.00815,
    "froot_mr_base": 0.00519,
    "livewood_mr_base": 0.00397
  },
  "DBF": {  # Deciduous Broadleaf Forest
    "LUEmax": 0.001165,
    "Tmin_min": -6.00,
    "Tmin_max": 9.94,
    "VPD_min": 650.0,
    "VPD_max": 1650.0,
    "SLA": 21.8,
    "Q10": 2.0,
    "froot_leaf_ratio": 1.1,
    "livewood_leaf_ratio": 0.203,
    "leaf_mr_base": 0.00778,
    "froot_mr_base": 0.00519,
    "livewood_mr_base": 0.00371
  },
  "MF": {  # Mixed forests
    "LUEmax": 0.001051,
    "Tmin_min": -7.00,
    "Tmin_max": 9.50,
    "VPD_min": 650.0,
    "VPD_max": 2400.0,
    "SLA": 21.5,
    "Q10": 2.0,
    "froot_leaf_ratio": 1.1,
    "livewood_leaf_ratio": 0.203,
    "leaf_mr_base": 0.00778,
    "froot_mr_base": 0.00519,
    "livewood_mr_base": 0.00371
  },
  "CShrub": {  # Closed Shrublands
    "LUEmax": 0.001281,
    "Tmin_min": -8.00,
    "Tmin_max": 8.61,
    "VPD_min": 650.0,
    "VPD_max": 4700.0,
    "SLA": 9.0,
    "Q10": 2.0,
    "froot_leaf_ratio": 1.0,
    "livewood_leaf_ratio": 0.079,
    "leaf_mr_base": 0.00869,
    "froot_mr_base": 0.00519,
    "livewood_mr_base": 0.00436
  },
  "OShrub": {  # Open Shrublands
    "LUEmax": 0.000841,
    "Tmin_min": -8.00,
    "Tmin_max": 8.80,
    "VPD_min": 650.0,
    "VPD_max": 4800.0,
    "SLA": 11.5,
    "Q10": 2.0,
    "froot_leaf_ratio": 1.3,
    "livewood_leaf_ratio": 0.040,
    "leaf_mr_base": 0.00519,
    "froot_mr_base": 0.00519,
    "livewood_mr_base": 0.00218
  },
  "WSavanna": {  # Woody Savannas
    "LUEmax": 0.001239,
    "Tmin_min": -8.00,
    "Tmin_max": 11.39,
    "VPD_min": 650.0,
    "VPD_max": 3200.0,
    "SLA": 27.4,
    "Q10": 2.0,
    "froot_leaf_ratio": 1.8,
    "livewood_leaf_ratio": 0.091,
    "leaf_mr_base": 0.00869,
    "froot_mr_base": 0.00519,
    "livewood_mr_base": 0.00312
  },
  "Savanna": {  # Savannas
    "LUEmax": 0.001206,
    "Tmin_min": -8.00,
    "Tmin_max": 11.39,
    "VPD_min": 650.0,
    "VPD_max": 3100.0,
    "SLA": 27.1,
    "Q10": 2.0,
    "froot_leaf_ratio": 1.8,
    "livewood_leaf_ratio": 0.051,
    "leaf_mr_base": 0.00869,
    "froot_mr_base": 0.00519,
    "livewood_mr_base": 0.00100
  },
  "Grass": {  # Grassland
    "LUEmax": 0.000860,
    "Tmin_min": -8.00,
    "Tmin_max": 12.02,
    "VPD_min": 650.0,
    "VPD_max": 5300.0,
    "SLA": 37.5,
    "Q10": 2.0,
    "froot_leaf_ratio": 2.6,
    "livewood_leaf_ratio": 0.000,
    "leaf_mr_base": 0.00980,
    "froot_mr_base": 0.00819,
    "livewood_mr_base": 0.00000
  },
  "Crop": {  # Croplands
    "LUEmax": 0.001044,  # IMPORTANT!
    "Tmin_min": -8.00,
    "Tmin_max": 12.02,
    "VPD_min": 650.0,
    "VPD_max": 4300.0,
    "SLA": 30.4,
    "Q10": 2.0,
    "froot_leaf_ratio": 2.0,
    "livewood_leaf_ratio": 0.000,
    "leaf_mr_base": 0.00980,
    "froot_mr_base": 0.00819,
    "livewood_mr_base": 0.00000
  }
}

crop_to_bplut_table = {
  "inactive_days_before_fallback_to_missing_data": 30,
  "missing_data": "Grass",  # The defaul / fallback entry for the dates where no data is provided.
  "Grass": "Grass",
  "Crop": "Crop",
  "weizen": "Crop",
  "kartoffeln": "Crop",
  "grassland": "Grass"
}
