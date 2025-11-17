# Request machen:

Ein Aufruf wird duch ein JSON request gemacht. Die Prameter hierzu sind:
    geometry, id, requests

## geometry:
"geometry" definiert einen Schlag geometrisch. Es hat die Parameter
- type,
    Hier soll "Polygon" stehen. Weiteres wurde noch nicht getestet, aber "Multipolygon"
    könnte später für bessere effizienz eingebaut werden.
- coordinates
    Eine List von Listen von geodeischen Koordinates. Hier wird crs 4326 verwendet
    (standard long/lat Koordinaten).
    Im unteren Beispiel wird dies genauer beschrieben.

## id:
Dies definiert den Sektor. Es besteht aus zwei Informationen (zwei Parameter):
    farm_id, sector_id
Beide sind beliebige strings. Zusammen jedoch müssen diese eindeutig sein.
Diese tragen dazu bei, dass Verzeichnisse erstellt werden, mit der Struktur:
    farm_id
    └── sector_id
        ├── ...
        └── ...
wobei im sector_id all die Informationen gespreichert werden, die man runtergeladen hat.
Die speicherung dient dazu, damit man nicht erneut Daten runterladen bzw. generieren
braucht, wenn man es schon hat.
Dies sollte wie "cache" behandelt werden: Bei Bedarf / Errors können diese gelöscht
werden, jedoch sind Download-zeiten, vor allem bei den Satellitendaten enorm hoch.

## requests:

Der Hauptteil des Aufrufs erfolgt hier. Die Parameter sind wie folgt:
    date_range, satellite_image, weathre_data, copernicus_data

### date_range:
Alles außer satellite_image ist Zeitsensitiv. date_range beeinflust eben diesen Faktor.
Es bestimmt, in welchen Zeitabständen man Daten herunterladen kann.
Es wird in der form [start_date, end_date] gegeben, wobei
    start_date, end_date
beide im Format "yyyy-mm-dd" gegeben werden (als string).

- Die Einschränkung hier ist, dass start_date nach 2017-01-01 sein muss! da Copernicus
  keine Satellitendaten zuvor hat.
  *Wenn man den Datum vorher ansetzt, wird es **nicht** zu einem Error kommen, da
  z. B., Wetter noch geladen werden kann usw. Man soll aber keine Satellitendaten
  erwarten.*

### satellite_image:

    satellite_image: {
      max_size: 500
    },

### weather_data:
Hier werden Wetter und Klima-daten heruntergeladen. Man kann dies sowohl als Stündliche,
wie auch Tägliche Daten runterladen. Die Parameter sind daher
    daily, hourly
In beiden Fällen, gibt man als eine List an, welche Informationen man will. Desdo mehr
man angibt, desdo mehr "kostet" es an API calls. Die Quelle hier ist *open-meteo.com*

*Zu open-meteo.com API's. Hier sei gesagt, dass "Free / Open-Access", unter dem ich es
bis jetzt benutze, **nicht** für kommerzielle Zwecke zugelassen ist. Auf deren Webseite
steht (Stand 14.11.2025):
    **During the evaluation, prototyping, and development stages, we encourage you to
    utilise the free tier. However, for your production system, we will provide you with
    an API key to ensure that your application reliably receives weather data.**
Sprich, für die erste Phase unter der ich es benutzt habe, war dies noch erlaubt. Weiter
muss man dies aber mit dehnen abklärren. Es gibt natürlich API-keys die man sich für
monatliche Gebüren kaufen/mieten kann.*

Die Parameter sind wie folgt:
- daily:
    - weather_code
    - temperature_2m_max
    - temperature_2m_min
    - apparent_temperature_max
    - apparent_temperature_min
    - sunrise
    - sunset
    - daylight_duration
    - sunshine_duration
    - uv_index_max
    - uv_index_clear_sky_max
    - rain_sum
    - showers_sum
    - snowfall_sum
    - precipitation_sum
    - precipitation_hours
    - precipitation_probability_max
    - wind_speed_10m_max
    - wind_gusts_10m_max
    - wind_direction_10m_dominant
    - shortwave_radiation_sum
    - et0_fao_evapotranspiration
    - vapour_pressure_deficit_max
    - wet_bulb_temperature_2m_min
    - wet_bulb_temperature_2m_max
    - wet_bulb_temperature_2m_mean
    - surface_pressure_mean
    - surface_pressure_max
    - leaf_wetness_probability_mean
    - precipitation_probability_mean
    - cape_min
    - cloud_cover_mean
    - relative_humidity_2m_max
    - relative_humidity_2m_min
    - wind_gusts_10m_mean
    - wind_speed_10m_mean
    - pressure_msl_max
    - pressure_msl_min
    - dew_point_2m_min
    - dew_point_2m_max
    - snowfall_water_equivalent_sum
    - winddirection_10m_dominant
    - visibility_max
    - visibility_min
    - updraft_max
    - surface_pressure_min
    - et0_fao_evapotranspiration_sum
    - growing_degree_days_base_0_limit_50
    - temperature_2m_mean
    - apparent_temperature_mean
    - cape_max
    - cape_mean
    - cloud_cover_max
    - cloud_cover_min
    - dew_point_2m_mean
    - pressure_msl_mean
    - relative_humidity_2m_mean
    - precipitation_probability_min
    - visibility_mean
    - wind_gusts_10m_min
    - wind_speed_10m_min
- hourly:
    - temperature_2m
    - relative_humidity_2m
    - dew_point_2m
    - apparent_temperature
    - precipitation_probability
    - rain
    - precipitation
    - showers
    - snowfall
    - snow_depth
    - weather_code
    - pressure_msl
    - cloud_cover
    - surface_pressure
    - cloud_cover_low
    - cloud_cover_mid
    - cloud_cover_high
    - visibility
    - evapotranspiration
    - et0_fao_evapotranspiration
    - vapour_pressure_deficit
    - temperature_180m
    - temperature_120m
    - temperature_80m
    - wind_direction_180m
    - wind_gusts_10m
    - wind_direction_120m
    - wind_direction_80m
    - wind_direction_10m
    - wind_speed_180m
    - wind_speed_120m
    - wind_speed_80m
    - wind_speed_10m
    - soil_temperature_0cm
    - soil_temperature_18cm
    - soil_temperature_6cm
    - soil_temperature_54cm
    - soil_moisture_0_to_1cm
    - soil_moisture_1_to_3cm
    - soil_moisture_3_to_9cm
    - soil_moisture_9_to_27cm
    - soil_moisture_27_to_81cm
    - uv_index
    - uv_index_clear_sky
    - is_day
    - sunshine_duration
    - wet_bulb_temperature_2m
    - total_column_integrated_water_vapour
    - cape
    - lifted_index
    - convective_inhibition
    - freezing_level_height
    - boundary_layer_height
    - shortwave_radiation
    - diffuse_radiation
    - direct_radiation
    - direct_normal_irradiance
    - global_tilted_irradiance
    - terrestrial_radiation
    - direct_radiation_instant
    - shortwave_radiation_instant
    - diffuse_radiation_instant
    - direct_normal_irradiance_instant
    - global_tilted_irradiance_instant
    - terrestrial_radiation_instant

## copernicus_data
This used the copernicus provided satellite data (Sentinel 2 mainly) and models trained
and derived by these satellite datas (the trained data is provided by the copernicus
marketplace).

*This required an API-key! The api-must be optained from copernicus. The legal aspects
of use, that is to say, in what capacity these may be used for commercial purposes, and
how to obtain more tokens (the currency of the api-calls) should be checked at the
copernicus webpage. **I take no warrenty for this!***

It has the following parameters:
    photosynthesis, nitrogen, canopy_water, vegetation_cover, leaf_area, bands, *composites*
Here, *composites* is different and will be explained below.

- For all of the others, the prameters
    spatial_mean, temporal_mean, all
  are working thus far.
    - spatial_mean will return a time-series, with each value being the mean of the sector at that time.
    - temporal_mean will return a single 2D raster-type dataset. Each value is the mean over time at that point.
    - all will return the full dataset -- a time-series of rasters / "images".

The given datasets represent the following:
- photosynthesis:
- nitrogen:
- canopy_water:
- vegetation_cover:
- leaf_area:
- bands:
    The raw data corresponding to the 10 bands
        **B02**, **B03**, **B04**, **B05**, **B06**, **B07**, **B08**, **B8A**, **B11**, **B12**
    which are values in given wavelengths.

- composites:
    **REQUIRES: bands**
    This denotes certain data, such as *normalised differential vegetaion index (ndvi)*.
    It is obtained by aggegating the raw band values.
    "composite" itself takes as parameters **not** spatial_mean etc, but itself names of
    data values, such as *ndvi*, which itself then takes values such as spatial_mean
    etc.

    The parameters are:
        ndvi, savi, evi, evi2, gndvi, msi, ndmi, nbr, psri, ari, arvi, sipi, ndyi, ags

    - These themselves now take parameters like the above ones, being
        spatial_mean, temporal_mean, all



# model call:
{
  geometry: {
    type: Polygon,
    coordinates: [
      [
        [11.354947914480557, 51.78796822768034],
        [11.355052085519442, 51.78796822768034],
        [11.355052085519442, 51.79015177231967],
        [11.352947914480557, 51.79015177231967],
        [11.352947914480557, 51.78796822768034]
      ]
    ]
  },
  id: {
    farm_id: test_farm_id,
    sector_id: test_sector_id
  },
  requests: {
    date_range: [2019-07-01, 2019-07-15],
    satellite_image: {
      max_size: 500
    },
    weather_data: {
      daily: [rain_sum, sunshine_duration],
      hourly: [snow_depth]
    },
    copernicus_data: {
      photosynthesis: [spatial_mean, temporal_mean],
      nitrogen: [spatial_mean, temporal_mean],
      canopy_water: [spatial_mean, temporal_mean],
      vegetation_cover: [spatial_mean, temporal_mean],
      leaf_area: [spatial_mean, temporal_mean],
      bands: [spatial_mean],
      composites: {
        ndvi:  [spatial_mean],
        savi:  [spatial_mean],
        evi:   [spatial_mean],
        evi2:  [spatial_mean],
        gndvi: [spatial_mean],
        msi:   [spatial_mean],
        ndmi:  [spatial_mean],
        nbr:   [spatial_mean],
        psri:  [spatial_mean],
        ari:   [spatial_mean],
        arvi:  [spatial_mean],
        sipi:  [spatial_mean],
        ndyi:  [spatial_mean],
        ags:   [spatial_mean]
      }
    }
  }
}


# model call:
{
  geometry: {
    type: Polygon,
    coordinates: [
      [
        [11.354947914480557, 51.78796822768034],
        [11.355052085519442, 51.78796822768034],
        [11.355052085519442, 51.79015177231967],
        [11.352947914480557, 51.79015177231967],
        [11.352947914480557, 51.78796822768034]
      ]
    ]
  },
  id: {
    farm_id: test_farm_id,
    sector_id: test_sector_id
  },
  requests: {
    date_range: [2019-07-01, 2019-07-15],
    satellite_image: {
      max_size: 500
    },
    weather_data: {
      daily: [rain_sum, sunshine_duration],
      hourly: [snow_depth]
    },
    copernicus_data: {
      photosynthesis: [spatial_mean, temporal_mean],
      nitrogen: [spatial_mean, temporal_mean],
      canopy_water: [spatial_mean, temporal_mean],
      vegetation_cover: [spatial_mean, temporal_mean],
      leaf_area: [spatial_mean, temporal_mean],
      bands: [spatial_mean],
      composites: {
        ndvi:  [spatial_mean],
        savi:  [spatial_mean],
        evi:   [spatial_mean],
        evi2:  [spatial_mean],
        gndvi: [spatial_mean],
        msi:   [spatial_mean],
        ndmi:  [spatial_mean],
        nbr:   [spatial_mean],
        psri:  [spatial_mean],
        ari:   [spatial_mean],
        arvi:  [spatial_mean],
        sipi:  [spatial_mean],
        ndyi:  [spatial_mean],
        ags:   [spatial_mean]
      }
    }
  }
}


# model call:
{
  geometry: {
    type: Polygon,
    coordinates: [
      [
        [11.354947914480557, 51.78796822768034],
        [11.355052085519442, 51.78796822768034],
        [11.355052085519442, 51.79015177231967],
        [11.352947914480557, 51.79015177231967],
        [11.352947914480557, 51.78796822768034]
      ]
    ]
  },
  id: {
    farm_id: test_farm_id,
    sector_id: test_sector_id
  },
  requests: {
    date_range: [2019-07-01, 2019-07-15],
    satellite_image: {
      max_size: 500
    },
    weather_data: {
      daily: [rain_sum, sunshine_duration],
      hourly: [snow_depth]
    },
    copernicus_data: {
      photosynthesis: [spatial_mean, temporal_mean],
      nitrogen: [spatial_mean, temporal_mean],
      canopy_water: [spatial_mean, temporal_mean],
      vegetation_cover: [spatial_mean, temporal_mean],
      leaf_area: [spatial_mean, temporal_mean],
      bands: [spatial_mean],
      composites: {
        ndvi:  [spatial_mean],
        savi:  [spatial_mean],
        evi:   [spatial_mean],
        evi2:  [spatial_mean],
        gndvi: [spatial_mean],
        msi:   [spatial_mean],
        ndmi:  [spatial_mean],
        nbr:   [spatial_mean],
        psri:  [spatial_mean],
        ari:   [spatial_mean],
        arvi:  [spatial_mean],
        sipi:  [spatial_mean],
        ndyi:  [spatial_mean],
        ags:   [spatial_mean]
      }
    }
  }
}


# model call:
{
  geometry: {
    type: Polygon,
    coordinates: [
      [
        [11.354947914480557, 51.78796822768034],
        [11.355052085519442, 51.78796822768034],
        [11.355052085519442, 51.79015177231967],
        [11.352947914480557, 51.79015177231967],
        [11.352947914480557, 51.78796822768034]
      ]
    ]
  },
  id: {
    farm_id: test_farm_id,
    sector_id: test_sector_id
  },
  requests: {
    date_range: [2019-07-01, 2019-07-15],
    satellite_image: {
      max_size: 500
    },
    weather_data: {
      daily: [rain_sum, sunshine_duration],
      hourly: [snow_depth]
    },
    copernicus_data: {
      photosynthesis: [spatial_mean, temporal_mean],
      nitrogen: [spatial_mean, temporal_mean],
      canopy_water: [spatial_mean, temporal_mean],
      vegetation_cover: [spatial_mean, temporal_mean],
      leaf_area: [spatial_mean, temporal_mean],
      bands: [spatial_mean],
      composites: {
        ndvi:  [spatial_mean],
        savi:  [spatial_mean],
        evi:   [spatial_mean],
        evi2:  [spatial_mean],
        gndvi: [spatial_mean],
        msi:   [spatial_mean],
        ndmi:  [spatial_mean],
        nbr:   [spatial_mean],
        psri:  [spatial_mean],
        ari:   [spatial_mean],
        arvi:  [spatial_mean],
        sipi:  [spatial_mean],
        ndyi:  [spatial_mean],
        ags:   [spatial_mean]
      }
    }
  }
}


# model call:
{
  geometry: {
    type: Polygon,
    coordinates: [
      [
        [11.354947914480557, 51.78796822768034],
        [11.355052085519442, 51.78796822768034],
        [11.355052085519442, 51.79015177231967],
        [11.352947914480557, 51.79015177231967],
        [11.352947914480557, 51.78796822768034]
      ]
    ]
  },
  id: {
    farm_id: test_farm_id,
    sector_id: test_sector_id
  },
  requests: {
    date_range: [2019-07-01, 2019-07-15],
    satellite_image: {
      max_size: 500
    },
    weather_data: {
      daily: [rain_sum, sunshine_duration],
      hourly: [snow_depth]
    },
    copernicus_data: {
      photosynthesis: [spatial_mean, temporal_mean],
      nitrogen: [spatial_mean, temporal_mean],
      canopy_water: [spatial_mean, temporal_mean],
      vegetation_cover: [spatial_mean, temporal_mean],
      leaf_area: [spatial_mean, temporal_mean],
      bands: [spatial_mean],
      composites: {
        ndvi:  [spatial_mean],
        savi:  [spatial_mean],
        evi:   [spatial_mean],
        evi2:  [spatial_mean],
        gndvi: [spatial_mean],
        msi:   [spatial_mean],
        ndmi:  [spatial_mean],
        nbr:   [spatial_mean],
        psri:  [spatial_mean],
        ari:   [spatial_mean],
        arvi:  [spatial_mean],
        sipi:  [spatial_mean],
        ndyi:  [spatial_mean],
        ags:   [spatial_mean]
      }
    }
  }
}


# model call:
{
  geometry: {
    type: Polygon,
    coordinates: [
      [
        [11.354947914480557, 51.78796822768034],
        [11.355052085519442, 51.78796822768034],
        [11.355052085519442, 51.79015177231967],
        [11.352947914480557, 51.79015177231967],
        [11.352947914480557, 51.78796822768034]
      ]
    ]
  },
  id: {
    farm_id: test_farm_id,
    sector_id: test_sector_id
  },
  requests: {
    date_range: [2019-07-01, 2019-07-15],
    satellite_image: {
      max_size: 500
    },
    weather_data: {
      daily: [rain_sum, sunshine_duration],
      hourly: [snow_depth]
    },
    copernicus_data: {
      photosynthesis: [spatial_mean, temporal_mean],
      nitrogen: [spatial_mean, temporal_mean],
      canopy_water: [spatial_mean, temporal_mean],
      vegetation_cover: [spatial_mean, temporal_mean],
      leaf_area: [spatial_mean, temporal_mean],
      bands: [spatial_mean],
      composites: {
        ndvi:  [spatial_mean],
        savi:  [spatial_mean],
        evi:   [spatial_mean],
        evi2:  [spatial_mean],
        gndvi: [spatial_mean],
        msi:   [spatial_mean],
        ndmi:  [spatial_mean],
        nbr:   [spatial_mean],
        psri:  [spatial_mean],
        ari:   [spatial_mean],
        arvi:  [spatial_mean],
        sipi:  [spatial_mean],
        ndyi:  [spatial_mean],
        ags:   [spatial_mean]
      }
    }
  }
}


# model call:
{
  geometry: {
    type: Polygon,
    coordinates: [
      [
        [11.354947914480557, 51.78796822768034],
        [11.355052085519442, 51.78796822768034],
        [11.355052085519442, 51.79015177231967],
        [11.352947914480557, 51.79015177231967],
        [11.352947914480557, 51.78796822768034]
      ]
    ]
  },
  id: {
    farm_id: test_farm_id,
    sector_id: test_sector_id
  },
  requests: {
    date_range: [2019-07-01, 2019-07-15],
    satellite_image: {
      max_size: 500
    },
    weather_data: {
      daily: [rain_sum, sunshine_duration],
      hourly: [snow_depth]
    },
    copernicus_data: {
      photosynthesis: [spatial_mean, temporal_mean],
      nitrogen: [spatial_mean, temporal_mean],
      canopy_water: [spatial_mean, temporal_mean],
      vegetation_cover: [spatial_mean, temporal_mean],
      leaf_area: [spatial_mean, temporal_mean],
      bands: [spatial_mean],
      composites: {
        ndvi:  [spatial_mean],
        savi:  [spatial_mean],
        evi:   [spatial_mean],
        evi2:  [spatial_mean],
        gndvi: [spatial_mean],
        msi:   [spatial_mean],
        ndmi:  [spatial_mean],
        nbr:   [spatial_mean],
        psri:  [spatial_mean],
        ari:   [spatial_mean],
        arvi:  [spatial_mean],
        sipi:  [spatial_mean],
        ndyi:  [spatial_mean],
        ags:   [spatial_mean]
      }
    }
  }
}


# model call:
{
  geometry: {
    type: Polygon,
    coordinates: [
      [
        [11.354947914480557, 51.78796822768034],
        [11.355052085519442, 51.78796822768034],
        [11.355052085519442, 51.79015177231967],
        [11.352947914480557, 51.79015177231967],
        [11.352947914480557, 51.78796822768034]
      ]
    ]
  },
  id: {
    farm_id: test_farm_id,
    sector_id: test_sector_id
  },
  requests: {
    date_range: [2019-07-01, 2019-07-15],
    satellite_image: {
      max_size: 500
    },
    weather_data: {
      daily: [rain_sum, sunshine_duration],
      hourly: [snow_depth]
    },
    copernicus_data: {
      photosynthesis: [spatial_mean, temporal_mean],
      nitrogen: [spatial_mean, temporal_mean],
      canopy_water: [spatial_mean, temporal_mean],
      vegetation_cover: [spatial_mean, temporal_mean],
      leaf_area: [spatial_mean, temporal_mean],
      bands: [spatial_mean],
      composites: {
        ndvi:  [spatial_mean],
        savi:  [spatial_mean],
        evi:   [spatial_mean],
        evi2:  [spatial_mean],
        gndvi: [spatial_mean],
        msi:   [spatial_mean],
        ndmi:  [spatial_mean],
        nbr:   [spatial_mean],
        psri:  [spatial_mean],
        ari:   [spatial_mean],
        arvi:  [spatial_mean],
        sipi:  [spatial_mean],
        ndyi:  [spatial_mean],
        ags:   [spatial_mean]
      }
    }
  }
}


# model call:
Create a file with the content below, e.g., request.json


---
request.json example
---

{
  geometry: {
    type: Polygon,
    coordinates: [
      [
        [11.354947914480557, 51.78796822768034],
        [11.355052085519442, 51.78796822768034],
        [11.355052085519442, 51.79015177231967],
        [11.352947914480557, 51.79015177231967],
        [11.352947914480557, 51.78796822768034]
      ]
    ]
  },
  id: {
    farm_id: test_farm_id,
    sector_id: test_sector_id
  },
  requests: {
    date_range: [2019-07-01, 2019-07-15],
    satellite_image: {
      max_size: 500
    },
    weather_data: {
      daily: [rain_sum, sunshine_duration],
      hourly: [snow_depth]
    },
    copernicus_data: {
      photosynthesis: [spatial_mean, temporal_mean],
      nitrogen: [spatial_mean, temporal_mean],
      canopy_water: [spatial_mean, temporal_mean],
      vegetation_cover: [spatial_mean, temporal_mean],
      leaf_area: [spatial_mean, temporal_mean],
      bands: [spatial_mean],
      composites: {
        ndvi:  [spatial_mean],
        savi:  [spatial_mean],
        evi:   [spatial_mean],
        evi2:  [spatial_mean],
        gndvi: [spatial_mean],
        msi:   [spatial_mean],
        ndmi:  [spatial_mean],
        nbr:   [spatial_mean],
        psri:  [spatial_mean],
        ari:   [spatial_mean],
        arvi:  [spatial_mean],
        sipi:  [spatial_mean],
        ndyi:  [spatial_mean],
        ags:   [spatial_mean]
      }
    }
  }
}

---
