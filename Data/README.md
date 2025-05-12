# Preprocessing of Data

## 1. Wind gust observations $fx$

The data files were accessed from the DWD's climate data center under https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/extreme_wind/historical/, last access on 31 October 2023.
Hours from 13-18 UTC, months from May to October and years between 2001 and 2018 were selected in using the `pandas` in Python. 

The full observational data set was reduced to a consistent set data wihtout `NaN`-values in several steps:

1. We excluded stations with less than 90% complete observational time series in our investigation window

2. An iterative data cleaning algorithm was used to generate the largest possible and most consistent dataset. 
The algorithm automatically identified and removed days with incomplete data (<4 observations from 13-18 UTC) at too many stations and stations with imcomplete data at too many of the remaining days.
Due to the automated nature of the process, the exact decision-making for individual days/stations cannot be fully retraced today. 
Given the small autocorrelation in wind gust data, this cleaning procedure is unlikely to significantly impact the modeling of underlying wind patterns. 
The step was taken to eliminate data errors and ensure a robust dataset, that is consistent in space and time.

After the creation of a consistent dataset, we caculated daily maxima using `pandas`. 
The final results can be found under `fx_training.csv`, containing the odd-numbered years, and `fx_evaluation.csv`, containing the even-numbered years.

## 1.1 Station metadata, coordinates & altitude

The station metadata used in this study is based on publicly available data from the German Weather Service (DWD, https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/extreme_wind/historical/FX_Stundenwerte_Beschreibung_Stationen.txt). The original file was adapted as follows:

- Special characters (e.g., German umlauts) were replaced with ASCII equivalents for compatibility reasons (e.g., “München” → “Muenchen”).
- The full list of stations was reduced in the data-cleaning procedure described under 1.1.

Geographic coordinates and elevation data were not modified. The final list of stations used is included in this repository under `used_stations.csv`

## 2. COSMO-REA6 variables

The predictor data from COSMO-REA6 used in this study is obtained from the DWD open data server https://opendata.dwd.de/climate_environment/REA/COSMO-REA6/. 
The data was downloaded as `.grb`-files and processed using `xarray` in Python.

### 2.1 VMAX_10M

The VMAX_10M data is obtained from  https://opendata.dwd.de/climate_environment/REA/COSMO_REA6/hourly/2D/VMAX_10M/, last access on 31 October 2023. 
For each SYNOP station, we retrieved a matching time series from the grid cell closest to the geographical coordinates provided in `used_stations.csv`, measured by the great-circle distance.
We took the hours 13-18 UTC and calculated caily maxima. The data are stored in a `.csv`-file using the same data structure as for the SYMOP-observations. 
The final preprocessed datasets can be found under `vmax_training.csv`, containing the odd-numbered years, and `vmax_evaluation.csv`, containing the even-numbered years.

### 2.2 VMEAN_10M

In order to obtain the values for VMEAN_2M, the horizontal wind components U_10M and V_10M were obtained from https://opendata.dwd.de/climate_environment/REA/COSMO_REA6/hourly/2D/U_10M/ 
& https://opendata.dwd.de/climate_environment/REA/COSMO_REA6/hourly/2D/V_10M/, last access: 31 October 2023.
For each SYNOP station, we retrieved a matching time series from the grid cell closest to the geographical coordinates provided in `used_stations.csv`, measured by the great-circle distance.
Then, we selected the hours 13-18 UTC and proceeded to calculate the mean wind via

$V_\textrm{m} = \sqrt{U^2+V^2}$

and calculated caily averages. The data are stored in a `.csv`-file using the same data structure as for the SYMOP-observations. 
The final preprocessed datasets can be found under `vmean_training.csv`, containing the odd-numbered years, and `vmean_evaluation.csv`, containing the even-numbered years.

### 2.3 HSURF

The surface elevation data from COSMO-REA6 is obtained from https://opendata.dwd.de/climate_environment/REA/COSMO_REA6/constant/COSMO_REA6_CONST_withOUTsponge.nc.bz2, last access on 31 October 2023. 
For each SYNOP station, we retrieved the variable `HSURF` from the grid cell closest to the geographical coordinates provided in `used_stations.csv`, measured by the great-circle distance.
The values used in this study are provided as `COSMO-REA6_HSURF_stations.csv`.
