# Wind-Gust-Postprocessing-from-COSMO-REA6-with-a-Spatial-Bayesian-Hierarchical-Extreme-Value-Model

The aim of this study is to provide a probabilistic gust analysis for the region of Germany that is calibrated with
station observations and with an interpolation to unobserved locations. To this end, we develop a spatial Bayesian hierarchical
model (BHM) for the post-processing of surface maximum wind gusts from the COSMO-REA6 reanalysis. Our approach uses
a non-stationary extreme value distribution for the gust observations at the top level, with parameters that vary according to a
linear model using COSMO-REA6 predictor variables. To capture spatial patterns in surface extreme wind gust behavior, the5
regression coefficients are modeled as 2-dimensional Gaussian random fields with a constant mean and an isotropic covariance
function that depends only on the distance between locations. In addition, we include an elevation offset in the distance metric
for the covariance function to account for differences in topography. This allows us to include data from mountaintop stations
in the training process and to utilize all available information. The training of the BHM is carried out with an independent
data set from which the data at the station to be predicted are excluded. We evaluate the spatial prediction performance at the10
withheld station using Brier score and quantile score, including their decomposition, and compare the performance of our BHM
to climatological forecasts and a non-hierarchical, spatially constant baseline model. This is done for 109 weather stations in
Germany.

All models are fitted using Hamiltonian Monte Carlo techniques and are implemented using Stan (https://mc-stan.org). The prediction process including the conditional simulation of the GRF in space at unobserved locations is implemented using GNU licensed free software from the R Project for Statistical Computing (http://www.r-project.org).

## Contents of this repository

- `Data`: contains the used data sets with a short description of sources and preprocessing methods
- `Stan_code`: contains model code for Stan
- `Model_training_examples.ipynb`: a JupyterNotebook containing python code examples for model training with pystan
- `Model_fits`: an output directory for the model fits in `Model_training_examples.ipynb`
- `Prediction`: contains R-scripts used for the spatiotemporal predictions from SpatBHM and the reference models

## Workflow notes & Computational requirements.

The entire pipeline was distributed across multiple tools and environments. The preprocessed data sets used for model training are found in the subdirectory `Data`, along with further details on the preprocessing. The model training was performed in a Jupyter notebook using the `pystan`-interface (Riddel et al. 2021). The model code is provided under `Stan_code` and a Jupyter Notebook with examples for its application is found under `Model_training_examples.ipynb`. Please note, that the execution of the model-fitting can take up to several weeks of computation time and creates an output directory `model_fits` with subdirectories according to model name. Prediction and verification were performed in R scripts via the shell. The prediction-scripts, including the spatial interpolation procedure, are provided in `Prediction`, along with further explanations in `Prediction/prediction.md`. Likewise, the running of `Prediction/spatial_interpolation.R` and `Prediction/predict_gusts.R` will create individual output directories.

To ensure reproducibility, all relevant scripts are available as separate files. A guide for switching between the different tools is provided in the accompanying `README.md` document.

The data is provided in standardized formats (e.g., CSV, NetCDF), which can be used seamlessly across all tools. Additional instructions on adjusting file paths and configuring the environment are also included in the `README.md` document.

The software versions that were used to run the analyses are the following:

Python (3.9.18)
- `numpy` (1.23.5)
- `pandas` (2.1.0)
- `xarray` (2024.7.0)
- `pystan` (3.9.0)
- `nest_asyncio` (1.5.6)

R (4.2.2)
- `MASS` (7.3.58.1)
- `geosphere` (1.5.18)
- `geoR` (1.9.4)
- `evd` (2.3.6.1)
- `stringr` (1.5.1)
- `reliabilitydiag` (0.2.1)
- `isotone` (1.1.1)

## Model names

Model names in this code differ from those in the manuscript. Please refer to the following table for correspondence:

| Name in paper     | Name in code          | Description                           |
|-------------------|-----------------------|---------------------------------------|
| ConstMod 1        | `Baseline0`           | spatially constant model, using $V_\text{max}$ as sole predictor |
| ConstMod 2        | `Baseline_mu2`        | as 1. + altitude predictor $\Delta_z$ for location $\mu$ |
| ConstMod 3        | `Baseline_vmean`      | but predicting $fx-V_\text{m}$        |
| ConstMod 4        | `Baseline_optimal`    | as 3. + predictor $V_\text{m}$        |
| LocMod            | `LocMod`              | as 4., without $\Delta z$, trained on individual stations |
| SpatBHM 1         | `SM_mu0_f`            | Model with spatial $\mu^0$, scaled altitude offset in distance metric |
| SpatBHM 2a        | `SM_mu0_mu1_f`        | Model with spatial $\mu^0, \mu^1$, scaled altitude offset in distance metric |
| SpatBHM 2b        | `SM_mu0_mu2_f`        | Model with spatial $\mu^0,\mu^2$, scaled altitude offset in distance metric |
| SpatBHM 2c        | `SM_mu0_sigma0_f`     | Model with spatial $\mu^0, \varsigma^0$, scaled altitude offset in distance metric |
| SpatBHM 3         | `SM_mu0_mu1_mu2_f`    | Model with spatial $\mu^0,\mu^1,\mu^2$, scaled altitude offset in distance metric |

For further details, please refer to the manuscript, Sect. 3. & 4.1

## R libraries for verification

We used the library `reliabilitydiag` (Dimitriadis et al. 2021) for the calculation of the Brier score. For the calculation and decomposition of the quantie score (Gneiting et al. 2023), we used MIT-licensed code shared by Daniel Wolffram under https://github.com/dwolffram/replication-ARSIA2023/tree/main/R.

## References
  
  Bollmeyer, C., Keller, J. D., Ohlwein, C., Wahl, S., Crewell, S., Friederichs, P., Hense, A., Keune, J., Kneifel, S., Pscheidt, I., Redl, S. and Steinke, S.: Towards a high-resolution regional reanalysis for the European CORDEX domain, Q. J. R. Meteorol. Soc., 141, 1–15, https://doi.org/10.1002/qj.2486, 2015.

  Climate Data Center: Observations Germany, https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/,
last access: 31 October 2023, 2023

  Dimitriadis, T., Gneiting, T., and Jordan, A. I.: Stable reliability diagrams for probabilistic classifiers, Proc. Natl. Acad. Sci. U.S.A., 118,
e2016191 118, https://doi.org/10.1073/pnas.2016191118, 2021

  Gneiting, T., Wolffram, D., Resin, J., Kraus, K., Bracher, J., Dimitriadis, T., Hagenmeyer, V., Jordan, A. I., Lerch, S., Phipps, K., and Schienle,
M.: Model diagnostics and forecast evaluation for quantiles, Ann. Rev. Stat. Appl., 10, 597–621, https://doi.org/10.1146/annurev-statistics-775
032921-020240, 2023.
  
  R Development Core Team: R: A language and environment for statistical computing, R Foundation for Statistical Computing, Vienna,
Austria, http://www.R-project.org, ISBN 3-900051-07-0, 2010.
  
  Riddell, A., Hartikainen, A., and Carter, M.: pystan (3.9.0), https://pystan.readthedocs.io/en/latest/, 2021.

  Stan Development Team: Stan: A probabilistic programming language, Stan Development Team, https://mc-stan.org/, version 2.33, 2022.
