# Wind-Gust-Postprocessing-from-COSMO-REA6-with-a-Spatial-Bayesian-Hierarchical-Extreme-Value-Model

We propose a probabilistic gust analysis for the region of Germany that can be interpolated to unobserved locations. To this end, we develop a spatial Bayesian hierarchical model (BHM) for the post-processing of surface maximum wind gusts from COSMO-REA6 (Bollmeyer et al. 2015). Our approach (SpatBHM) uses a non-stationary extreme value distribution for the gust observations at the top level, with parameters that vary according to a linear model using predictor variables from the COSMO-REA6 reanalysis. To capture spatial patterns in surface extreme wind gust behavior, the regression coefficients are modeled as 2-dimensional Gaussian random fields with a constant mean and an isotropic covariance function that depends only on the distance between locations. In addition, we include a scaled elevation offset in the distance calculation for the covariance function to account for differences in topography, which allows to include data from mountain top stations in the training process and utilize all available information. The model is compared against a model with spatially constant coefficients (ConstMod, LocMod). 

All models are fitted using Hamiltonian Monte Carlo techniques and are implemented using Stan (https://mc-stan.org). The prediction process including the conditional simulation of the GRF in space at unobserved locations is implemented using GNU licensed free software from the R Project for Statistical Computing (http://www.r-project.org).

## Workflow notes:

The entire pipeline was distributed across multiple tools and environments, as each phase required specific libraries. The data was first downloaded via shell scripts and preprocessed in a Jupyter notebook using `xarray`, `pandas`, and `numpy`. The preprocessed data sets used for model training are found in the subfolder `Data`, along with further details on the preprocessing in `Data/data_preprocessing.md`. The model training was also performed in a Jupyter notebook using the `pystan`-interface (Riddel et al. 2021). The model code is provided under `Stan_code`, along with a description and a Jupyter Notebook with examples for model training. Prediction and verification were performed in R scripts via the shell. The prediction-scripts, including the spatial interpolation procedure, are provided in `Prediction`, along with further explanations in `Prediction/prediction.md`. Please note, that the execution of the model-fitting can take up to several weeks of computation time and creates an output directory `model_fits` with subdirectories according to model name. Likewise, the running of `Prediction/spatial_interpolation.R` and `Prediction/predict_gusts.R` will create individual output directories.

To ensure reproducibility, all relevant scripts are available as separate files. A detailed guide for setting up the respective environments, installing the necessary libraries, and switching between the different tools is provided in the accompanying `README.md` document.

The data is provided in standardized formats (e.g., CSV, NetCDF), which can be used seamlessly across all tools. Additional instructions on adjusting file paths and configuring the environment are also included in the `README.md` document.

## Model names

Model names in this code differ from those in the manuscript. Please refer to the following table for correspondence:

| Name in paper     | Name in code          | Description                           |
|-------------------|-----------------------|---------------------------------------|
| ConstMod 1        | `Baseline0`           | most simple spatially constant model  |
| ConstMod 2        | `Baseline_mu2`        | as 1. + altitude predictor $\Delta_z$ |
| ConstMod 3        | `Baseline_vmean`      | but predicting $fx-V_\text{m}$        |
| ConstMod 4        | `Baseline_optimal`    | as 3. + predictor $V_\text{m}$        |
| LocMod            | `LocMod`              | as 4., without $\Delta_z$, trained on individual stations |
| SpatBHM 1         | `SM_mu0_f`            | Model with spatial $\mu^0$, scaled altitude offset in distance metric |
| SpatBHM 2a        | `SM_mu0_mu1_f`        | Model with spatial $\mu^0, \mu^1$, scaled altitude offset in distance metric |
| SpatBHM 2b        | `SM_mu0_mu2_f`        | Model with spatial $\mu^0,\mu^2$, scaled altitude offset in distance metric |
| SpatBHM 2c        | `SM_mu0_sigma0_f`     | Model with spatial $mu^0, \varsigma^0$, scaled altitude offset in distance metric |
| SpatBHM 3         | `SM_mu0_mu1_mu2_f`    | Model with spatial $mu^0,\mu^1,\mu2$, scaled altitude offset in distance metric |

## References:
  
  Bollmeyer, C., Keller, J. D., Ohlwein, C., Wahl, S., Crewell, S., Friederichs, P., Hense, A., Keune, J., Kneifel, S., Pscheidt, I., Redl, S. and Steinke, S.: Towards a high-resolution regional reanalysis for the European CORDEX domain, Q. J. R. Meteorol. Soc., 141, 1–15, https://doi.org/10.1002/qj.2486, 2015.

  Climate Data Center: Observations Germany, https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/,
last access: 31 October 2023, 2023
  
  R Development Core Team: R: A language and environment for statistical computing, R Foundation for Statistical Computing, Vienna,
Austria, http://www.R-project.org, ISBN 3-900051-07-0, 2010.
  
  Riddell, A., Hartikainen, A., and Carter, M.: pystan (3.0.0), PyPI, 2021.

  Stan Development Team: Stan: A probabilistic programming language, Stan Development Team, https://mc-stan.org/, version 2.33, 2022.
