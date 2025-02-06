# Wind-Gust-Postprocessing-from-COSMO-REA6-with-a-Spatial-Bayesian-Hierarchical-Extreme-Value-Model

We propose a probabilistic gust analysis for the region of Germany that can be interpolated to unobserved locations. To this end, we develop a spatial Bayesian hierarchical model (BHM) for the post-processing of surface maximum wind gusts from COSMO-REA6. Our approach (SpatBHM) uses a non-stationary extreme value distribution for the gust observations at the top level, with parameters that vary according to a linear model using predictor variables from the COSMO-REA6 reanalysis. To capture spatial patterns in surface extreme wind gust behavior, the regression coefficients are modeled as 2-dimensional Gaussian random fields with a constant mean and an isotropic covariance function that depends only on the distance between locations. In addition, we include a scaled elevation offset in the distance calculation for the covariance function to account for differences in topography, which allows to include data from mountain top stations in the training process and utilize all available information. The model is compared against a model with spatially constant coefficients (ConstMod, LocMod). 

All models are fitted using Hamilonian Monte Carlo techniques and are implemented using Stan (https://mc-stan.org). The prediction process including the conditional simulation of the GRF in space at unobserved locations is implemented using GNU licensed free software from the R Project for Statistical Computing (http://www.r-project.org).



## References:
***
Stan Development Team. 2025. Stan Reference Manual, 2.33. https://mc-stan.org

## Disclaimer:
This software uses R packages that are licensed under GPLv3. These packages are not included in this repository and must be installed separately. The code in this repository is licensed under the MIT license.
