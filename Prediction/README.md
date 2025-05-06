# Prediction notes

For all ConstMod versions and and LocMod, apply `predict_gusts.R` directly, to predict the model in cross-validation. The script will create an output directory with the provided model name in this directoy, that contains the prediction samples.

For all SpatBHM versions, first run `spatial_interpolation.R`, to interpolate the model to the prediction locations. Then run `predict_gusts.R`. The interpolation script will store the output under `Kriging/<model_name>/`

The scripts are run from the linux shell via
```shell
Rscript predict_gusts.R <model_names>
```
and

```shell
Rscript spatial_interpolation.R <model_names>
```
The scripts can predict multiple models in succession.
