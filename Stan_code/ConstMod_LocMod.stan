// Author: Philipp Ertz, Institute of Geosciences, Meteorology Section, Bonn University (pertz@uni-bonn.de)
// For further information see: ... insert paper ...

// Stan code for the spatially ConstBHM model without a variable number of covariates. Prior parameters and covariates are specified from the outside.
data {
    // model specifications
    int<lower=1> N;              // observation/forecast sample size
    vector[N] y;                 // predictand, either fx or fx-VMEAN
    int<lower=1> Mmu;            // number of covariates for lcoation
    int<lower=1> Msigma;         // number of covariates for scale
    matrix[N, Mmu] xmu;          // predictors for location, e.g. VMAX or VMEAN or station height
    matrix[N, Msigma] xsigma;    // predictors for scale

    // prior parameters (normal)
    vector[Mmu] mumeans;          // mean of location regression coefficients
    vector[Mmu] muvariances;      // variance of location regression coefficients
    vector[Msigma] sigmameans;    // mean of scale regression coefficients
    vector[Msigma] sigmavariances;// variance of scale regression coefficients
}
parameters {
    //regression parameters for location
    vector[Mmu] mu;
    
    // regression parameters for scale
    vector[Msigma] sigma;
}
model {
    // priors 
    for (i in 1:Mmu){mu[i] ~ normal(mumeans[i],muvariances[i]);}
    for (j in 1:Msigma){sigma[j] ~ normal(sigmameans[j],sigmavariances[j]);}

    // actual model
    y ~ gumbel(xmu*mu, exp(xsigma*sigma));    
}
