// Stan code for the spatially ConstBHM model without a variable number of covariates. Prior parameters and covariates are specified from the outside.
data {
    int<lower=1> N; // number of data pairs of covariate and predictand
    vector[N] y; //predictand, either fx or fx-VMEAN
    int<lower=1> Mmu; // number of covariates for 
    int<lower=1> Msigma; // number of covariates for scale
    matrix[N, Mmu] xmu; // covariates for location, e.g. VMAX or VMEAN or station height
    matrix[N, Msigma] xsigma; // covariates for scale
    vector[Mmu] mumeans;
    vector[Mmu] muvariances;
    vector[Msigma] sigmameans;
    vector[Msigma] sigmavariances;
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
