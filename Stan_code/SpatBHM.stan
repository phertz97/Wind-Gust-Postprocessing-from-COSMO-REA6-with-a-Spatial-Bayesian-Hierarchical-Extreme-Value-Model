// Stan Code for the SpatBHM
// Author: Philipp Ertz, Institute of Geosciences, Meteorology Section, Bonn University (pertz@uni-bonn.de)

// This model can include other predictors and switch between spatial and non-spatial parameters, specified via input data
functions {
    // haversine formula for great circle distance including scaled elevation offset
    real gr_circle_d(vector x, vector y, real hx, real hy, real f) {
        real r = 6371; // earth's radius
        vector[2] x_rad = x*pi()/180;
        vector[2] y_rad = y*pi()/180;
        real z_diff = abs(hy-hx)/1000*f;  //difference in altitude, scaled with f to contribute to distance
        real d_gc = 2 * r * asin(sqrt( sin( (y_rad[2]-x_rad[2]) /2 )^2 + cos(x_rad[2]) * cos(y_rad[2]) * sin( (y_rad[1]-x_rad[1]) /2)^2 )) + z_diff;
        real d =sqrt(d_gc^2 + z_diff^2):
        return d;
    }
    // matern 3/2 covariance function
    real matern32_general(real distance, real sigma, real length_scale) {
        real k = square(sigma) * (1 + sqrt(3) * distance / length_scale) * exp( - sqrt(3) * distance / length_scale );
        return k;
    }
}
data {
    // sample size
    int<lower=1> N; // number of observations
    int<lower=1> M; // number of stations

    // weather station data
    array[M] vector[2] coord; // coordinates (lat/lon)
    vector[M] z_stat; // station altitude
    vector[M] z_grid; // model topography

    // predictand
    vector[N] y;

    // predictors
    int<lower=1> nxmu; // number of covariates for mu
    int<lower=1> nxsig; // number of covariates for sigma
    matrix[nxmu,N] Xmu; // covariates for mu, one column per predictor variable
    matrix[nxsig,N] Xsigma; // covariates for sigma

    // assignment vector data point to weather station
    array[N] int<lower=1, upper=M> nn; // vector assigning data point to station

    // select spatial parameters, binary/boolean vector
    array[nxmu] int<lower=0, upper=1> mu_s; // which mu coefficients are modelled spatially? 1=yes, 0=no
    array[nxsig] int<lower=0, upper=1> sigma_s; // which sigma coefficients are modelled spatially? 1=yes, 0=0

    // parameters for GRF mean priors (Gaussian)
    array[nxmu] vector[2] mu_prior; // one vector for each covariate for mu, vector contains location and scale
    array[nxsig] vector[2] sigma_prior; // one vector for each covariate for sigma, vector contains location and scale

    // covariance parameter priors (inverse Gamma)
    array[sum(mu_s)+ sum(sigma_s)] vector[2] sills_prior; // one vector for each spatial field, vector contains shape and scale
    array[sum(mu_s) + sum(sigma_s)] vector[2] ranges_prior; // one vector for each spatial field, vector contains shape and scale

    // prior for the altitude scaling factor f_z (Gamma)
    vector[2] f_prior; // shape and inverse scale
}
transformed data {
    array[M] vector[2] x; // transform coordinate array to vectors
    for (m in 1:M){
        x[m] = to_vector(coord[m,]);
    }
    vector[M] z_diff = z_stat-z_grid; // calculate Delta z
}
parameters {
    // factor for scaling altitude difference in distance metric
    real<lower=0> f;

    // means of spatially constant fields
    vector[nxmu] expect_mu;
    vector[nxsig] expect_sigma;

    // GRF parameters (sills and ranges)
    vector<lower=0>[sum(mu_s) + sum(sigma_s)] sills;
    vector<lower=0>[sum(mu_s) + sum(sigma_s)] ranges;

    // store representations of spatial fields, array of vectors for each spatial field
    array[sum(mu_s) + sum(sigma_s)] vector[M] spat_pars;
}
model {
    // define helper variables
    real dist;

    // mu vectors
    matrix[M,nxmu] mu;

    // sigma vectors
    matrix[M,nxsig] sigma;

    // start modeling
    // prior for the altitude factor
    f ~ gamma(f_prior[1],f_prior[2]);

    // counter for spatial fields
    int spat_count;
    spat_count = 0;

    // model mu^j
    for (k in 1:nxmu){
        // define prior for GRF mean/spatially constant parameter
        expect_mu[k] ~ normal(mu_prior[k][1], mu_prior[k][2]);
        if (mu_s[k]){
            spat_count += 1;
            // draw GRF parameters from their priors
            sills[spat_count] ~ inv_gamma(sills_prior[spat_count][1], sills_prior[spat_count][2]);
            ranges[spat_count] ~ inv_gamma(ranges_prior[spat_count][1], ranges_prior[spat_count][2]);

            // store covariance matrix and its Cholesky decomposition
            matrix[M,M] K;
            matrix[M,M] L_K;

            // define predictive K matrix
            for (i in 1:M) {
                K[i,i] = square(sills[spat_count]);
                for (j in (i+1):M) {
                    dist = gr_circle_d(to_vector(x[i]),to_vector(x[j]), z_stat[i], z_stat[j], f);
                    K[i,j] = matern32_general(dist, sills[spat_count], ranges[spat_count]); // range goes in in km
                    K[j,i] = K[i,j];
                }
            }
            L_K = cholesky_decompose(K);

            vector[M] mu_vec = rep_vector(expect_mu[k], M);
            spat_pars[spat_count] ~ multi_normal_cholesky(mu_vec,L_K);
            mu[,k] = spat_pars[spat_count];

        } else {
            mu[,k] = rep_vector(expect_mu[k], M);
        }
    }

    // model sigma^j
    for (k in 1:nxsig){
        // define prior for GRF mean/spatially constant parameter
        expect_sigma[k] ~ normal(sigma_prior[k][1], sigma_prior[k][2]);
        if (sigma_s[k]){
            spat_count += 1;
            // draw GRF parameters from their priors
            sills[spat_count] ~ inv_gamma(sills_prior[spat_count][1], sills_prior[spat_count][2]);
            ranges[spat_count] ~ inv_gamma(ranges_prior[spat_count][1], ranges_prior[spat_count][2]);

            // store covariance matrix and its Cholesky decomposition
            matrix[M,M] K;
            matrix[M,M] L_K;

            // define predictive K matrix
            for (i in 1:M) {
                K[i,i] = square(sills[spat_count]);
                for (j in (i+1):M) {
                    dist = gr_circle_d(to_vector(x[i]),to_vector(x[j]), z_stat[i], z_stat[j], f);
                    K[i,j] = matern32_general(dist, sills[spat_count], ranges[spat_count]); // range goes in in km
                    K[j,i] = K[i,j];
                }
            }
            L_K = cholesky_decompose(K);

            vector[M] mu_vec = rep_vector(expect_sigma[k], M);
            spat_pars[spat_count] ~ multi_normal_cholesky(mu_vec, L_K);
            sigma[,k] = spat_pars[spat_count];

        } else {
            sigma[,k] = rep_vector(expect_sigma[k],M);
        }
    }

    // construct parameter vectors of length N of mu and sigma vectors for sampling
    vector[N] mu_vector;
    vector[N] sigma_vector;

    // work with station-data assignment vector nn and loop over all data points to construct the Gumbel-parameter vectors
    for (i in 1:N){
        mu_vector[i] = to_row_vector(mu[nn[i],])*to_vector(Xmu[,i]);
        sigma_vector[i] = to_row_vector(sigma[nn[i],])*to_vector(Xsigma[,i]);
    }

    // model wind gust distribution
    y ~ gumbel(mu_vector, exp(sigma_vector));
}
