// Fit a Gaussian random field to the fields
functions {
    // euclidean distance metric, compsed of haversinde formla and scaled elevation offset
    real gr_circle_d(vector x, vector y, real hx, real hy, real f) {
        real r = 6371; // in km,  mean earth radius
        vector[2] x_rad = x*pi()/180;
        vector[2] y_rad = y*pi()/180;
        real z_diff = abs(hy-hx)/1000*f;  //difference in altitude, scaled with f to contribute to distance
        real d_gc = 2 * r * asin(sqrt( sin( (y_rad[2]-x_rad[2]) /2 )^2 + cos(x_rad[2]) * cos(y_rad[2]) * sin( (y_rad[1]-x_rad[1]) /2)^2 ));
        real d = sqrt(d_gc^2 + z_diff^2);
        return d;
    }
    // matern 3/2 kernel for general distance measures d; here I will use the great circle distance above
    real matern32_general(real distance, real sigma, real length_scale) {
        real k = square(sigma) * (1 + sqrt(3) * distance / length_scale) * exp( - sqrt(3) * distance / length_scale );
        return k;
    }
}
data {
    int<lower=1> N; // number of observations
    vector[N] y; // values to be fitted
    array[N] vector[2] coord; // coordinates of observations
    vector[N] z_stat; // station altitude
}
parameters {
    real mu; //mean of GRF
    real<lower=0> sill; //sill/process variance of covariance function
    real<lower=0> range; //range/length scale of covariance function

    // scaling factor for altitude difference in distance metric
    real<lower=0> f;
}
model {
    // calculate covariance matrix
    real dist;
    matrix[N,N] K;
    matrix[N,N] L_K;

    // define predictive K matrix
    for (i in 1:N) {
        K[i,i] = square(sill);
        for (j in (i+1):N) {
            dist = gr_circle_d(coord[i],coord[j], z_stat[i], z_stat[j], f);
            K[i,j] = matern32_general(dist, sill, range); // range goes in in km
            K[j,i] = K[i,j];
        }

    }
    L_K = cholesky_decompose(K);

    // priors for GP parameters
    target += normal_lpdf(mu|8,3);
    target += inv_gamma_lpdf(sill|2,2);
    target += inv_gamma_lpdf(range|1.05,50);
    target += gamma_lpdf(f|15,0.15);

    row_vector[N] mu_vec;
    mu_vec = rep_row_vector(mu,N);
    target += multi_normal_cholesky_lpdf(y|mu_vec, L_K);
}
