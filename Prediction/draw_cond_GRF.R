############################################################################################################################
# Author: Philipp Ertz, Institue of Geosciences, Meteorology Section, Bonn University (pertz@uni-bonn.de)
# Last updated: July 2024
# For further information see: ... insert paper...
# 
# Functions required for the prediction process fro SpatBHM:
# predict.gusts: function for sampling from Gumbel-distribution, provided covariate vectors and regression coefficients.
# draw_cond_GRF: drawing one representation for the prediction locations from one Gaussian random field.
# draw_cond_large_GRF: iterative application of draw_cond_GRF, for drawing large fields from a conditional GRF (Chapter 5.5)
#############################################################################################################################

library("MASS") # for drawing from MV-norm distribution
library("geosphere") # for Haversine formula
library("geoR") # for matern function. Be careful to scale phi with 1/sqrt(3) and set kappa to 3/2
library("evd") # for sampling from Gumbel distribution

# Function for sampling from conditional Gumbel distribution.
predict.gusts <- function(mus, # regression coefficients for location; numerical vector, 1 element per covariate
                            sigmas, # regression coefficients for scale; numerical vector 1 element per covariate
                            xmu, # covariate vectors for location, matrix with columns for each covariate
                            xsigma  # covariate vectors for scale, matrix with columns for each covaraite
                            )
{
    # Vectorized over covariate values, apply for each draw from the conditional distribution of the regression coefficient to obtain a sample.
    # Be careful to add VMEAN again, when fitting fx-VMEAN, as the function only returns the direct predictant
    if (dim(xmu)[1]!= dim(xsigma)[1]){stop("Covariate vectors must have the same length.")}
    N <- dim(xmu)[1]

    # calculate Gumbel parameters
    mu <- mus %*% xmu
    sigmax <- exp(sigmas %*% xsigma)

    # sammple from Gumbel distribution
    fx_predicted <- rgumbel(n=N, loc=mu, scale=sigma)
    return(fx_predicted)
}

# draw one representation of the conditional values at the prediction locations ("kriging")
draw_cond_GRF <- function(theta, # observed values
                  xnew, # new locations
                  xgiven, # old locations
                  alpha=0, # mean
                  sigma=1, # sill
                  rho=100, # range
                  f=0, # scaling factor for elevation offset
                  znew = NULL, # altitude of new locations
                  zgiven = NULL # altitude of given locations
                  ){
    # retrieve number of given and prediction locations
    if(is.vector(xnew)){xnew <- matrix(xnew, ncol=2)}
    N <- dim(xgiven)[1]     # number of observations
    J <- dim(xnew)[1]       # number of predictions
    nu <- 3/2               # smoothness parameter of Matern-function (fixed value)
    R_e <- 6371000          # Earth's radius

    # calculate distance matrix for given locations
    dist_given <- matrix(0,nrow=N, ncol=N)
    for (i in 1:N){
        dist_given[,i] <- distHaversine(xgiven, xgiven[i,], r=R_e)/1000 # use great circle distance
    }
    if (f){
        d_update <- matrix(0,nrow=N, ncol=N)
        for (i in 1:N){
            d_update[,i] <- f * abs(zgiven-zgiven[i]) /1000 # calculate altitude contribution
        }
        dist_given <- dist_given + d_update
    }

    # evaluate covariance function
    K_given <- (sigma^2) * matern(dist_given, rho/sqrt(2*nu), nu) # be careful with definition of parameters, as this is a general matern-class function
    #diag(K_given) <- diag(K_given) + 1e-5 # regularize
    K_inv <- solve(K_given) # invert

    # calculate covariance for new stations
    dist_new <- matrix(0,nrow=J, ncol=J)
    for (i in 1:J){
        dist_new[,i] <- distHaversine(xnew, xnew[i,], r=R_e)/1000
    }
    if (f){
        d_update <- matrix(0,nrow=J, ncol=J)
        for (i in 1:J){
            d_update[,i] <- f * abs(znew-znew[i]) /1000
        }
        dist_new <- dist_new + d_update
    }
    K_new <- sigma^2*matern(dist_new, rho/sqrt(nu*2), nu)

    # calculate cross-covariance matrix between given and new locations
    dist_cross <- matrix(0,nrow=N,ncol=J)
    for (i in 1:J){
        dist_cross[,i] <- distHaversine(xgiven, xnew[i,], r=R_e)/1000
    }
    if (f){
        d_update <- matrix(0,nrow=N,ncol=J)
        for (i in 1:J){
            d_update[,i] <- f * abs(zgiven-znew[i]) /1000
        }
        dist_cross <- dist_cross + d_update
    }
    K_cross <- sigma^2*matern(dist_cross,rho/sqrt(2*nu), nu)

    # calculate kriging estimate
    alpha_vec <- rep(alpha,max(N,J)) # expectation vectors
    theta_new_mean <- alpha_vec[1:J] + t(K_cross) %*% K_inv %*% (theta - alpha_vec[1:N])

    # calculate conditional covariance
    theta_cov <- K_new - t(K_cross) %*% K_inv %*% K_cross

    # draw new representations
    theta_new <- mvrnorm(mu=theta_new_mean, Sigma=theta_cov, tol=1e-1)

    # return predicted vector
    return(theta_new)
}

## Function for iteratively draw large samples. Has to be performed several times to obtain an estimate of the mean and variance at each location.
## Code breaks down numerically, when the elevation offset is included. We tested regularizations of the covariance matrices to no avail.
draw_cond_large_GRF <- function(xgiven,
                                theta,
                                xnew,
                                alpha,
                                sigma,
                                rho,
                                f=0,
                                zgiven=NULL,
                                znew=NULL,
                                nbatch=500,
                                resample=FALSE){
    N <- dim(xnew)[1]
    M <- dim(xgiven)[1]
    nbins <- ceiling(N/nbatch)

    # resample coordinates
    if (resample){
        print("resampling...")
        i_res <- sample(1:N, N, replace=FALSE)
        xnew <- xnew[i_res,]
        if(f){znew <- znew[i_res]}
    }

    # save returned field
    theta_new <- vector("numeric", length=N)

    # loop over all bins
    print(paste("Start iterative drawing procedure with ", nbins, "iterations and a batch size of ", nbatch))
    pb = txtProgressBar(min=0, max=nbins, initial=1)
    for (i in 1:nbins){
        print(i)
        # select indices for bins
        if (i==nbins){
            ind <- ((i-1)*nbatch+1):N # make sure the index does not exceed the end
        } else {
            ind <- ((i-1)*nbatch+1):(i*nbatch)
        }

        # parse proper input data
        if (i==1){
            # take original values in the first round
            xg_bin <- xgiven
            theta_bin <- theta
            if(f){zg_bin <- zgiven} else{zg_bin<-NULL} # only when f is provided
        } else {
            # add results from last fit in all subsequent rounds
            xg_bin <- rbind(xnew[ind-nbatch,], xgiven)
            theta_bin <- c(theta_new[ind-nbatch], theta)
            if (f){zg_bin <- c(znew[ind-nbatch], zgiven)} else {zg_bin<-NULL}# only when f is provided
        }
        xnew_bin <- xnew[ind,]
        if(f){znew_bin <- znew[ind]} else {znew_bin<-NULL}

        # draw from conditional GRF
        theta_new[ind] <- draw_cond_GRF(theta=theta_bin,
                                        xgiven=xg_bin,
                                        xnew=xnew_bin,
                                        alpha=alpha,
                                        sigma=sigma,
                                        rho=rho,
                                        f=f,
                                        zgiven=zg_bin,
                                        znew=znew_bin)

        # adjust progress progress bar
        setTxtProgressBar(pb,i)
    }
    # close progress bar
    close(pb)
    if(resample){theta_new <- theta_new[order(i_res)]} # resample back for output
    return(theta_new)
}
