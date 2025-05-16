########################################################################################################
# Author: Philipp Ertz, July 2024. Last edit: May 2025
#
# Functions required for the spatial interpolation process:
# draw_cond_GRF: drawing one conditional representation for the prediction locations from one Gaussian random field.
# draw_cond_large_GRF: iterative application of draw_cond_GRF, for drawing large fields from a conditional GRF (Sect. 5.6)
########################################################################################################

library("MASS") # for drawing from MV-norm distribution
library("geosphere") # for Haversine formula
library("geoR") # for matern function. Be careful to scale phi with 1/sqrt(3) and set kappa to 3/2
library("evd") # for sampling from Gumbel distribution

###############################
# Simple interpolation function
# #############################
draw_cond_GRF <- function(theta, # vector of observed values
                  r.new, # coordinates for predicted location, matrix of dim [J,2] for J prediction locations
                  r.given, # coordinated for given locations, matrixof dim [Mx2] for M known locations
                  alpha=0, # mean 
                  sigma=1, # sill
                  rho=100, # range
                  f=0, # scaling factor for altitude difference
                  z.new = NULL, # altitude of new locations, vector [J]
                  z.given = NULL # altitude of given locationsm vector [M]
                  ){
    # retrieve number of given and prediction locations
    if(is.vector(r.new)){r.new <- matrix(r.new, ncol=2)}
    N <- dim(r.given)[1]
    J <- dim(r.new)[1] 
    nu <- 3/2 # smoothness parameter of Matern-function
    R_e <- 6371000 # Earth's radius

    # calculate covariance matrix for given stations
    dist.given <- matrix(0,nrow=N, ncol=N)
    for (i in 1:N){
        dist.given[,i] <- sqrt( (distHaversine(r.given, r.given[i,], r=R_e)/1000)^2 + (f * abs(z.given-z.given[i]) /1000)^2)
    }
    K.given <- (sigma^2) * matern(dist.given, rho/sqrt(2*nu), nu)
    diag(K.given) <- diag(K.given) + 1e-5
    K.inv <- solve(K.given)

    # calculate covariance for new stations
    dist.new <- matrix(0,nrow=J, ncol=J)
    for (i in 1:J){
        dist.new[,i] <- sqrt((distHaversine(r.new, r.new[i,], r=R_e)/1000)^2 + (f * abs(z.new-z.new[i]) /1000)^2)
    }
    K.new <- sigma^2*matern(dist.new, rho/sqrt(nu*2), nu)

    # calculate cross-distance matrix
    dist.cross <- matrix(0,nrow=N,ncol=J)
    for (i in 1:J){
        dist.cross[,i] <- sqrt((distHaversine(r.given, r.new[i,], r=R_e)/1000)^2 + (f * abs(z.given-z.new[i]) /1000)^2)
    }
    K.cross <- sigma^2*matern(dist.cross,rho/sqrt(2*nu), nu)

    # calculate kriging estimate
    alpha.vec <- rep(alpha,max(N,J)) # expectation vectors
    theta.new.mean <- alpha.vec[1:J] + t(K.cross) %*% K.inv %*% (theta - alpha.vec[1:N])

    # calculate conditional covariance
    theta.cov <- K.new - t(K.cross) %*% K.inv %*% K.cross

    # draw new representations
    theta.new <- mvrnorm(mu=theta.new.mean, Sigma=theta.cov, tol=1e-1)

    # return predicted vector
    return(theta.new)
}

##############################################################
# Iterative interpolation procedure (Sect. 5.6 of manuscript)
##############################################################
## Function for iteratively draw large samples. Has to be performed several times to obtain an estimate of the mean and variance at each location
draw_cond_large_GRF <- function(r.given,
                                theta,
                                r.new,
                                alpha,
                                sigma,
                                rho,
                                f=0,
                                z.given=NULL,
                                z.new=NULL,
                                nbatch=500 # controls the number of prediction locations per iteration){
    N <- dim(r.new)[1]
    M <- dim(r.given)[1]
    nbins <- ceiling(N/nbatch)

    # save returned field
    theta.new <- vector("numeric", length=N)

    # loop over all bins
    print(paste("Start iterative drawing procedure with", nbins, "iterations and a batch size of", nbatch))
    pb = txtProgressBar(min=0, max=nbins, initial=1)
    for (i in 1:nbins){
        # select indices for bins
        if (i==nbins){
            ind <- ((i-1)*nbatch+1):N # make sure the index does not exceed the length of the vector
        } else {
            ind <- ((i-1)*nbatch+1):(i*nbatch)
        }

        # parse proper input data
        if (i==1){
            # take original values in the first round
            r.given.bin <- r.given
            theta.bin <- theta
            z.given.bin <- z.given
        } else {
            # add results from last fit in all subsequent rounds
            r.given.bin <- rbind(r.new[ind-nbatch,], r.given)
            theta.bin <- c(theta.new[ind-nbatch], theta)
            z.given.bin <- c(z.new[ind-nbatch], z.given)
        }
        r.new.bin <- r.new[ind,]
        z.new.bin <- z.new[ind]

        # draw from conditional GRF
        theta.new[ind] <- draw_cond_GRF(theta=theta.bin,
                                        r.given=r.given.bin,
                                        r.new=rn.bin,
                                        alpha=alpha,
                                        sigma=sigma,
                                        rho=rho,
                                        f=f,
                                        z.given=z.given.bin,
                                        z.new=zn.bin)

        # adjust progress progress bar
        setTxtProgressBar(pb,i)
    }
    # close progress bar
    close(pb)
                                
    return(theta.new)
}
