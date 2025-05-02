######################################################################################################
# Author: Philipp Ertz, July 2024
#
# Some licence
#
# This script performs the spatial interpolation of a spatial bayesian hierarchical wind gust model
######################################################################################################

# Load kriging function: draw_cond_GRF() draws one conditional representation from a GRF specified in my way (including altitude scaling, matern-covariance, distance along great circles)
source("/automount/user/s6phertz/Dokumente/PhD/Code/Gust_code/Gust_code_pv/Prediction/draw_cond_GRF_euclidified.R")
#source("/automount/user/s6phertz/Dokumente/PhD/Code/Gust_code/gust_code.git/Prediction/draw_cond_GRF.R")


# tidyverse library for string manipulation (for retrieving execution info from model name)
library(stringr)

# pass model names as command line argument
models <- commandArgs(trailingOnly=TRUE)

# Project directory on the JupyterHub
hub_dir <- "/automount/hubhome/s6phertz/Coming_decade/Gust_paper/"

# list of eligible models, non-eligible models will be skipped
eligible_models <- c("Baseline_mu2",
                        "Baseline_vmean",
                        "Baseline_vmean_pred",
                        "Baseline_flat",
                        "SM0",
                        "SM_mu0",
                        "SM_mu0_mount",
                        "SM_mu0_sigma0",
                        "SM_mu0_sigma0_mount",
                        "SM_mu0_mu1",
                        "SM_mu0_mu1_mount",
                        "SM_mu0_sigma1",
                        "SM_mu0_mu2",
                        "SM_mu0_mu2_mount",
                        "SM_mu0_sigma2",
                        "SM_mu0_mu1_sigma0",
                        "SM_mu0_mu1_mu2",
                        "SM_mu0_mu1_mu2_mount",
                        "Baseline0",
                        "Baseline_mu2_mount",
                        "Baseline_optimal",
                        "Baseline_vmean_mount",
                        "SM_mu0_f",
                        "SM_mu0_mu1_f",
                        "SM_mu0_mu2_f",
                        "SM_mu0_sigma0_f",
                        "SM_mu0_mu1_mu2_f",
                        "SM_mu0_mu2_sigma0_f",
                        "SM_mu0_mt")

# translate parameter names between model name and Stan output files
pd <- data.frame(name=c("mu0","mu1","mu2","mu3","sigma0", "sigma1", "sigma2"), col = c("mu.1", "mu.2", "mu.3", "mu.4", "sigma.1", "sigma.2", "sigma.3"), row.names=T )
pd_back <- data.frame( col = c("mu.1", "mu.2", "mu.3", "mu.4", "sigma.1", "sigma.2", "sigma.3"),name=c("mu0","mu1","mu2","mu3","sigma0", "sigma1", "sigma2"), row.names=T )

for (model_name in models){
    if (!(model_name %in% eligible_models)){
        print(paste(model_name,"is not a valid model."))
        next
    }

    ###############
    # Preparation
    ###############

    # read station list
    labels <- c(unlist(unname(read.table("/automount/hubhome/s6phertz/Coming_decade/Gust_paper/data/used_station_ids.csv", skip=1, colClasses = c("character")))))

    # read meta data of stations (coordinates and altitude)
    stat_info <- read.table("/automount/hubhome/s6phertz/Coming_decade/Gust_paper/data/used_stations_109.csv", sep=",", header=TRUE)

     # give update as command line out put
    print(paste("Interpolating ", model_name, ". Prediction files will be saved in ", hub_dir, "Kriging/", sep=""))

    # Exclude mountain top stations, for models where it is necessary
    if (model_name %in% c("SM_mu0",
                        "SM_mu0_sigma0",
                        "SM_mu0_mu1",
                        "SM_mu0_sigma1",
                        "SM_mu0_mu2",
                        "SM_mu0_sigma2",
                        "SM_mu0_mu1_sigma0",
                        "SM_mu0_mu1_mu2")){
        flat_flag <- TRUE
    } else if(model_name %in% c("SM_mu0_f",
                                "SM_mu0_mount",
                                "SM_mu0_mu1_f",
                                "SM_mu0_mu1_mount",
                                "SM_mu0_mu2_f",
                                "SM_mu0_mu2_mount",
                                "SM_mu0_sigma0_f",
                                "SM_mu0_sigma0_mount",
                                "SM_mu0_mu1_mu2_f",
                                "SM_mu0_mu1_mu2_mount",
                                "SM_mu0_mu2_sigma0_f",
                                "SM_mu0_mt")){
        flat_flag <- FALSE
    } else {
        print(paste("Invalid model encountered for", model_name))
        next
    }

    if (flat_flag==TRUE){
        fit_ind = which(stat_info$Stationshoehe<=800)
        stat_info = stat_info[fit_ind,]
        labels = labels[fit_ind]
    }

    M <- length(labels)
    mrange <- 1:M

    # get spatial parameters from model name, exclude "_mount" and "_f"
    if (str_sub(model_name, start=str_length(model_name)-4)=="mount"){end_parstring <- str_length(model_name)-6} else if (str_sub(model_name, start=str_length(model_name))=="f"){end_parstring <- str_length(model_name)-2} else if (str_sub(model_name, start=str_length(model_name)-1)=="f") {} else {end_parstring<-str_length(model_name)-3}

    # read and parse parameter names
    model_pars <- pd[unlist(str_split(str_sub(model_name, start=4, end=end_parstring), "_")),]

    # get altitude inclusion from model name
    include_f <- str_sub(model_name, start= str_length(model_name))=="f" # check whether last letter is f

    # directory for loading fit results ans storing kriged values
    mod_dir <- paste(hub_dir, "Model_fits/",model_name, sep="")
    krig_dir <- paste(hub_dir, "Kriging/", sep="")

    # list of constant parameters
    const_pars <- pd$col[which(!(pd$col %in% model_pars))]

    # list of column names to read the spatial parameter representations
    nstat <- length(labels)
    npar <- length(model_pars)
    #spat_cols <- c()
    #for (i in 1:npar){
    #    for (j in 1:nstat){
     #       spat_cols <- c(spat_cols, paste("spat_pars.",i,".",j, sep=""))
    #    }
    #}

    # get column names for the GRF parameters
    grf_params <- c()
    for (p in 1:npar){
        grf_params <- c(grf_params, c(paste("expect_",pd_back[model_pars[p],], sep=""),
                                        paste("sills.",p, sep=""),
                                        paste("ranges.",p, sep="")))
    }

    #################################
    # Start looping over all CV-fits
    #################################

    # create progress bar
    pb = txtProgressBar(min=0, max=M, initial=1)

    for (i in 1:35){
        id <- labels[i]
        ind <- mrange[mrange!=i]
        J <- 1

        # find observation coordinates
        coord <- matrix(0, ncol=2, nrow=M-J)
        coord[,1] <- stat_info$geoLaenge[ind]
        coord[,2] <- stat_info$geoBreite[ind]
        z <- stat_info$Stationshoehe[ind]

        # find prediciton coordinates
        eval_coord <- matrix(0, ncol=2, nrow=J)
        eval_coord[,1] <- stat_info$geoLaenge[i]
        eval_coord[,2] <- stat_info$geoBreite[i]
        z_pred <- stat_info$Stationshoehe[i]

        # read model fit data
        fit_file <- paste(mod_dir, "/", model_name, "_", labels[i],"_fit.csv", sep="")

        # catch naming inconsistencies for the fit files
        if (model_name %in% c("SM_mu0", "SM_mu0_mu1", "SM_mu0_sigma0")){fit_file <- paste(mod_dir, "/", model_name, "_fit_", labels[i],".csv", sep="")}
        if (model_name=="SM_mu0_sigma1"){fit_file <- paste(mod_dir, "/", model_name, "_", labels[i],".csv", sep="")}

        # read model fit
        fit <- read.table(fit_file, sep=",", header=T)
        N <- dim(fit)[1]

        if (include_f){f_scale <- fit$f}
        kriging_values <- data.frame()

        ########################################
        # start looping over spatial parameters
        #######################################

        par_count <- 0 # counter for number of iteration

        for (par in 1:length(model_pars)){
            par_count <- par_count + 1

            # create parameter vectors for easier looping
            alpha <- fit[,paste("expect_", model_pars[par], sep="")]
            sigma <- fit[, paste("sills.", par, sep="")]
            rho <- fit[, paste("ranges.", par, sep="")]

            ###########################################
            # Start looping over Markov chains
            ###########################################

            for (n in 1:N){
                # get column names for fitted values
                spat_cols <- c()
                for (j in 1:(nstat-1)){
                    spat_cols <- c(spat_cols, paste("spat_pars.",par_count,".",j, sep=""))
                }

                # obtain fitted values
                theta_fit <- unname(unlist(fit[n,spat_cols]))

                # draw a conditional sample from the GRF
                if (include_f){
                    kriging_values[n,pd_back[model_pars[par],]] <- draw_cond_GRF(theta=theta_fit,
                                                                 xnew=eval_coord,
                                                                 xgiven=coord,
                                                                 alpha=alpha[n],
                                                                 sigma=sigma[n],
                                                                 rho=rho[n],
                                                                 f=f_scale[n],
                                                                 znew=z_pred,
                                                                 zgiven=z)
                } else {
                    kriging_values[n,pd_back[model_pars[par],]] <- draw_cond_GRF(theta=theta_fit,
                                                                 xnew=eval_coord,
                                                                 xgiven=coord,
                                                                 alpha=alpha[n],
                                                                 sigma=sigma[n],
                                                                 rho=rho[n])
                }
            }
        }
        for (par in const_pars){
            kriging_values[pd_back[par,]] = fit[paste("expect_", par, sep="")]
        }

        ###################
        # write to file
        ###################

        filename <- paste(krig_dir, "kriged_values_",id,"_", model_name,".csv" ,sep="")
        write.table(kriging_values, file=filename, col.names=TRUE, sep=",", row.names=FALSE)
        setTxtProgressBar(pb,i)
    }
    close(pb)
}
