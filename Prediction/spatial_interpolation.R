######################################################################################################
# Author: Philipp Ertz, July 2024; Laste edited: May 2025
#
# This script performs the spatial interpolation of SpatBHM in cross-validation
######################################################################################################

# Load kriging-function: draw_cond_GRF() draws one conditional representation from a GRF
source("Prediction/draw_cond_GRF.R")

# tidyverse library for string manipulation (for retrieving execution info from model name)
library(stringr)

# pass model names as command line argument
models <- commandArgs(trailingOnly=TRUE)

# project directory
dir.path <- "./"

# list of eligible models, non-eligible models will be skipped
eligible_models <- c("Baseline_mu2",
                        "Baseline_vmean",
                        "Baseline0",
                        "Baseline_optimal",
                        "SM_mu0_f",
                        "SM_mu0_mu1_f",
                        "SM_mu0_mu2_f",
                        "SM_mu0_sigma0_f",
                        "SM_mu0_mu1_mu2_f")

# translate parameter names between model name and Stan output files
pd <- data.frame(name=c("mu0","mu1","mu2","mu3","sigma0", "sigma1", "sigma2"), col = c("mu.1", "mu.2", "mu.3", "mu.4", "sigma.1", "sigma.2", "sigma.3"), row.names=T )
pd_back <- data.frame( col = c("mu.1", "mu.2", "mu.3", "mu.4", "sigma.1", "sigma.2", "sigma.3"),name=c("mu0","mu1","mu2","mu3","sigma0", "sigma1", "sigma2"), row.names=T )

for (model_name in models){
  if (!(model_name %in% eligible_models)){
    print(paste(model_name,"is not a valid model."))
    next
  }

  # directory for loading fit results ans storing kriged values
  fit.dir <- paste(dir.path, "Model_fits/",model_name, sep="")
  output.dir <- paste(hub_dir, "Kriging/",model_name, sep="")
  if (!file.exists(out.dir)){
    dir.create(output.dir, recursive=TRUE)
  }
  print(paste("Interpolating ", model_name, ". Prediction files will be saved in ", output.dir, sep=""))

  ###############
  # Preparation
  ###############

  # get station IDs, coordinates and elevations
  stat_info <- read.table(file = paste(dir.path, "Data/used_stations.csv", sep=""), sep=",", header=TRUE, colClasses= c(rep("character",3), rep("numeric",3), rep("character",2)))
  stat.ids <- str_pad(stat_info$Stations_id, width = 5, side = "left", pad = "0")
  z_station <- stat_info$Stationshoehe

  # get spatial parameters from model name, exclude "_f"
  end_parstring <- str_length(model_name)-2
  model_pars <- pd[unlist(str_split(str_sub(model_name, start=4, end=end_parstring), "_")),]

  # list of constant parameters
  const_pars <- pd$col[which(!(pd$col %in% model_pars))]

  # list of column names to read the spatial parameter representations
  M <- length(stat.ids)
  npar <- length(model_pars)
  mrange <- 1:M

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

  for (i in mrange){
    # select fitted station ids
    ind_cv <- mrange[mrange!=i] # given stations
    J <- 1 # number of removed staions

    # find observation coordinates
    coord <- matrix(0, ncol=2, nrow=M-J)
    coord[,1] <- stat_info$geoLaenge[ind_cv]
    coord[,2] <- stat_info$geoBreite[ind_cv]
    z <- stat_info$Stationshoehe[ind_cv]

    # find prediciton coordinates
    eval_coord <- matrix(0, ncol=2, nrow=J)
    eval_coord[,1] <- stat_info$geoLaenge[i]
    eval_coord[,2] <- stat_info$geoBreite[i]
    z_pred <- stat_info$Stationshoehe[i]

    # read model fit data
    fit_file <- paste(fit.dir, "/", model_name, "_", stat.ids[i],"_fit.csv", sep="")

    # read model fit
    fit <- read.table(fit_file, sep=",", header=T)
    N <- dim(fit)[1]
    f_scale <- fit$f
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
        for (j in 1:(length(ind_cv))){
          spat_cols <- c(spat_cols, paste("spat_pars.",par_count,".",j, sep=""))
        }

        # obtain sample of fitted values for the iteration
        theta_fit <- unname(unlist(fit[n,spat_cols]))

        # draw a conditional sample from the GRF
        kriging_values[n,pd_back[model_pars[par],]] <- draw_cond_GRF(theta=theta_fit,
                                                         r.new=eval_coord,
                                                         r.given=coord,
                                                         alpha=alpha[n],
                                                         sigma=sigma[n],
                                                         rho=rho[n],
                                                         f=f_scale[n],
                                                         z.new=z_pred,
                                                         z.given=z)
      }
    }
    for (par in const_pars){
      kriging_values[pd_back[par,]] = fit[paste("expect_", par, sep="")]
    }

    ####################
    # write output file
    ####################

    filename <- paste(output.dir, "/kriged_values_",stat.ids[i], "_", model_name,".csv" ,sep="")
    write.table(kriging_values, file=filename, col.names=TRUE, sep=",", row.names=FALSE)
    setTxtProgressBar(pb,i)
  }
  close(pb)
}
