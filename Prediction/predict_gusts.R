
####################################################################################################
# Author: Philipp Ertz, July 2024. Last edited: May 2025
#
# This is an R script to perform the wind gust prediction process for all different gust models.
# The models are selected via command line arguments from the shell.
####################################################################################################
library(dplyr)
library(tidyr)
library(evd)

# pass model name as command line argument
args <- commandArgs(trailingOnly=TRUE)
models <- args
store.q_thr <- FALSE # save quantiles and threshold predictions instead of complete samples
dir.path <- "./"
eligible_models <- c("Baseline0",
                      "Baseline_mu2",
                      "Baseline_vmean",
                      "Baseline_optimal",
                      "SM_mu0_f",
                      "SM_mu0_mu1_f",
                      "SM_mu0_mu2_f",
                      "SM_mu0_sigma0_f",
                      "SM_mu0_mu1_mu2_f",
                      "LocMod")

for (model_name in models){
  if (!(model_name %in% eligible_models)){
      print(paste(model_name,"is not a valid model."))
      next
  }

  # get station IDs, coordinates and elevations
  stat_info <- read.table(file = paste(dir.path, "Data/used_stations.csv", sep=""), sep=",", header=TRUE, colClasses= c(rep("character",3), rep("numeric",3), rep("character",2)))
  stat.ids <- str_pad(stat_info$Stations_id, width = 5, side = "left", pad = "0")
  z_station <- stat_info$Stationshoehe
  z_grid <- unname(unlist(read.table(paste(dir.path,"Data/COSMO-REA6_HSURF_stations.csv",sep=""), header=TRUE, sep=",")$HSURF))

  # output directory
  output.dir <- paste(dir.path, "Prediction/", model_name, sep="")
  if (!dir.exists(output.dir)){
    dir.create(output.dir, recursive=TRUE)
  }
  print(paste("Predicting from ", model_name, ". Prediction files will be saved in ", output.dir, sep=""))

  # read data for fx, VMEAN and VMAX
  fx <- read.table(paste(dir.path,"/Data/fx_evaluation.csv",sep=""), header=TRUE, sep=",")
  vmax <- read.table(paste(dir.path,"/Data/vmax_evaluation.csv",sep=""), header=TRUE, sep=",")
  vmean <- read.table(paste(dir.path,"/Data/vmean_evaluation.csv",sep=""), header=TRUE, sep=",")

  M <- length(stat.ids)
  dates <- fx[,1]
  fx <- fx[,2:(M+1)]
  vmax <- vmax[,2:(M+1)]
  vmean <- vmean[,2:(M+1)]
  rownames(fx) <- dates
  rownames(vmax) <- dates
  rownames(vmean) <- dates

  # get elevation difference between REA6 and weather stations
  elevation <- z_station - z_grid

  normalize_predictor <- function(x, stat_id, stat_list){
    # select all stations but the current one
    ind <- which(stat_list==stat_id)[1]
    range <- 1:length(stat_list)
    stat_index <- range[range!=ind]

    # normalize predictor
    mu <- mean(unlist(x[,stat_index]))
    sd <- sd(unlist(x[,stat_index]))

    return((x[,ind] - mu)/sd)
  }

  print("Finished reading observations. Start prediction of wind gusts.")
  ######################################
  # perform actual prediction process
  ######################################
  
  n.sample = 10000

  # create progress bar
  pb = txtProgressBar(min=0, max=M, initial=1)

  # loop over all stations
  for (s in 1:length(stat.ids)){
    # get standardized predictor vectors
    if (model_name=="LocMod"){
        # individual station normalization
        vmax_pred <- (vmax[,s] - mean(vmax[,s]))/sd(vmax[,s])
        vmean_pred <- (vmean[,s] - mean(vmean[,s]))/sd(vmean[,s])
        N <- length(vmax_pred)
    } else {
        # whole data set normalization
        vmax_pred <- normalize_predictor(vmax, stat.ids[s], stat.ids)
        vmean_pred <- normalize_predictor(vmean, stat.ids[s], stat.ids)
        N <- length(vmax_pred)
    }

    # get vmean_offset values
    vmean_loc <- unlist(vmean[,s])

    # get parameter samples
    sampling_ind <- sample(1:1000,n.sample, replace=TRUE) # generate a vector sampling from the Markov chains
    if (substring(model_name,1,1)=="B"){
        MC <- read.table(paste(dir.path,paste("/Model_fits/",model_name,"/",model_name,"_",stat.ids[s],"_fit.csv", sep=""),sep=""), sep=",", header=TRUE)[sampling_ind,]
    } else if(model_name=="LocMod"){
        MC <- read.table(paste(dir.path,paste("/Model_fits/",model_name,"/",model_name,"_",stat.ids[s],"_fit.csv", sep=""),sep=""), sep=",", header=TRUE)[sampling_ind,]
    } else if (substring(model_name,1,1)=="S"){
        MC <- read.table(paste(dir.path,paste("Kriging/",model_name,"/kriged_values_",stat.ids[s],"_",model_name,".csv", sep=""),sep=""), sep=",", header=TRUE)[sampling_ind,]
    } else {stop("Invalid model")}

    rownames(MC) <- 1:10000
    # loop over all dates
    fx_pred <- matrix(0,nrow=n.sample,ncol=N)
    for (n in 1:N){
      # this is the absolute baseline model
      if (model_name %in% c("Baseline0",
                          "SM_mu0_f",
                          "SM_mu0_mu1_f",
                          "SM_mu0_mu2_f",
                          "SM_mu0_mu1_mu2_f"
                          "SM_mu0_sigma0_f")){
        mu <- MC$mu0 + vmax_pred[n] * MC$mu1
        sigma_pred <- MC$sigma0 + vmax_pred[n] * MC$sigma1
      } else {
        # ConstMod uses a differnet naming convention for the parameters.
        mu <- MC$mu.1 + vmax_pred[n] * MC$mu.2
        sigma_pred <- MC$sigma.1 + vmax_pred[n] * MC$sigma.2
      }
      
      # increase mu with altitude predictor for ConstMod versions without VMEAN-predictor
      if (model_name %in% c("Baseline_mu2",
                            "Baseline_vmean")){
        mu <- mu + MC$mu.3 * elevation[s]/100
      }      

      # include VMEAN as predictor
      if (model_name %in% c("Baseline_optimal")){
        mu <- mu + vmean_pred[n] * MC$mu.3 + MC$mu.4 * elevation[s]/100
        sigma_pred <- sigma_pred + vmean_pred[n] * MC$sigma.3
      }

      # only include VMEAN for LocMod
      if (model_name == "LocMod"){
        mu <- mu + vmean_pred[n] * MC$mu.3
        sigma_pred <- sigma_pred + vmean_pred[n] * MC$sigma.3
      }

      # include VMEAN and altitude for the spatial models
      if (model_name %in% c("SM_mu0_f",
                            "SM_mu0_mu1_f",
                            "SM_mu0_mu2_f",
                            "SM_mu0_mu1_mu2_f",
                            "SM_mu0_sigma0_f")){
        mu <- mu + vmean_pred[n] * MC$mu2 + MC$mu3 * elevation[s]/100
        sigma_pred <- sigma_pred + vmean_pred[n] * MC$sigma2
      }

      # sampple from Gumbel distribution
      fx_pred[,n] <- rgumbel(n=n.sample,loc=mu, scale=exp(sigma_pred))

      # increment fx_pred by VMEAN if necessary
      if (model_name %in% c("Baseline_vmean",
                            "Baseline_optimal",
                            "LocMod",
                            "SM_mu0_f",
                            "SM_mu0_mu1_f",
                            "SM_mu0_mu2_f",
                            "SM_mu0_mu1_mu2_f",
                            "SM_mu0_sigma0_f")){
        fx_pred[,n] <- fx_pred[,n] + vmean_loc[n]
      }
    }
    
    # store data in data.frame
    fx_pred <- data.frame(fx_pred)
    colnames(fx_pred) <- dates

    if (store.q_thr==FALSE){
      # write output file
      write.table(fx_pred,file=paste(output.dir, "/fx_predicted_",stat.ids[s], ".csv", sep=""), sep=",", row.names=T, col.names=T)
    }

    # store quantiles and threshold predictions **instead**
    if (store.q_thr==TRUE){
      thr <- c(14,18,25,29,33)
      qs <- c(0.0,0.001, seq(0.01,0.99,0.01), 0.999,1)
      
      # calculate quantiles for saving
      q_pred <- data.frame(apply(fx_pred, c(2), quantile, probs=qs,type=8))
      colnames(q_pred) <- dates
  
      # calculate threshold excesses
      thr_pred = matrix(0,nrow=length(thr), ncol=N)
      for (t in 1:length(thr)){
        for (d in 1:length(dates)){
          thr_pred[t,d] = mean(fx_pred[,d]>thr[t])
        }
      }
      thr_pred <- data.frame(thr_pred)
      colnames(thr_pred) <- dates
      rownames(thr_pred) <- thr

      # write output files for quantiles and thresholds
      write.table(q_pred,file=paste(output.dir,"/q_predicted_",stat.ids[s], ".csv", sep=""), sep=",", row.names=T, col.names=T)
      write.table(thr_pred,file=paste(output.dir, "/thr_predicted_",stat.ids[s], ".csv", sep=""), sep=",", row.names=T, col.names=T)
    }
    setTxtProgressBar(pb,s)
  }
  close(pb)
}
    
    
          



