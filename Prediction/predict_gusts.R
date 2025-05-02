
####################################################################################################
# Author: Philipp Ertz, July 2024
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
hub_dir <- "/automount/hubhome/s6phertz/Coming_decade/Gust_paper/"
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
                        "Baseline0",
                        "Baseline_mu2_mount",
                        "Baseline_optimal",
                        "Baseline_vmean_mount",
                        "Baseline_maxime", # includes suggestion of altitude handling by Maxime Taillardat, Meteo-France
                        "SM_mu0_f",
                        "SM_mu0_mu1_f",
                        "SM_mu0_mu2_f",
                        "SM_mu0_sigma0_f",
                        "SM_mu0_mu1_mu2_f",
                        "SM_mu0_mu2_sigma0_f",
                        "SM_mu0_mu1_mu2_mount",
                        "SM_mu0_mt",
                        "LocMod")

for (model_name in models){
    if (!(model_name %in% eligible_models)){
        print(paste(model_name,"is not a valid model."))
        next
    }

    labels <- c(unlist(unname(read.table("/automount/hubhome/s6phertz/Coming_decade/Gust_paper/data/used_station_ids.csv", skip=1, colClasses = c("character")))))

    # check directory
    dir_path <- paste(hub_dir, "Predictions/", model_name, sep="")
    if (!dir.exists(dir_path)){
        dir.create(dir_path, recursive=TRUE)
    }
    print(paste("Predicting from ", model_name, ". Prediction files will be saved in ", hub_dir, "Predictions/", model_name, sep=""))

    # exclude mountain top stations. There is a clear mapping from the model names
    if (model_name %in% c("Baseline_mu2",
                        "Baseline_vmean",
                        "Baseline_vmean_pred",
                        "Baseline_flat",
                        "SM0",
                        "SM_mu0",
                        "SM_mu0_sigma0",
                        "SM_mu0_mu1",
                        "SM_mu0_sigma1",
                        "SM_mu0_mu2",
                        "SM_mu0_sigma2",
                        "SM_mu0_mu1_sigma0",
                        "SM_mu0_mu1_mu2")){
        flat_flag <- TRUE
    } else if(model_name %in% c("Baseline0",
                                "Baseline_mu2_mount",
                                "Baseline_optimal",
                                "Baseline_vmean_mount",
                                "Baseline_maxime",
                                "SM_mu0_f",
                                "SM_mu0_mt",
                                "SM_mu0_mount",
                                "SM_mu0_mu1_f",
                                "SM_mu0_mu1_mount",
                                "SM_mu0_mu2_f",
                                "SM_mu0_mu2_mount",
                                "SM_mu0_sigma0_f",
                                "SM_mu0_sigma0_mount",
                                "SM_mu0_mu1_mu2_f",
                                "SM_mu0_mu2_sigma0_f",
                                "SM_mu0_mu1_mu2_mount",
                                "LocMod")){
        flat_flag <- FALSE
    } else {
        stop("invalid model")
    }

    # read data for fx, VMEAN and VMAX
    fx <- read.table(paste(hub_dir,"data/fx_for_evaluation.csv",sep=""), header=TRUE, sep=",")
    vmax <- read.table(paste(hub_dir,"data/vmax_for_evaluation.csv",sep=""), header=TRUE, sep=",")
    vmean <- read.table(paste(hub_dir,"data/vmean_for_evaluation.csv",sep=""), header=TRUE, sep=",")

    M <- length(labels)
    dates <- fx[,1]
    fx <- fx[,2:(M+1)]
    vmax <- vmax[,2:(M+1)]
    vmean <- vmean[,2:(M+1)]
    rownames(fx) <- dates
    rownames(vmax) <- dates
    rownames(vmean) <- dates

    # get station IDs, coordinates and elevations
    stat_info <- read.table(paste(hub_dir,"data/used_stations_109.csv", sep=""), sep=",", header=TRUE)
    z_station <- stat_info$Stationshoehe
    z_grid <- unname(unlist(read.table(paste(hub_dir,"data/rea6_elevation.txt",sep=""))))

    # get elevation difference between REA6 and weather stations
    elevation <- z_station - z_grid # use height difference between REA6-Grid and Synop-Station as further predictor
    if (model_name %in% c("Baseline_maxime", "SM_mu0_mt")){elevation <- elevation * (elevation>0)}

    # exclude mountain top stations where necessary
    if (flat_flag==TRUE){
        fit_ind = which(z_station<800)
        stat_info = stat_info[fit_ind,]
        labels = labels[fit_ind]

        # Read wind gust data
        fx = fx[fit_ind]
        vmax = vmax[fit_ind]
        vmean = vmean[fit_ind]

        # make a height predictor vector
        elevation = elevation[fit_ind] # use height difference between REA6-Grid and Synop-Station as further predictor for the location parameter
        M <- length(labels)
    }

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

    # perform actual prediction process
    # specify sample size
    size = 10000

    # thresholds and quantiles to save:
    thr <- c(14,18,25,29,33)
    qs <- c(0.0,0.001, seq(0.01,0.99,0.01), 0.999,1)
    range <- 1:35#length(labels)

    # create progress bar
    pb = txtProgressBar(min=0, max=M, initial=1)

    # loop over all stations
    for (s in range){
        # get standardized predictor vectors
        if (model_name=="LocMod"){
            vmax_pred <- (vmax[,s] - mean(vmax[,s]))/sd(vmax[,s])
            vmean_pred <- (vmean[,s] - mean(vmean[,s]))/sd(vmean[,s])
            N <- length(vmax_pred)
        } else {
            vmax_pred <- normalize_predictor(vmax, labels[s], labels)
            vmean_pred <- normalize_predictor(vmean, labels[s], labels)
            N <- length(vmax_pred)
        }

        # get vmean_offset values
        vmean_loc <- unlist(vmean[,s])

        # get parameter samples
        sampling_ind <- sample(1:1000,size, replace=TRUE) # generate a vector sampling from the Markov chains
        if (substring(model_name,1,1)=="B"){
            MC <- read.table(paste(hub_dir,paste("Model_fits/",model_name,"/",tolower(model_name),"_fit_",labels[s],".csv", sep=""),sep=""), sep=",", header=TRUE)[sampling_ind,]
        } else if(model_name=="LocMod"){
            MC <- read.table(paste(hub_dir,paste("Model_fits/",model_name,"/",model_name,"_fit_",labels[s],".csv", sep=""),sep=""), sep=",", header=TRUE)[sampling_ind,]
        } else if (substring(model_name,1,1)=="S"){
            MC <- read.table(paste(hub_dir,paste("Kriging/kriged_values_",labels[s],"_",model_name,".csv", sep=""),sep=""), sep=",", header=TRUE)[sampling_ind,]
        } else {stop("Invalid model")}

        rownames(MC) <- 1:10000
        # loop over all dates
        fx_pred <- matrix(0,nrow=size,ncol=N)
        for (n in 1:N){
            # this is the absolute baseline model
            if (model_name %in% c("Baseline0",
                                "SM_mu0",
                                "SM_mu0_mu1",
                                "SM_mu0_mu2",
                                "SM_mu0_sigma0",
                                "SM_mu0_sigma1",
                                "SM_mu0_sigma2",
                                "SM_mu0_mu1_mu2",
                                "SM_mu0_mu1_sigma0",
                                "SM_mu0_mount",
                                "SM_mu0_mu1_mount",
                                "SM_mu0_mu2_mount",
                                "SM_mu0_sigma0_mount",
                                "SM_mu0_mu1_mu2_mount",
                                "SM_mu0_f",
                                "SM_mu0_mu1_f",
                                "SM_mu0_mu2_f",
                                "SM_mu0_mu1_mu2_f",
                                "SM_mu0_mu2_sigma0_f",
                                "SM_mu0_sigma0_f",
                                "SM_mu0_mt")){
                mu <- MC$mu0 + vmax_pred[n] * MC$mu1
                sigma_pred <- MC$sigma0 + vmax_pred[n] * MC$sigma1
            } else {
                # later models use a different naming convention for the parameters, as I started using the adaptive code.
                mu <- MC$mu.1 + vmax_pred[n] * MC$mu.2
                sigma_pred <- MC$sigma.1 + vmax_pred[n] * MC$sigma.2
            }
            # increase mu with altitude predictor (models with different numbering of coefficients)
            if (model_name %in% c("Baseline_mu2",
                                "Baseline_mu2_mount",
                                "Baseline_vmean",
                                "Baseline_vmean_mount")){

                mu <- mu + MC$mu.3 * elevation[s]/100
            }

            # include VMEAN as predictor
            if (model_name %in% c("Baseline_vmean_pred",
                                "Baseline_optimal",
                                "Baseline_maxime")){
                mu <- mu + vmean_pred[n] * MC$mu.3 + MC$mu.4 * elevation[s]/100
                sigma_pred <- sigma_pred + vmean_pred[n] * MC$sigma.3
            }

            # only include VMEAN for LocMod
            if (model_name == "LocMod"){
                mu <- mu + vmean_pred[n] * MC$mu.3
                sigma_pred <- sigma_pred + vmean_pred[n] * MC$sigma.3
            }

            # include VMEAN and altitude for the spatial models
            if (model_name %in% c("SM_mu0",
                                "SM_mu0_mu1",
                                "SM_mu0_mu2",
                                "SM_mu0_sigma0",
                                "SM_mu0_sigma1",
                                "SM_mu0_sigma2",
                                "SM_mu0_mu1_mu2",
                                "SM_mu0_mu1_sigma0",
                                "SM_mu0_mount",
                                "SM_mu0_mu1_mount",
                                "SM_mu0_mu2_mount",
                                "SM_mu0_sigma0_mount",
                                "SM_mu0_mu1_mu2_mount",
                                "SM_mu0_f",
                                "SM_mu0_mu1_f",
                                "SM_mu0_mu2_f",
                                "SM_mu0_mu1_mu2_f",
                                "SM_mu0_sigma0_f",
                                "SM_mu0_mu2_sigma0_f",
                                "SM_mu0_mt")){
                mu <- mu + vmean_pred[n] * MC$mu2 + MC$mu3 * elevation[s]/100
                sigma_pred <- sigma_pred + vmean_pred[n] * MC$sigma2
            }

            # sampple from Gumbel distribution
            fx_pred[,n] <- rgumbel(n=size,loc=mu, scale=exp(sigma_pred))

            # increment fx_pred by VMEAN if necessary
            if (model_name %in% c("Baseline_vmean",
                                "Baseline_vmean_pred",
                                "Baseline_vmean_mount",
                                "Baseline_optimal",
                                "Baseline_maxime",
                                "LocMod",
                                "SM_mu0",
                                "SM_mu0_mu1",
                                "SM_mu0_mu2",
                                "SM_mu0_sigma0",
                                "SM_mu0_sigma1",
                                "SM_mu0_sigma2",
                                "SM_mu0_mu1_mu2",
                                "SM_mu0_mu1_sigma0",
                                "SM_mu0_mount",
                                "SM_mu0_mu1_mount",
                                "SM_mu0_mu2_mount",
                                "SM_mu0_sigma0_mount",
                                "SM_mu0_mu1_mu2_mount",
                                "SM_mu0_f",
                                "SM_mu0_mu1_f",
                                "SM_mu0_mu2_f",
                                "SM_mu0_mu1_mu2_f",
                                "SM_mu0_sigma0_f",
                                "SM_mu0_mu2_sigma0_f",
                                "SM_mu0_mt")){
                fx_pred[,n] <- fx_pred[,n] + vmean_loc[n]
            }
        }
        # store data to data.frame
        fx_pred <- data.frame(fx_pred)
        colnames(fx_pred) <- dates

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

        # write data to file
        write.table(q_pred,file=paste(hub_dir,"Predictions/",model_name, "/q_predicted_",labels[s], ".csv", sep=""), sep=",", row.names=T, col.names=T)
        write.table(thr_pred,file=paste(hub_dir,"Predictions/",model_name, "/thr_predicted_",labels[s], ".csv", sep=""), sep=",", row.names=T, col.names=T)
        setTxtProgressBar(pb,s)
    }
    close(pb)
}



