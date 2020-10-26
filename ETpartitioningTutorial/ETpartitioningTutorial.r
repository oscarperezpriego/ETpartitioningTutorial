## ET partitioning algorithm to eddy covariance data
## Author: Oscar Perez-Priego
## Propose:This script contains a tutorial of the main functions and steps of the algorithm



# Introduction

# Here we will go through application of the Pérez-Priego et al (2018) partitioning algorithm to eddy covariance data.
#   The script is designed to run on FLUXNET2015 .csv files directly, which ensures consistent variable names, processing and units.
#     The tutorial will use data from the AU-Cum in Australia, but can be applied to any FLUXNET2015 dataset.
#
# Some experience in R will make things easy, but I will try to explain the process step by step so as to be accessible to all backgrounds.
#
# This example has been adapted from the original example:
#
#   https://github.com/oscarperezpriego/ETpartitioning/blob/master/inst/main_ETpartitioning.r

# first things first
#
# The firs step is to import all needed packages:

# library(ETpartitioning) # The package containing the partitioning code

# Here we will source the folder with main functions that I have already customized to save running times
source("D:/scratch/4people/4Belinda/ETpartitioningTutorial/R/cost_function.r")
source("D:/scratch/4people/4Belinda/ETpartitioningTutorial/R/optimize_function.r")
source("D:/scratch/4people/4Belinda/ETpartitioningTutorial/R/calculate_chi_o.r")
source("D:/scratch/4people/4Belinda/ETpartitioningTutorial/R/calculate_WUE_o.r")
source("D:/scratch/4people/4Belinda/ETpartitioningTutorial/R/photos_model.r")
source("D:/scratch/4people/4Belinda/ETpartitioningTutorial/R/transpiration_model.r")
source("D:/scratch/4people/4Belinda/ETpartitioningTutorial/R/gc_model.r")



library(FME)            # Package for parameter optimization

library(bigleaf)        # Some useful tools


## loading the dataset
# Here we open the Rda file. This will load all variables.
load("D:/scratch/4people/4Belinda/ETpartitioningTutorial/data/AU-Cum.HH.2012.2014.Rda")
ds <- Data.F


#-- Estimating long-term effective chi_o (i.e. Ci/Ca) and WUE_o parameters

# The first step is to calculate chi_o via: Ln(chi_o/(1-chi_o)) = 0.0545(Temp-25)-0.58Ln(D)-0.0815+C

#     Temp is the mean seasonal temperature
#     D is the mean seasonal VPD
# This is accomplished using the calculate_chi_o function below.
# The function takes the dataset ds
# as well as designations for the GPP (ColPhotos), VPD (ColVPD), and air temperature (ColTair) columns.

Chi_o <- calculate_chi_o(data= ds
                         ,ColPhotos = "GPP_NT_VUT_MEAN"
                         ,ColVPD = "VPD_F"
                         ,ColTair = "TA_F"
                         ,C = 1.189 ##<< Empirical coeficient for C3 species (see Wang et al., 2017; Plant Nature).
                         ,Z=0.026) ##<< altitude (km)

# Here Z is the site altitude in kilometers
# and C is an empirical coefficient for C3 species (Wang et al., 2017). Here it is assumed as a constant but it varies and is something I am currently working on!!


# Second we calculate WUE_o via WUE_o = Ca Pa (1-Chi_o/1.6 D)
#   Where Ca and Pa are set to 390 ppm and 96 kPa respectively.
#     This is accomplished using the function below

WUE_o <- calculate_WUE_o(data= Data.F
                         ,ColPhotos = "GPP_NT_VUT_MEAN"
                         ,ColVPD = "VPD_F"
                         ,ColTair = "TA_F"
                         ,C = 1.189 ##<< Empirical coeficient for C3 species (see Wang et al., 2017; Plant Nature).
                         ,Z=0.27)


# Prepare for the partitioning

# The partitioning step optimizes parameters for each day using a 5 day moving window.
# Here we convert the timestamp string in the original file into both the individual date and time components,
# giving each day an identifyer to use within the loop (ds$loop).

ds$Chi_o <- Chi_o
ds$WUE_o <- WUE_o
ds$year <- (substr(ds$TIMESTAMP_END, 1, 4))
ds$Month <- as.numeric(substr(ds$TIMESTAMP_END, 5, 6))
ds$DD <- as.numeric(substr(ds$TIMESTAMP_END, 7, 8))
ds$Hour <- (substr(ds$TIMESTAMP_END, 9, 10))
ds$Min <-(substr(ds$TIMESTAMP_END, 11, 12))
ds$rDate <- strptime(paste0(ds$year,"/",ds$Month,"/",ds$DD," ", ds$Hour,":",ds$Min) , format="%Y/%m/%d %H:%M", tz='GMT')
ds$date <- strptime(paste0(paste0(ds$year,"/",ds$Month,"/",ds$DD)) , format="%Y/%m/%d", tz='GMT')
ds$loop <- as.numeric(as.Date(ds$date))-min(as.numeric(as.Date(ds$date)), na.rm=T)+1


#-- Running a subset
#
# The 4 model parameters (a1, Do, Topt and beta, see Perez-Priego et al., 2018) are estimated
#   using a multi-constraint Markov Chain Monte Carlo (MCMC) over a given time period (e.g. 5 days).
#     In the present form, the number of iterations was set to 20000, and the first half of the chains was discarded.
#        We updated the proposed distribution every 500 iterations. Please note that this setup is computationally expensive.
#   Note that the MCMC can be customized as required in the optimal_parameters function to possibly reduce runtimes.
#     In this example we have set the number of iterations to 1000.
#       We will run the model during one month

days_to_run <- subset(ds, year == 2013 & Month == 10)


#-- Partitioning

# For each day in days_to_run we make a temperary dataset (tmp) containing the 5 day window.
#   The optimized model parameters are then estimated (optimal_parameters),
#     which are then passed on the the transpiration model(transpiration_mod).
#       The resulting ET, transpiration, and evaporation are then saved in the estimates dataset.


var_list2 <- unique(days_to_run$loop)

Liston <-  list()

for (i in var_list2){
  print(i)
  # i <- 990
  tmp <- subset(days_to_run,loop %in% c(i-2,i-1, i,i+1, i+2)) ##<< Defining 5 days window in a loop

  tmpp <- subset(tmp, NIGHT == 0)

  #-- optimzazing model parameters

  ##<< we will skip the optimization when important variables are missing
  ## when uncertinties are not provided we assummed them as 10% of the fluxes
  tmpp$NEE_VUT_JOINTUNC <- ifelse(is.na(tmpp$NEE_VUT_USTAR50_JOINTUNC)==TRUE, tmpp$NEE_VUT_MEAN*0.1,  tmpp$NEE_VUT_USTAR50_JOINTUNC)
  tmpp$PPFD_IN <- ifelse(is.na(tmpp$PPFD_IN)==TRUE, tmpp$SW_IN_F,  tmpp$PPFD_IN)
  if (is.na(mean(tmpp$PPFD_IN, na.rm=TRUE))==TRUE) {
    next
  }
  if (is.na(mean(tmpp$CO2_F_MDS, na.rm=TRUE))==TRUE ) {
    next
  }
  if (is.na(mean(tmpp$USTAR, na.rm=TRUE))==TRUE ) {
    next
  }


  ans <-  optimal_parameters(par_lower= c(0,0, 10, 0)
                                  ,par_upper = c(400,0.4, 30, 1)
                                  ,data=tmpp
                                  ,ColPhotos="GPP_NT_VUT_MEAN"
                                  ,ColPhotos_unc="NEE_VUT_JOINTUNC"
                                  ,ColH="H_F_MDS"
                                  ,ColVPD="VPD_F"
                                  ,ColTair="TA_F"
                                  ,ColPair="PA_F"
                                  ,ColQ="PPFD_IN"
                                  ,ColCa="CO2_F_MDS"
                                  ,ColUstar="USTAR"
                                  ,ColWS="WS_F"
                                  ,ColSW_in="SW_IN_F"
                                  ,Chi_o = Chi_o
                                  ,WUE_o= WUE_o)

  par <- as.numeric(ans)

  #-- Estimating transpiration rates
  transpiration_mod <- transpiration_model(
    par=par
    ,data=tmpp
    ,ColPhotos="GPP_NT_VUT_MEAN"
    ,ColH="H_F_MDS"
    ,ColVPD="VPD_F"
    ,ColTair="TA_F"
    ,ColPair="PA_F"
    ,ColQ="PPFD_IN"
    ,ColCa="CO2_F_MDS"
    ,ColUstar="USTAR"
    ,ColWS="WS_F"
    ,ColSW_in="SW_IN_F"
    ,Chi_o = Chi_o
    ,WUE_o= WUE_o)

  #-- Estimating evaporation rates


  lamda <- (3147.5-2.37*(tmpp$TA_F+273.15))*1000 # Latent heat of evaporisation [J kg-1]
  ET <- tmpp$LE_F_MDS/lamda*1000000/18 # from Wm-2 to mmol m-2 s-1
  evaporation_mod <- ET-transpiration_mod

  #-- ET partitioning

  tmpp$ET <- ET
  tmpp$transpiration_mod <- transpiration_mod
  tmpp$evaporation_mod <- evaporation_mod

  #-- Saving parameters

  tmpp$parameter_a1 <- par[1]
  tmpp$parameter_D0 <- par[2]
  tmpp$parameter_T0 <- par[3]
  tmpp$parameter_beta <- par[4]

  #--

  tmpp <- tmpp[tmp$loop == i,] ## selecting the central day
  Liston[[i]] <- tmpp

}

out <- do.call(rbind, Liston)

#-- Converting units and estimating evaporation rates

out$transpiration_mmh <-out$transpiration_mod * (18.01528/1e6)*3600  # from mmol m-2 s-2 to mm per hour
out$ET_mmh <- LE.to.ET(out$LE_F_MDS,out$TA_F) * 3600 # from Wm-2 to mm per hour
out$evaporation_mmh <- out$ET_mmh-out$transpiration_mmh

##-- Plotting diurnal cycles

plot(out$rDate, out$ET_mmh, type="l", ylab=expression(water~flux~(mm~h^-1)), xlab="Time")
lines(out$rDate, out$transpiration_mmh, ylab=expression(water~flux~(mm~h^-1)), xlab="Time", col="green")
lines(out$rDate, out$evaporation_mmh, ylab=expression(water~flux~(mm~h^-1)), xlab="Time", col="red")
legend("topleft", legend = c("ET", "plant transpiration",  "evaporation"), bty = "n",
      lwd=c(1,1,1), cex = 1, col = c( "black","green","red"), pch = c(1, 1,1))

##-- Plotting montly mean diurnal cycle

ET <- aggregate(out$ET_mmh, by=list(out$Hour, out$Month), FUN=mean, na.rm=TRUE)
transpiration <- aggregate(out$transpiration_mmh, by=list(out$Hour, out$Month), FUN=mean, na.rm=TRUE)
evaporation <- aggregate(out$evaporation_mmh, by=list(out$Hour, out$Month), FUN=mean, na.rm=TRUE)
par(oma=c(1,1,1,1),mar=c(5, 5, 5, 4))
plot(ET$Group.1, ET$x, type="l", ylab=expression(water~flux~(mm~h^-1)), xlab="hour")
lines(transpiration$Group.1, transpiration$x, col="green")
lines(evaporation$Group.1, evaporation$x, col="red")
legend("topleft", legend = c("ET", "plant transpiration",  "evaporation"), bty = "n",
       lwd=c(1,1,1), cex = 1, col = c( "black","green","red"), pch = c(1, 1,1))



