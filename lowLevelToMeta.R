####################
#### TITLE:     Simulate null data from low subject level to meta-analysis. Calculate ES, weighted average, CI and coverage based on the simulations.
#### Contents:
####
#### Source Files: HPC - Version
#### First Modified: 15/01/2016
#### Notes:
#################



##
###############
### Notes
###############
##

# Simulating simple null images. Start with blocked design for individual subjects.
# Location and noise varies over subjects (only white, temporal and spatial noise, but magnitude differs).
# Two conditions, contrast is 1 -1.
# These N subjects are pooled using simple OLS pooling.
# The resulting images are converted to Hedges' g and pooled using fixed effects meta-analysis.


########################################
## CI for group maps after low level simulation of fMRI data using neuRosim.
#   * Create 16x16x16 null-images for N subjects and K studies.
#   * Each image is created using the same design.
#   * No between study effects.
#   * Pool each study with ordinary OLS pooling method: T-maps.
#   * Transform each study to ES with formula used in FixRan study.
#   * Aggregate studies using fixed effects meta-analysis.
#   * Construct 3 types of CI in each voxel.
#   * Check coverage in each voxel over all simulations


##
###############
### Preparation
###############
##

# Reset working directory
rm(list=ls())
gc(verbose = FALSE)

# Take argument from master file
input <- commandArgs(TRUE)
  # K'th simulation
  K <- as.numeric(as.character(input)[1]); cat(K, '\n')
  # Which scenario
  SCEN <- as.numeric(as.character(input)[2])

# Set starting seed: it is the product of the amount of voxels, the number of studies and the number of subjects!
starting.seed <- 36864*K
set.seed(starting.seed)

# Set WD
wd <- '/user/scratch/gent/gvo000/gvo00022/vsc40728/Simulation'
setwd(wd)


# Load in libraries
library(AnalyzeFMRI)
library(fmri)
library(lattice)
library(gridExtra)
library(oro.nifti)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(Hmisc)
library(devtools)
library(neuRosim)



# Load in functions from FixRan study: THIS HAS TO COME AFTER ALL LIBRARIES ARE LOADED AS WE SOMETIMES FIX FUNCTIONS THAT ARE BUGGED IN THE PACKAGES
source_url('https://raw.githubusercontent.com/HBossier/FixRanStudyGit/master/Development/functions.R',sha1='c4c3b98288ab8a9bdf0d081f2ace902d5cd13e18')




##
###############
### Simulation steps
###############
##

####************####
#### Global options
####************####
nstud <- 15
nsub <- 15
TR <- 2
nscan <- 200
total <- TR*nscan
on1 <- seq(1,total,40)
on2 <- seq(20,total,40)
onsets <- list(on1,on2)
duration <- list(20,20)
effect.null <- list(0,0)                              ## No effect
effect <- list(1,1) 			                            ## Effect of 1 for designmatrix
DIM <- c(16,16,16)


####************####
#### Scenario specific simulation details
####************####
Noise <- list(
  'S1' = c(1,0,0,0,0,0),
  'S2' = c(1,0,0,0,0,0),
  'S3' = c(0.84,0.05,0.02,0.02,0.02,0.05),
  'S4' = c(0.64,0.15,0.02,0.02,0.02,0.15),
  'S5' = c(0.45,0.25,0.02,0.02,0.02,0.25),
  'S6' = c(0.24,0.35,0.02,0.02,0.02,0.35),
  'S7' = c(0.04,0.45,0.02,0.02,0.02,0.45)
  )



####************####
#### Subject/Study specific simulation details
####************####
# Subject parameters
TrueLoc1 <- c(4,4,4)
TrueLoc2 <- c(10,10,10)
TrueLocations <- rbind(TrueLoc1,TrueLoc2)
TrueWhiteNoise <- Noise[SCEN]
TrueRadius <- 1
COPE <- VARCOPE <- TMAP <- array(NA,dim=c(DIM,nsub))


####************####
#### Study parameters
####************####
SWEIGHTS <- SHEDGE <- SCOPE <- SVARCOPE <- STMAP <- array(NA,dim=c(DIM,nstud))



####************####
#### Design matrices
####************####
# Design Matrices via neuRosim:
#     * We need three design vectors:
#     * The first two have an intercept (needed for analysis).
#        * These will be the two columns (1 -1 contrast) of the design matrix in the analysis.
#     * The third one is used to generate data with a NULL effect.
designC1 <- simprepTemporal(onsets = list(on1), durations = list(duration[[1]]),
                         hrf = "double-gamma", TR = TR, totaltime = total,
                         effectsize = list(effect[[1]]))

designC2 <- simprepTemporal(onsets = list(on2), durations = list(duration[[2]]),
                        hrf = "double-gamma", TR = TR, totaltime = total,
                        effectsize = list(effect[[1]]))

design.null <- simprepTemporal(regions = 2, onsets = onsets, durations = duration,
                         hrf = "double-gamma", TR = TR, totaltime = total,
                         effectsize = effect.null)

# X-matrix in order to fit the model later on (combination of C1 and C2).
x <- fmri.design(matrix(c(simTSfmri(designC1, nscan=nscan, TR=TR, noise="none"),
        simTSfmri(designC2, nscan=nscan, TR=TR, noise="none")),ncol=2),0)



####************####
#### GENERATE DATA
####************####
# For loop over studies
for(t in 1:nstud){
  print(paste('------------------------- STUDY ', t,' -------------------------',sep=''))
  locations1 <- locations2 <- weights <- radius <- c()
  # For loop over nsub
  for(s in 1:nsub){
    print(paste('At study ', t, ', subject ', s,', scenario ',SCEN, ' in simulation ', K, sep=''))
    # First define the locations, weights and radius of the subjects for this study
    LocSeed <- (t*nstud) + (starting.seed + s)
    WeightSeed <- ((nstud*nstud) + (starting.seed + nsub)) + ((t-1)*nstud) + (s)
      # Location
      locations <- subjectfMRIlocation(TrueLocations,seed=LocSeed,DIM=DIM)
      locations1 <- rbind(locations1,locations[1,])
      locations2 <- rbind(locations2,locations[2,])
        coordinates <- list(c(locations1[s,]),c(locations2[s,]))
      # Weights
      weights <- rbind(weights,subjectfMRINoise(Noise[[SCEN]],seed = WeightSeed))
      # Radius
      radius.tmp <- TrueRadius + sample(c(0,1),size=2)
      radius <- rbind(radius,'rad' = radius.tmp)

    # Define two regions (which does nothing as there is no effect, )
    regions <- simprepSpatial(regions = 2, coord = coordinates, radius = list(radius[s,1],radius[s,2]), form ="cube", fading = 0)
      rm(coordinates)

    # Weighting structure: white, temporal and spatial noise.
    #   * Order = white, temporal, low-frequency, physyiological, task related and spatial.
    w <- round(c(weights[s,]),2)

    # Base value
    base <- 5

    # Actual simulated data
    sim.data <- simVOLfmri(design=design.null, image=regions, base=base, dim=DIM, SNR=0.5,
                 type ="gaussian", noise= "mixture", spat="gaussRF", FWHM=2, weights=w, verbose = TRUE)
        rm(w)

    # Not smoothing in scenario 1
    if(SCEN!=1){
      # 3D Gaussian Kernel over the 4D data
      fwhm <- 3
      sigma <- fwhm/(sqrt(8)*log(2))
      smoothint <- -round(2*sigma):round(2*sigma)

      for(i in 1:nscan){
       sim.data[,,,i] <- GaussSmoothArray(sim.data[,,,i], voxdim=c(1,1,1), ksize = length(smoothint), sigma = diag(sigma, 3))
      }
    }

    # Reshape the data for fMRI analysis and make it the correct class
    datafmri <- list(ttt=writeBin(c(sim.data), raw(),4), mask=array(1,dim=DIM), dim = c(DIM, nscan))
    class(datafmri) <- "fmridata"
      rm(sim.data)
    ####************####
    #### ANALYZE DATA
    ####************####


    # Fitting GLM model: estimated AR(1)-coefficients are used to whiten data, may produce warnings because data is pre-smoothed.
    model <- fmri.lm(datafmri,x, actype = "accalc", keep="all",contrast=c(1,-1))

    # Estimated contrast of parameter beta's from model
    COPE.sub <- model$cbeta
      COPE[,,,s] <- COPE.sub
    VARCOPE.sub <- model$var
      VARCOPE[,,,s] <- VARCOPE.sub
    # Constructing t-map
    TMAP.sub <- array(c(COPE.sub)/sqrt(c(VARCOPE.sub)), dim=c(DIM))
      TMAP[,,,s] <- TMAP.sub

      rm(COPE.sub,VARCOPE.sub,TMAP.sub)

  }

  ####************####
  #### GROUP ANALYSIS
  ####************####

  # Group COPE (average)
  GCOPE <- apply(COPE,c(1,2,3),mean)

  # Now we will do the OLS estimation of the variance
  GVARCOPE <- apply(COPE,c(1,2,3),var)

  # TMAP
  GTMAP <- GCOPE/sqrt(GVARCOPE/(nsub))


  ####************####
  #### TRANSFORM TO ES
  ####************####
  # Transform to an ES using hedgeG function, for each study
  HedgeG <- apply(matrix(GTMAP,ncol=1),1,FUN=hedgeG,N=nsub)
  # Calculate variance of ES
  VarianceHedgeG <- apply(matrix(HedgeG,ncol=1),1,FUN=varHedge,N=nsub)
  # Weights of this study
  weigFix <- 1/VarianceHedgeG

    # Put GCOPE, GVARCOPE, GTMAP, hedge's G and weights in vector
    SCOPE[,,,t] <- GCOPE
    SVARCOPE[,,,t] <- GVARCOPE
    STMAP[,,,t] <- GTMAP
    SHEDGE[,,,t] <- HedgeG
    SWEIGHTS[,,,t] <- weigFix

    rm(GCOPE,GVARCOPE,HedgeG,weigFix)
}


# Re-format arrays
SHEDGE.mat <- matrix(SHEDGE,ncol=nstud)
SWEIGHTS.mat <- matrix(SWEIGHTS,ncol=nstud)

# Now calculate weighted average.
WeightedAvg <- (apply((SHEDGE.mat*SWEIGHTS.mat),1,sum))/(apply(SWEIGHTS.mat,1,sum))

# Calculate variance of weighted average
varWeightAvg <- 1/apply(SWEIGHTS.mat,1,sum)

####************####
#### CALCULATE CI'S
####************####

# CI for weighted average based on normal distribution
CI.upper.norm <- matrix(WeightedAvg,ncol=1) + (1.96 * sqrt(matrix(varWeightAvg,ncol=1)))
CI.lower.norm <- matrix(WeightedAvg,ncol=1) - (1.96 * sqrt(matrix(varWeightAvg,ncol=1)))

# CI for weighted average based on t-distribution
CI.upper.t <- matrix(WeightedAvg,ncol=1) + (qt(0.975,df=nstud-1) * sqrt(matrix(varWeightAvg,ncol=1)))
CI.lower.t <- matrix(WeightedAvg,ncol=1) - (qt(0.975,df=nstud-1) * sqrt(matrix(varWeightAvg,ncol=1)))

# CI for weighted average based on variance CI
CI.weightedVariance <- (apply((SWEIGHTS.mat*(SHEDGE.mat - WeightedAvg)^2),c(1),sum))/((nstud - 1) * apply(SWEIGHTS.mat,1,sum))
CI.upper.weightVar <- matrix(WeightedAvg,ncol=1) + (qt(0.975,df=nstud-1) * sqrt(matrix(CI.weightedVariance,ncol=1)))
CI.lower.weightVar <- matrix(WeightedAvg,ncol=1) - (qt(0.975,df=nstud-1) * sqrt(matrix(CI.weightedVariance,ncol=1)))



##
###############
### Save objects
###############
##

save(SCOPE, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/SCOPE_',K,sep=''))
save(SVARCOPE, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/SVARCOPE_',K,sep=''))
save(STMAP, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/STMAP_',K,sep=''))
save(SHEDGE, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/SHEDGE_',K,sep=''))
save(WeightedAvg, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/WeightedAvg_',K,sep=''))
save(varWeightAvg, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/varWeightAvg_',K,sep=''))
save(CI.upper.norm, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/CI.upper.norm_',K,sep=''))
save(CI.lower.norm, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/CI.lower.norm_',K,sep=''))
save(CI.upper.t, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/CI.upper.t_',K,sep=''))
save(CI.lower.t, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/CI.lower.t_',K,sep=''))
save(CI.upper.weightVar, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/CI.upper.weightVar_',K,sep=''))
save(CI.lower.weightVar, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/CI.lower.weightVar_',K,sep=''))







