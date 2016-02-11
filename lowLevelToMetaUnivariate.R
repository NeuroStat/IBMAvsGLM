####################
#### TITLE:     Simulate null data from low subject level to meta-analysis for one voxel only (UNIVARIATE APPROACH).
#### Contents:
####
#### Source Files: HPC - Version
#### First Modified: 03/02/2016
#### Notes:
#################



##
###############
### Notes
###############
##

# Here, we try to see if the coverages in the univariate approach are indeed what we suspect them to be.



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
DIM <- c(1,1,1)


####************####
#### Scenario specific simulation details
####************####
Noise <- list(
  'S1' = c(1,0,0,0,0,0)
  )



####************####
#### Subject/Study specific simulation details
####************####
# Subject parameters
TrueLocation <- c(1,1,1)
TrueWhiteNoise <- Noise[SCEN]
TrueRadius <- 1
COPE <- array(NA,dim=nsub)


####************####
#### Study parameters
####************####
SWEIGHTS <- SHEDGE <- SCOPE <- SVARCOPE <- STMAP <- array(NA,dim=nstud)



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

design.null <- simprepTemporal(regions = 1, onsets = onsets, durations = duration,
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
  # For loop over nsub
  for(s in 1:nsub){
    print(paste('At study ', t, ', subject ', s,', scenario ',SCEN, ' in simulation ', K, sep=''))
    # First define the locations, weights and radius of the subjects for this study
      # Location
      coordinates <- list(TrueLocation)
      # Radius
      radius <- TrueRadius

    # Define region
    regions <- simprepSpatial(regions = 1, coord = coordinates, radius = list(radius), form ="cube", fading = 0)
      rm(coordinates)

    # Weighting structure: white, temporal and spatial noise.
    #   * Order = white, temporal, low-frequency, physyiological, task related and spatial.
    w <- c(1,0,0,0,0,0)

    # Base value
    base <- 5

    # Actual simulated data
    sim.data <- simVOLfmri(design=design.null, image=regions, base=base, dim=DIM, SNR=0.5,
                 type ="gaussian", noise= "mixture", spat="gaussRF", FWHM=2, weights=w, verbose = TRUE)
        rm(w)


  	####************####
    #### ANALYZE DATA
    ####************####
	LM.sim.data <- array(sim.data,dim=nscan)
	fit <- lm(LM.sim.data~x[,-3])
		b1 <- summary(fit)$coefficients[2]
		b2 <- summary(fit)$coefficients[3]
	BETAS <- c(b1,b2)
	CONTRAST <- c(1,-1)
    # Estimated contrast of parameter beta's
	COPE.sub <- CONTRAST %*% BETAS
      COPE[s] <- COPE.sub

      rm(COPE.sub,sim.data,LM.sim.data)
  }

  ####************####
  #### GROUP ANALYSIS
  ####************####

  # Group COPE (average)
  GCOPE <- mean(COPE,na.rm=TRUE)

  # Now we will do the OLS estimation of the variance
  GVARCOPE <- var(COPE,na.rm=TRUE)

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
    SCOPE[t] <- GCOPE
    SVARCOPE[t] <- GVARCOPE
    STMAP[t] <- GTMAP
    SHEDGE[t] <- HedgeG
    SWEIGHTS[t] <- weigFix

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

# CI for weighted average based on weighted variance CI
CI.weightedVariance <- (apply((SWEIGHTS.mat*(SHEDGE.mat - WeightedAvg)^2),c(1),sum))/((nstud - 1) * apply(SWEIGHTS.mat,1,sum))
CI.upper.weightVar <- matrix(WeightedAvg,ncol=1) + (qt(0.975,df=nstud-1) * sqrt(matrix(CI.weightedVariance,ncol=1)))
CI.lower.weightVar <- matrix(WeightedAvg,ncol=1) - (qt(0.975,df=nstud-1) * sqrt(matrix(CI.weightedVariance,ncol=1)))




##
###############
### Save objects
###############
##

save(SCOPE, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/UNI_SCOPE_',K,sep=''))
save(SVARCOPE, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/UNI_SVARCOPE_',K,sep=''))
save(STMAP, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/UNI_STMAP_',K,sep=''))
save(SHEDGE, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/UNI_SHEDGE_',K,sep=''))
save(WeightedAvg, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/UNI_WeightedAvg_',K,sep=''))
save(varWeightAvg, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/UNI_varWeightAvg_',K,sep=''))
save(CI.upper.norm, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/UNI_CI.upper.norm_',K,sep=''))
save(CI.lower.norm, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/UNI_CI.lower.norm_',K,sep=''))
save(CI.upper.t, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/UNI_CI.upper.t_',K,sep=''))
save(CI.lower.t, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/UNI_CI.lower.t_',K,sep=''))
save(CI.upper.weightVar, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/UNI_CI.upper.weightVar_',K,sep=''))
save(CI.lower.weightVar, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/UNI_CI.lower.weightVar_',K,sep=''))





