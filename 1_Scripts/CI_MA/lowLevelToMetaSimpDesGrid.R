####################
#### TITLE:     Simulate null data from low subject level to meta-analysis. Simple design, 4096 voxels.
#### Contents:
####
#### Source Files: HPC - Version
#### First Modified: 16/02/2016
#### Notes:
#################



##
###############
### Notes
###############
##


# We let the sample size and amount of studies increase.
# The amount of studies in the meta-analysis are 2, 5 or 10. The sample size goes from 10 - 100.
# The grid of voxels is 16 x 16 x 16



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
  K <- try(as.numeric(as.character(input)[1]),silent=TRUE)
  # Which scenario
  SCEN <- try(as.numeric(as.character(input)[2]),silent=TRUE)
  # Which machine
  MACHINE <- try(as.character(input)[3],silent=TRUE)
    # If no machine is specified, then it has to be this machine!
    if(is.na(MACHINE)){
      MACHINE <- 'MAC'
      K <- 1
      SCEN <- 10
    }

# Set starting seed: it is the product of the amount of voxels, the number of studies and the number of subjects!
starting.seed <- 36865*K
set.seed(starting.seed)

# Set WD
if(MACHINE=='HPC'){
  wd <- '/user/scratch/gent/gvo000/gvo00022/vsc40728/Simulation'
}
if(MACHINE=='MAC'){
  wd <- '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/Take6'
}
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
if(MACHINE == 'MAC'){
  source('~/Dropbox/PhD/PhDWork/Meta\ Analysis/R\ Code/Studie_FixRan/FixRanStudyGit.git/Development/functions.R')
}
if(MACHINE == 'HPC'){
  source('/user/scratch/gent/gvo000/gvo00022/vsc40728/Simulation/functions.R')
}



##
###############
### Simulation steps
###############
##

####************####
#### Global options
####************####
TR <- 2
nscan <- 200
total <- TR*nscan
on1 <- seq(1,total,40)
onsets <- list(on1)
duration <- list(20)
effect.null <- list(0)                              ## No effect
effect <- list(1) 			                            ## Effect of 1 for designmatrix
DIM <- c(16,16,16)


####************####
#### Scenario specific simulation details
####************####


# First make the data frame with the combinations of subjects and studies
OverView <- data.frame('Subjects' = rep(seq(10,100,by=10), 3),
                      'Studies' = rep(c(2,5,10), each= 10))

# Now take the correct amount of subjects and studies (according to SCEN)
nsub <- OverView[SCEN,'Subjects']
nstud <- OverView[SCEN,'Studies']


# At the moment, only consider one noise structure:
	Noise <- list(
	  'S1' = c(1,0,0,0,0,0)
	  )

####************####
#### Subject/Study specific simulation details
####************####
# Subject parameters
TrueLocations <- c(4,4,4)
TrueWhiteNoise <- Noise[1]						# MIND THE INDEX HERE!
TrueRadius <- 1
COPE <- VARCOPE <- TMAP <- array(NA,dim=c(prod(DIM),nsub))


####************####
#### Study parameters
####************####
SWEIGHTS <- SHEDGE <- SCOPE <- SVARCOPE <- STMAP <- array(NA,dim=c(prod(DIM),nstud))


####************####
#### Design matrices
####************####
# Design Matrices via neuRosim:
#     * We need two design vectors:
#     * The first one have an intercept (needed for analysis).
#        * This will be the column of the design matrix in the analysis.
#     * The second one is used to generate data with a NULL effect.
design.Cond1 <- simprepTemporal(onsets = list(on1), durations = list(duration[[1]]),
                       hrf = "double-gamma", TR = TR, totaltime = total,
                       effectsize = list(effect[[1]]))

design.null <- simprepTemporal(regions = 1, onsets = onsets, durations = duration,
                       hrf = "double-gamma", TR = TR, totaltime = total,
                       effectsize = effect.null)

# X-matrix in order to fit the model later on.
x <- matrix(c(simTSfmri(design.Cond1, nscan=nscan, TR=TR, noise="none")),ncol=1)



####************####
#### GENERATE DATA
####************####
# For loop over studies
for(t in 1:nstud){
  print(paste('At study ', t, ', scenario ',SCEN, ' in simulation ', K, sep=''))
  # For loop over nsub
  for(s in 1:nsub){
    # Define two regions (which does nothing as there is no effect, )
    regions <- simprepSpatial(regions = 1, coord = TrueLocations, radius = list(TrueRadius), form ="cube", fading = 0)

    # Weighting structure.                                                                    MIND THE INDEX OF NOISE HERE!!!!!
    #   * Order = white, temporal, low-frequency, physyiological, task related and spatial.
    w <- Noise[[1]]

    # Base value
    base <- 5

    # Actual simulated data
    sim.data <- simVOLfmri(design=design.null, image=regions, base=base, dim=DIM, SNR=0.5,
                 type ="gaussian", noise= "mixture", spat="gaussRF", FWHM=2, weights=w, verbose = TRUE)
      # Transform it to correct dimension (Y = t x V)
      Y.data <- t(matrix(sim.data,ncol=nscan))

      rm(w, sim.data)
    ####************####
    #### ANALYZE DATA: 1e level
    ####************####

    # Fitting GLM model.
    model.lm <- lm(Y.data ~ x)
    b1 <- coef(model.lm)['x',]
    COPE[,s] <- b1
  }

  ####************####
  #### GROUP ANALYSIS: 2e level
  ####************####

  # Group COPE (average)
  GCOPE <- apply(COPE,1,mean,na.rm=TRUE)

  # Now we will do the OLS estimation of the variance
  GVARCOPE <- apply(COPE,1,var,na.rm=TRUE)

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
    SCOPE[,t] <- GCOPE
    SVARCOPE[,t] <- GVARCOPE
    STMAP[,t] <- GTMAP
    SHEDGE[,t] <- HedgeG
    SWEIGHTS[,t] <- weigFix

    rm(GCOPE,GVARCOPE,HedgeG,weigFix)
}

# Now calculate weighted average.
WeightedAvg <- (apply((SHEDGE*SWEIGHTS),1,sum))/(apply(SWEIGHTS,1,sum))

# Calculate variance of weighted average
varWeightAvg <- 1/apply(SWEIGHTS,1,sum)

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
CI.weightedVariance <- (apply((SWEIGHTS*(SHEDGE - WeightedAvg)^2),c(1),sum))/((nstud - 1) * apply(SWEIGHTS,1,sum))
CI.upper.weightVar <- matrix(WeightedAvg,ncol=1) + (qt(0.975,df=nstud-1) * sqrt(matrix(CI.weightedVariance,ncol=1)))
CI.lower.weightVar <- matrix(WeightedAvg,ncol=1) - (qt(0.975,df=nstud-1) * sqrt(matrix(CI.weightedVariance,ncol=1)))



##
###############
### Save objects: all in one list
###############
##

AllObjects <- list(
  'SCOPE' = SCOPE,
  'SVARCOPE' = SVARCOPE,
  'STMAP' = STMAP,
  'SHEDGE' = SHEDGE,
  'WeightedAvg' = WeightedAvg,
  'varWeightAvg' = varWeightAvg,
  'CI.upper.norm' = CI.upper.norm,
  'CI.lower.norm' = CI.lower.norm,
  'CI.upper.t' = CI.upper.t,
  'CI.lower.t' = CI.lower.t,
  'CI.upper.weightVar' = CI.upper.weightVar,
  'CI.lower.weightVar' = CI.lower.weightVar
  )

save(AllObjects, file=paste(wd,'/Results/',K,'/S',SCEN,'_SDG_AllObjects_K',K,sep=''))




















