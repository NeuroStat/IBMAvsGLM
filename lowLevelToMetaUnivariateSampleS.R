####################
#### TITLE:     Simulate null data from low subject level to meta-analysis for one voxel only (UNIVARIATE APPROACH): increasing sample size and studies
#### Contents:
####
#### Source Files: HPC - Version
#### First Modified: 011/02/2016
#### Notes:
#################



##
###############
### Notes
###############
##

# Here, we check the univariate approach for increasing sample size and k studies.

K <- 333
SCEN <- 16


[6,]   22  422
[7,]   25  330
[8,]   25  333
[9,]   27  330
[10,]   27  422

# Reset working directory
rm(list=ls())
gc(verbose = FALSE)

# Take argument from master file
input <- commandArgs(TRUE)
  # K'th simulation
  K <- try(as.numeric(as.character(input)[1]), silent=TRUE)
  # Which scenario
  SCEN <- try(as.numeric(as.character(input)[2]), silent=TRUE)
  # Which machine
  MACHINE <- try(as.character(input)[3],silent=TRUE)
  if(!exists(MACHINE)){
    MACHINE <- 'MAC'
  }

# Set starting seed: it is the product of the amount of voxels, the number of studies and the number of subjects!
starting.seed <- 36864*K
set.seed(starting.seed)

# Set WD
if(MACHINE=='HPC'){
  wd <- '/user/scratch/gent/gvo000/gvo00022/vsc40728/Simulation'
}
if(MACHINE=='MAC'){
  wd <- '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/Take4'
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
on2 <- seq(20,total,40)
onsets <- list(on1,on2)
duration <- list(20,20)
effect.null <- list(0,0)                              ## No effect
effect <- list(1,1) 			                            ## Effect of 1 for designmatrix
DIM <- c(1,1,1)



# First make the data frame with the combinations of subjects and studies
is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
ss <- rep(seq(10,100,by=10),9)
k <- rep(seq(2,10),each=10)
ID <- is.wholenumber(ss/k)
OverView <- data.frame('Subjects' = ss, 'Studies' = k, 'Selection' = ID)
	OverView.Sel <- OverView[ID,-3]


# Now take the correct amount of subjects and studies (according to SCEN)
nsub <- OverView.Sel[SCEN,'Subjects']
nstud <- OverView.Sel[SCEN,'Studies']


####************####
#### Noise details
####************####
Noise <- list(
  'S1' = c(1,0,0,0,0,0)
  )



####************####
#### Subject/Study specific simulation details
####************####
# Subject parameters
TrueLocation <- c(1,1,1)
TrueWhiteNoise <- Noise[1]
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

save(SCOPE, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/UNI_SK_SCOPE_',K,sep=''))
save(SVARCOPE, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/UNI_SK_SVARCOPE_',K,sep=''))
save(STMAP, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/UNI_SK_STMAP_',K,sep=''))
save(SHEDGE, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/UNI_SK_SHEDGE_',K,sep=''))
save(WeightedAvg, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/UNI_SK_WeightedAvg_',K,sep=''))
save(varWeightAvg, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/UNI_SK_varWeightAvg_',K,sep=''))
save(CI.upper.norm, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/UNI_SK_CI.upper.norm_',K,sep=''))
save(CI.lower.norm, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/UNI_SK_CI.lower.norm_',K,sep=''))
save(CI.upper.t, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/UNI_SK_CI.upper.t_',K,sep=''))
save(CI.lower.t, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/UNI_SK_CI.lower.t_',K,sep=''))
save(CI.upper.weightVar, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/UNI_SK_CI.upper.weightVar_',K,sep=''))
save(CI.lower.weightVar, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/UNI_SK_CI.lower.weightVar_',K,sep=''))





