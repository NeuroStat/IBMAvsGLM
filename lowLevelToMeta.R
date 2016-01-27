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
#   * 500 simulations
#   * Create 16x16x16 null-images for N subjects and K studies.
#   * Each image is created using the same design.
#   * No between study effects.
#   * Pool each study with ordinary OLS pooling method: T-maps.
#   * Transform each study to ES with formula used in FixRan study.
#   * Aggregate studies using fixed effects meta-analysis.
#   * Construct CI in each voxel.
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

# Set starting seed!!!!!!
set.seed(80*K)

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
  # Function to be fixed in neuRosim
  stimfunction<-function (totaltime, onsets, durations, accuracy)
{
  if (max(onsets) > totaltime) {
    stop("Mismatch between onsets and totaltime")
  }
  s <- rep(0, totaltime/accuracy)
  os <- onsets/accuracy
  dur <- durations/accuracy
  if (length(durations) == 1) {
    dur <- dur * rep(1, length(onsets))
  }
  else if (length(durations) != length(onsets)) {
    stop("Mismatch between number of onsets and number of durations.")
  }
  for (i in (1:length(onsets))) {
    if ((os[i] + dur[i]) <= totaltime/accuracy) {
      s[c(os[i]:(os[i] + dur[i]-1))] <- 1
    }
    else {
      s[c(os[i]:(totaltime/accuracy))] <- 1
    }
  }
  return(s)
}



# Load in functions from FixRan study
source_url('https://raw.githubusercontent.com/HBossier/FixRanStudyGit/master/Development/functions.R')


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
effect <- list(3,3) 			                            ## Effect of 1 for designmatrix
DIM <- c(16,16,16)


####************####
#### Scenario specific simulation details
####************####
whiteNoise <- c(1,1,rev(seq(0.1,0.9,by=0.1)))

whiteNoise <- list(
  'S1' = c(1,0,0,0,0,0),
  'S2' = c(1,0,0,0,0,0),
  'S3' = c(0.84,0.05,0.02,0.02,0.02,0.05),
  'S4' = c(0.64,0.15,0.02,0.02,0.02,0.15),
  'S5' = c(0.45,0.25,0.02,0.02,0.02,0.25),
  'S6' = c(0.24,0.35,0.02,0.02,0.02,0.35),
  'S7' = c(0.04,0.45,0.02,0.02,0.02,0.45)
  )



TrueWeights <- c(0.9, 0.05, 0, 0, 0, 0.05)

subjectfMRINoise <- function(TrueWeights,seed){
  # Set seed
  set.seed(seed)
  # Small amount of noise to low-frequency, physiological and task related
  LF <- phys <- task <- 0.02

  # Bounds
  UpperBound <- 0.94
  UpperWhite <- TrueWeights[1] + 0.2
  LowerBound <- 0
  LowerWhite <- TrueWeights[1] - 0.2
  # Start with white noise component
  white <- round(TrueWeights[1] + rnorm(1, mean = 0, sd = 0.2),2)
    # Now check if it is between the bounds
    while(white > UpperWhite || white > UpperBound || white < LowerWhite|| white < LowerBound){
      white <- round(TrueWeights[1] + rnorm(1, mean = 0, sd = 0.2),2)
    }
  # Now to temporal component: again create bounds
  UpperTemporal <- 0.94-white
  temporal <- round(TrueWeights[2] + rnorm(1, mean = 0, sd = 0.125),2)
    # Now check if it is between the bounds
    while(temporal > UpperTemporal || temporal > UpperBound || temporal < LowerBound){
      temporal <- round(TrueWeights[2] + rnorm(1, mean = 0, sd = 0.125),2)
    }
  # Spatial noise
  spatial <- 0.94 - white - temporal
  return(c(white,temporal,LF,phys,task,spatial))
}





####************####
#### Subject/Study specific simulation details
####************####
# Subject parameters
TrueLoc1 <- c(4,4,4)
TrueLoc2 <- c(10,10,10)
TrueWhiteNoise <- whiteNoise[SCEN]
TrueRadius <- 1
locations1 <- locations2 <- weights <- radius <- c()

COPE <- VARCOPE <- TMAP <- array(NA,dim=c(DIM,nsub))

# Loop over nsub to get the weights and the locations
for(s in 1:nsub){
  # Random locations
  loc.tmp1 <- TrueLoc1 + round(rnorm(3,0,2),0)
    while(any(loc.tmp1<1)){
      loc.tmp1 <- TrueLoc1 + round(rnorm(3,0,2),0)
    }
  loc.tmp2 <- TrueLoc2 + round(rnorm(3,0,2),0)
    while(any(loc.tmp2>16)){
      loc.tmp2 <- TrueLoc2 + round(rnorm(3,0,2),0)
    }
  locations1 <- rbind(locations1,loc.tmp1)
  locations2 <- rbind(locations2,loc.tmp2)

  # Radius
  radius.tmp <- TrueRadius + sample(c(0,1),size=2)
  radius <- rbind(radius,'rad' = radius.tmp)

  # Random noise components
  whiteNoise <- round(TrueWhiteNoise + rnorm(1,0,0.5),2)
  while(whiteNoise > 1 || whiteNoise < 0.5){
    whiteNoise <- round(TrueWhiteNoise + rnorm(1,0,0.5),2)
  }
  temporalNoise <- round((1-whiteNoise)/2 + rnorm(1,0,0.25),2)
  while(temporalNoise > c(1-whiteNoise) || temporalNoise < 0){
    temporalNoise <- round(TrueWhiteNoise + rnorm(1,0,0.5),2)
  }
  spatialNoise <- 1-whiteNoise-temporalNoise

  weights <- rbind(weights, c(whiteNoise, temporalNoise,0 ,0, 0, spatialNoise))
  rm(loc.tmp1,loc.tmp2,whiteNoise,temporalNoise,spatialNoise)
}

# Study parameters
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
#### GENERATE DATA: INCLUDE SUBJECTS --> STUDIES
####************####
# For loop over studies
for(t in 1:nstud){
  print(paste('------------------------- STUDY ', t,' -------------------------',sep=''))
  # For loop over nsub
  for(s in 1:nsub){
    print(paste('At subject, ', s, sep=''))
    coordinates <- list(c(locations1[s,]),c(locations2[s,]))
    # Define two regions (which does nothing as there is no effect, )
    regions <- simprepSpatial(regions = 2, coord = coordinates, radius = list(radius[s,1],radius[s,2]), form ="cube", fading = 0)
      rm(coordinates)

    # Weighting structure: white, temporal and spatial noise.
    #   * Order = white, temporal, low-frequency, physyiological, task related and spatial.
    #w <- weights[s,]
    w <- c(1,0,0,0,0,0)
    # Base value
    base <- 5

    # Actual simulated data
    sim.data <-  simVOLfmri(design=design.null, image=regions, base=base, dim=DIM, SNR=0.5,
                 type ="gaussian", noise= "mixture", spat="gaussRF", FWHM=2, weights=w, verbose = TRUE)
        rm(w)

    # 3D Gaussian Kernel over the 4D data
    # fwhm <- 2
    # sigma <- fwhm/(sqrt(8)*log(2))
    # smoothint <- -round(2*sigma):round(2*sigma)

    # for(i in 1:nscan){
    #  sim.data[,,,i] <- GaussSmoothArray(sim.data[,,,i], voxdim=c(1,1,1), ksize = length(smoothint), sigma = diag(sigma, 3))
    # }

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

      # Need a plot to check?
      PLOT <- FALSE
      if(isTRUE(PLOT)){
        levelplot(TMAP[,,,1])
        PVAL <- array(1-pt(TMAP,df=nscan-1),dim=DIM)
        levelplot(PVAL[,,,1])
        IDsign <- PVAL[,,,1] < 0.05
        PVAL[IDsign] <- 1
        PVAL[!IDsign] <- 0
        levelplot(PVAL)
      }
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
    gc(verbose = FALSE)
}


# Re-format arrays
SHEDGE.mat <- matrix(SHEDGE,ncol=nstud)
SWEIGHTS.mat <- matrix(SWEIGHTS,ncol=nstud)

# Now calculate weighted average.
WeightedAvg <- (apply((SHEDGE.mat*SWEIGHTS.mat),1,sum))/(apply(SWEIGHTS.mat,1,sum))

# Calculate variance of weighted average
varWeightAvg <- 1/apply(SWEIGHTS.mat,1,sum)

# Now calculate confidence intervals for weighted average based on assumption of normal distributed ES
CI.upper.norm <- matrix(WeightedAvg,ncol=1) + (1.96 * sqrt(matrix(varWeightAvg,ncol=1)))
CI.lower.norm <- matrix(WeightedAvg,ncol=1) - (1.96 * sqrt(matrix(varWeightAvg,ncol=1)))



##
###############
### Save objects
###############
##

save(SCOPE, file=paste(wd,'/Results/',K,'/WNSmSCOPE_',K,sep=''))
save(SVARCOPE, file=paste(wd,'/Results/',K,'/WNSmSVARCOPE_',K,sep=''))
save(STMAP, file=paste(wd,'/Results/',K,'/WNSmSTMAP_',K,sep=''))
save(SHEDGE, file=paste(wd,'/Results/',K,'/WNSmSHEDGE_',K,sep=''))
save(WeightedAvg, file=paste(wd,'/Results/',K,'/WNSmWeightedAvg_',K,sep=''))
save(varWeightAvg, file=paste(wd,'/Results/',K,'/WNSmvarWeightAvg_',K,sep=''))
save(CI.upper.norm, file=paste(wd,'/Results/',K,'/WNSmCI.upper.norm_',K,sep=''))
save(CI.lower.norm, file=paste(wd,'/Results/',K,'/WNSmCI.lower.norm_',K,sep=''))







