####################
#### TITLE:     Some checks for the GLM in simulation study.
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
library(MASS)



# Load in functions from FixRan study: THIS HAS TO COME AFTER ALL LIBRARIES ARE LOADED AS WE SOMETIMES FIX FUNCTIONS THAT ARE BUGGED IN THE PACKAGES
  source('~/Dropbox/PhD/PhDWork/Meta\ Analysis/R\ Code/Studie_FixRan/FixRanStudyGit.git/Development/functions.R')




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
DIM <- c(2,2,2)


####************####
#### Scenario specific simulation details
####************####


# Amount of subjects and studies
nsub <- 100
nstud <- 10

# Noise structure
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
  # Manual calculate COPE values:
    # Design matrix: extended with 1's for the intercept
    xIn <- cbind(1,x)
    # Contrast: 0,1 as we are not interested in b0
    CONTRAST <- c(0,1)
  manb1 <- CONTRAST %*% (solve(t(xIn) %*% xIn )) %*% t(xIn) %*% Y.data
  VisualCheckSub <- cbind(c(b1),c(manb1))


}

####************####
#### GROUP ANALYSIS: 2e level
####************####

# Group COPE (average)
GCOPE <- apply(COPE,1,mean,na.rm=TRUE)
  # Manual calculate group COPE values:
  Xg <- matrix(1,nrow=nsub)
    # We need the pseudo inverse of Xg, use function ginv from package MASS
      # pseudo inverse can also be calculated as:
      PseudInver <- solve(t(Xg) %*% Xg) %*% t(Xg)
  manGCOPE <- ginv(Xg) %*% t(COPE)
VisualCheckGCOPE <- cbind(c(GCOPE),c(manGCOPE))

# Now we will do the OLS estimation of the variance
GVARCOPE <- apply(COPE,1,var,na.rm=TRUE)

# TMAP
GTMAP <- GCOPE/sqrt(GVARCOPE/(nsub))

# Now check this T-map
manT.test <- apply(COPE,1,t.test)
  manT.val.tmp <- unlist(manT.test)
manT.val <- as.numeric(manT.val.tmp[names(manT.val.tmp)=='statistic.t'])
  VisualCheckTval <- data.frame(GTMAP,manT.val)








#### Some output from FSL's FLAMEO (FLAME1, mixed effects modeling at group stage)
  # We can calucalate CI around t-value by using the t-value and the VARCOPE.
  # However, we also get two files from FSL (zflame1uppertstat1 and zflame1lowertstat1). Not sure what these do.
    # CI around T-value
    TVAL <- matrix(readNIfTI(paste(DataWrite,"/MA_stats/tstat1.nii",sep=""), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,],ncol=1)
    SE <- sqrt(matrix(readNIfTI(paste(DataWrite,"/MA_stats/varcope1.nii",sep=""), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,],ncol=1))
      # Degrees of freedom:
      tdof_t1 <- readNIfTI(paste(DataWrite,"/MA_stats/tdof_t1.nii",sep=""), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[1,1,1]
    CI.upper.t <- TVAL +  (qt(0.975,df=tdof_t1) * SE)
    CI.lower.t <- TVAL -  (qt(0.975,df=tdof_t1) * SE)

    ZVAL <- matrix(readNIfTI(paste(DataWrite,"/MA_stats/zstat1.nii",sep=""), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,],ncol=1)
    upper <- ZVAL + (1.96 * SD)
    upper.zflame <- matrix(readNIfTI(paste(DataWrite,"/MA_stats/zflame1uppertstat1.nii",sep=""), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,],ncol=1)
      cbind(upper,upper.zflame)

    lower.zflame <- matrix(readNIfTI(paste(DataWrite,"/MA_stats/zflame1lowertstat1.nii",sep=""), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,],ncol=1)
    upper.zflame - ((upper.zflame-lower.zflame)/2)
    ZVAL
    matrix(TVAL,ncol=1)





