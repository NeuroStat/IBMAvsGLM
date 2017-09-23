####################
#### TITLE:     MA vs IBMA using FSL's FLAME on the COPE and VARCOPES: activation.
#### Contents:
####
#### Source Files:
#### First Modified: 24/02/2016
#### Notes:
#################



##
###############
### Notes
###############
##

# In this scenario, we will compare the outcome of transforming the second level
# GLM to an ES and execute a meta-analsysis with the scenario in which you use a
# third level GLM.

# MEASURES:
#   * CI coverage
#   * Standardized bias
#   * Average CI length


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
      SCEN <- 1
    }
  # DataWrite directory: where all files are written to
  DataWrite <- try(as.character(input)[4],silent=TRUE)

# Set starting seed: it is the product of the amount of voxels, the number of studies and the number of subjects!
starting.seed <- 36865*K
set.seed(starting.seed)

# Set WD
if(MACHINE=='HPC'){
  wd <- '/user/scratch/gent/gvo000/gvo00022/vsc40728/IBMAvsMA'
}
if(MACHINE=='MAC'){
  wd <- '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/MAvsIBMA'
}
setwd(wd)

# Give path to FSL
if(MACHINE=='HPC'){
  fslpath <- ''
}
if(MACHINE=='MAC'){
  fslpath <- '/usr/local/fsl/bin/'
}

# DataWrite if machine = Mac
if(MACHINE=='MAC'){
  DataWrite <- '~/Desktop/IBMA2'
}


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

# Load in functions from FixRan study: THIS HAS TO COME AFTER ALL 
# LIBRARIES ARE LOADED AS WE SOMETIMES FIX FUNCTIONS THAT ARE BUGGED IN THE PACKAGES
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

###################
#### Global options
###################

# Image characteristics
DIM <- c(9,9,9)
voxdim <- c(3.5, 3.5, 3.51) # Voxelsize
ext <- 1                    #  Extend
nregio <- 1

# Signal characteristics
TR <- 2
nscan <- 200
total <- TR*nscan
on1 <- seq(1,total,40)
onsets <- list(on1)
duration <- list(20)

# Effect size
es1 <- 0.02
base <- 100
actes1 <- es1*base

# Spatial smoothing of signal
fwhm <- 8
sigma <- fwhm/sqrt(8*log(2))
width <- 5

# Sigma of white noise
whiteSigma <- 7

##############################
#### Scenario specific details
##############################

# At this moment: one amount of subjects and studies
nsub <- 100
nstud <- 5

# At the moment, only consider one noise structure:
Noise <- list(
  'S1' = c(1,0,0,0,0,0)
)

###################################
#### Subject/Study specific details
###################################
# Subject parameters
TrueLocations <- c(4,4,5)
TrueWhiteNoise <- Noise[1]						# MIND THE INDEX HERE!
TrueRadius <- 1
COPE <- VARCOPE <- array(NA,dim=c(prod(DIM),nsub))

#####################
#### Study parameters
#####################
STHEDGE <- STWEIGHTS <- STCOPE <- STVARCOPE <- array(NA,dim=c(prod(DIM),nstud))

###########################
###### GROUND TRUTH #######
###########################

# We generate a temporary design for getting a true signal
truthdesign <- simprepTemporal(1,1,onsets=1,effectsize = 1, durations=1, TR=1, acc=0.1)

# Now use this to get a sphere shaped area
area1 <- simprepSpatial(regions=1, coord=list(TrueLocations), radius=ext, form="sphere", fading=0)
truth1 <- simVOLfmri(design=truthdesign, image=area1, dim=DIM, SNR=1,noise="none")[,,,1]
GroundTruth <- ifelse(truth1 > 0, 1, 0)

#######################################
#### DESIGN AND SIGNAL TIME SERIES ####
#######################################

# Generating a design matrix
X <- simprepTemporal(total,1,onsets=onsets,effectsize = 100, durations=duration,
                     TR = TR, acc=0.1, hrf="double-gamma")

# Generate time series for ONE active voxel: predicted signal
pred <- simTSfmri(design=X, base=100, SNR=1, noise="none", verbose=FALSE)

# Now we can convert to % BOLD signal changes
design1 <- es1 * (pred-base) + base

# Smooth the GT and put it into the map
SmGT1 <- GaussSmoothArray(truth1, voxdim = voxdim ,ksize = width, sigma = diag(sigma,3))
SmoothGT <- SmGT1

# Now get the smoothed signal
Rawsignal1 <- SmGT1 %o% design1
Rawsignal <- Rawsignal1

##################
#### GENERATE DATA
##################

# For loop over studies
for(t in 1:nstud){
  print(paste('At study ', t, ', scenario ',SCEN, ' in simulation ', K, sep=''))
  # For loop over nsub
  for(s in 1:nsub){
    # Make white noise
    whiteNoise <- array(rnorm(n = (prod(DIM) * nscan), mean = 0,
                              sd = whiteSigma), dim = c(DIM, nscan))
    
    # And smooth
    smoothNoise <- GaussSmoothArray(whiteNoise, voxdim = voxdim,
                                    ksize = width, sigma = diag(sigma,3))
    
    # Create image for this subject
    SubjData <- Rawsignal + smoothNoise
    
    # Transform it to correct dimension (Y = t x V)
    Y.data <- t(matrix(SubjData,ncol=nscan))
    
    ####************####
    #### ANALYZE DATA: 1e level GLM
    ####************####
    # COPE (beta 1) --> fit GLM
    model.lm <- lm(Y.data ~ pred)
    b1 <- coef(model.lm)['pred',]
    COPE[,s] <- b1
    
    # VARCOPE --> estimate residual (we need to extend the design matrix with an intercept)
    xIN <- cbind(1,pred)
    BETA <- coef(model.lm)
    res <- (t(Y.data - xIN %*% BETA) %*% (Y.data - xIN %*% BETA))/(nscan - 2)
    res <- diag(res)
    # Contrast: not interested in intercept
    CONTRAST <- matrix(c(0,1),nrow=1)
    # Calculate varcope
    VARCOPE[,s] <- CONTRAST %*% (solve(t(xIN) %*% xIN )) %*% t(CONTRAST) %*% res
    
    # Clean objects
    rm(model.lm, b1, xIN, BETA,res,CONTRAST)
    
  }
  ####************####
  #### GROUP ANALYSIS: 2e level using FLAME
  ####************####
  
  # Write auxiliarly files to DataWrite. We need:
  # GRCOPE in nifti
  # GRVARCOPE in nifti
  # 4D mask
  # design.mat file
  # design.grp file
  # design.con file
  
  #----- 1 ----#
  ### Design.mat
  fileCon <- paste(DataWrite,"/design.mat",sep="")
  # Text to be written to the file
  cat('/NumWaves\t1
      /NumPoints\t',paste(nsub,sep=''),'
      /PPheights\t\t1.000000e+00
      
      /Matrix
      ',rep("1.000000e+00\n",nsub),file=fileCon)
  
  #----- 2 ----#
  ### Design.con
  fileCon <- file(paste(DataWrite,"/design.con", sep=""))
  writeLines('/ContrastName1	Group Average
             /NumWaves	1
             /NumContrasts	1
             /PPheights		1.000000e+00
             /RequiredEffect		5.034
             
             /Matrix
             1.000000e+00
             ',fileCon)
    close(fileCon)

      #----- 3 ----#
      ### Design.grp
    fileCon <- paste(DataWrite,"/design.grp",sep="")
    # Text to be written to the file
    cat('/NumWaves\t1
        /NumPoints\t',paste(nsub,sep=''),'
        
        /Matrix
        ',rep("1\n",nsub),file=fileCon)

  #----- 4 ----#
  ### COPE.nii
  GRCOPE4D <- nifti(img=array(COPE,dim=c(DIM,nsub)),dim=c(DIM,nsub),datatype = 16)
  writeNIfTI(GRCOPE4D, filename = paste(DataWrite,'/GRCOPE',sep=''),gzipped=FALSE)
  
  #----- 5 ----#
  ### VARCOPE.nii
  GRVARCOPE4D <- nifti(img=array(VARCOPE,dim=c(DIM,nsub)),dim=c(DIM,nsub),datatype = 16)
  writeNIfTI(GRVARCOPE4D, filename = paste(DataWrite,'/GRVARCOPE',sep=''),gzipped=FALSE)
  
  #----- 6 ----#
  ### mask.nii
  mask <- nifti(img=array(1, dim=c(DIM,nsub)), dim=c(DIM,nsub), datatype=2)
  writeNIfTI(mask, filename = paste(DataWrite,'/mask',sep=''),gzipped=FALSE)
  
  # FSL TIME!
  setwd(DataWrite)
  command <- paste(fslpath, 'flameo --cope=GRCOPE --vc=GRVARCOPE --mask=mask --ld=study',t,'_stats --dm=design.mat --cs=design.grp --tc=design.con --runmode=flame1', sep='')
  Sys.setenv(FSLOUTPUTTYPE="NIFTI")
  system(command)
  
  # Put the result of pooling subjects in a vector for the COPE and VARCOPE
  STCOPE[,t] <- readNIfTI(paste(DataWrite,"/study",t,"_stats/cope1.nii",sep=""), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
  STVARCOPE[,t] <- readNIfTI(paste(DataWrite,"/study",t,"_stats/varcope1.nii",sep=""), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
  
  ## WE WILL NEED TO HAVE THE ES WITH ITS VARIANCE FOR THE FIRST APPROACH:
  # Load in T-map
  STMAP <- readNIfTI(paste(DataWrite,"/study",t,"_stats/tstat1.nii",sep=""), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
  # Transform to an ES using hedgeG function, for each study
  HedgeG <- apply(matrix(STMAP,ncol=1),1,FUN=hedgeG,N=nsub)
  # Calculate variance of ES
  VarianceHedgeG <- apply(matrix(HedgeG,ncol=1),1,FUN=varHedge,N=nsub)
  # Weights of this study
  weigFix <- 1/VarianceHedgeG
  # Now put in a vector
  STHEDGE[,t] <- HedgeG
  STWEIGHTS[,t] <- weigFix
  
  # Clean up objects
  rm(GRCOPE4D,GRVARCOPE4D,command,weigFix,HedgeG,STMAP)
}
########################################################################################################################################################################
########################################################################################################################################################################

####************####
#### META-ANALYSIS: classical approach
####************####
# Calculate weighted average.
MA.WeightedAvg <- (apply((STHEDGE*STWEIGHTS),1,sum))/(apply(STWEIGHTS,1,sum))

# CI for weighted average based on weighted variance CI
CI.MA.weightedVariance <- (apply((STWEIGHTS*(STHEDGE - MA.WeightedAvg)^2),c(1),sum))/((nstud - 1) * apply(STWEIGHTS,1,sum))
CI.MA.upper.weightVar <- matrix(MA.WeightedAvg,ncol=1) + (qt(0.975,df=nstud-1) * sqrt(matrix(CI.MA.weightedVariance,ncol=1)))
CI.MA.lower.weightVar <- matrix(MA.WeightedAvg,ncol=1) - (qt(0.975,df=nstud-1) * sqrt(matrix(CI.MA.weightedVariance,ncol=1)))


########################################################################################################################################################################
########################################################################################################################################################################

####************####
#### IBMA: 3e level using FLAME
####************####

# Write auxiliarly files to DataWrite. We need:
# STCOPE in nifti
# STVARCOPE in nifti
# 4D mask
# design.mat file
# design.grp file
# design.con file

#----- 1 ----#
### Design.mat
fileCon <- paste(DataWrite,"/STdesign.mat",sep="")
# Text to be written to the file
cat('/NumWaves\t1
    /NumPoints\t',paste(nstud,sep=''),'
    /PPheights\t\t1.000000e+00
    
    /Matrix
    ',rep("1.000000e+00\n",nstud),file=fileCon)

#----- 2 ----#
### Design.con
fileCon <- file(paste(DataWrite,"/STdesign.con", sep=""))
writeLines('/ContrastName1	Group Average
           /NumWaves	1
           /NumContrasts	1
           /PPheights		1.000000e+00
           /RequiredEffect		5.034
           
           /Matrix
           1.000000e+00
           ',fileCon)
close(fileCon)

#----- 3 ----#
### Design.grp
fileCon <- paste(DataWrite,"/STdesign.grp",sep="")
# Text to be written to the file
cat('/NumWaves\t1
    /NumPoints\t',paste(nstud,sep=''),'
    
    /Matrix
    ',rep("1\n",nstud),file=fileCon)

#----- 4 ----#
### STCOPE.nii
STCOPE4D <- nifti(img=array(STCOPE,dim=c(DIM,nstud)),dim=c(DIM,nstud),datatype = 16)
writeNIfTI(STCOPE4D, filename = paste(DataWrite,'/STCOPE',sep=''),gzipped=FALSE)

#----- 5 ----#
### VARCOPE.nii
STVARCOPE4D <- nifti(img=array(STVARCOPE,dim=c(DIM,nstud)),dim=c(DIM,nstud),datatype = 16)
writeNIfTI(STVARCOPE4D, filename = paste(DataWrite,'/STVARCOPE',sep=''),gzipped=FALSE)

#----- 6 ----#
### mask.nii
mask <- nifti(img=array(1, dim=c(DIM,nstud)), dim=c(DIM,nstud), datatype=2)
writeNIfTI(mask, filename = paste(DataWrite,'/mask',sep=''),gzipped=FALSE)


# FSL TIME!
setwd(DataWrite)
command <- paste(fslpath, 'flameo --cope=STCOPE --vc=STVARCOPE --mask=mask --ld=MA_stats --dm=STdesign.mat --cs=STdesign.grp --tc=STdesign.con --runmode=flame1', sep='')
Sys.setenv(FSLOUTPUTTYPE="NIFTI")
system(command)


### Now CI around COPE
IBMA.COPE <- matrix(readNIfTI(paste(DataWrite,"/MA_stats/cope1.nii",sep=""), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,],ncol=1)
IBMA.SE <- sqrt(matrix(readNIfTI(paste(DataWrite,"/MA_stats/varcope1.nii",sep=""), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,],ncol=1))
# Degrees of freedom:
tdof_t1 <- readNIfTI(paste(DataWrite,"/MA_stats/tdof_t1.nii",sep=""), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[1,1,1]

CI.IBMA.upper.t <- IBMA.COPE +  (qt(0.975,df=tdof_t1) * IBMA.SE)
CI.IBMA.lower.t <- IBMA.COPE -  (qt(0.975,df=tdof_t1) * IBMA.SE)

########################################################################################################################################################################
########################################################################################################################################################################



##
###############
### Save objects
###############
##
ObjectsMAvsIBMA <- list(
  'CI.MA.upper.weightVar' = CI.MA.upper.weightVar,
  'CI.MA.lower.weightVar' = CI.MA.lower.weightVar,
  'MA.WeightedAvg' = MA.WeightedAvg,
  'CI.IBMA.upper.t' = CI.IBMA.upper.t,
  'CI.IBMA.lower.t' = CI.IBMA.lower.t,
  'IBMA.COPE' = IBMA.COPE,
  'CI.MA.weightedVariance' = CI.MA.weightedVariance,
  'STHEDGE' = STHEDGE,
  'STWEIGHTS' = STWEIGHTS
)
save(ObjectsMAvsIBMA, file = paste(wd,'/Results/',K,'/SCEN_',SCEN,'/ObjectsMAvsIBMA_',K,sep=''))



    
    
    
    





