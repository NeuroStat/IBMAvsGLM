####################
#### TITLE:   Vector simulation from time series to MA.
#### Contents:
####
#### Source Files:
#### First Modified: 23/01/2018
#### Notes:
#################

##
###############
### Notes
###############
##

# Simulate fMRI time series per subject, then combine using OLS
# Finally run MA. 
# Goal is to define true value per level and hence get the population effect true value.

# I will simulate on grid, but only save middle voxel (no smoothing)

##
###############
### Preparation
###############
##

# Take argument from master file
input <- commandArgs(TRUE)
# K'th simulation
hpcID <- try(as.numeric(as.character(input)[1]),silent=TRUE)
# Which machine
MACHINE <- try(as.character(input)[2],silent=TRUE)
# If no machine is specified, then it has to be this machine!
if(is.na(MACHINE)){
  MACHINE <- 'MAC'
  hpcID <- 1
}

# Implement for loop over r iterations here: hpcID goes from 1 to 100 in master file
rIter <- 10
startIndex <- try(1 + rIter * (hpcID - 1), silent = TRUE)
endIndex <- try(startIndex + (rIter - 1), silent = TRUE)


# Set WD: this is location where results are written
if(MACHINE=='HPC'){
  wd <- '/user/scratch/gent/gvo000/gvo00022/vsc40728/VectorMA'
}
if(MACHINE=='MAC'){
  wd <- '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/VectorMA'
}

# Load in libraries
library(AnalyzeFMRI)
library(lattice)
library(gridExtra)
library(oro.nifti)
library(ggplot2)
library(dplyr)
library(tibble)
library(reshape2)
library(RColorBrewer)
library(mvmeta)
library(metafor)
library(Hmisc)
library(devtools)
library(neuRosim)
library(NeuRRoStat)
library(fMRIGI)

# Data frame with results:
MAvec <- tibble(sim = integer(),
                      Wavg = numeric(),
                      sigma = numeric(),
                      nstud = numeric())
##
###############
### Functions
###############
##

##
###############
### Simulation parameters
###############
##


###################
#### Global variables
###################

# Number of subject: median sample size at 2018 = 28.5 (Poldrack et al., 2017)
nsub <- 29

# Number of studies
nstud <- 50

###################
#### Data characteristics
###################

# Signal characteristics
TR <- 2
nscan <- 200
total <- TR*nscan
on1 <- seq(1,total,40)
onsets <- list(on1)
duration <- list(20)

# Image characteristics
DIM <- c(9,9,9)
voxdim <- c(3.5, 3.5, 3.51) # Voxelsize
ext <- 1                    #  Extend
nregio <- 1

###################
#### Generate a design: GROUND TRUTH DESIGN
###################

# %BOLD change => fixed quantity
#   We will change the amount of noise to change effect size, Cohen's d
# See MAvsIBMA_Act_true_values.R on how we obtained values for Cohen's d
#   and the amount of noise within subjects to achieve these ES.
BOLDC <- 3

# Base of signal
base <- 100

# Spatial smoothing of signal
fwhm <- 8
sigma <- fwhm/sqrt(8*log(2))
width <- 5

# Subject parameters
TrueLocations <- c(5,5,5)
# In 1D-vector location
TrueLocVec_tmp <- array(NA, dim = DIM)
TrueLocVec_tmp[TrueLocations[1],
               TrueLocations[2],
               TrueLocations[3]] <- 1
TrueLocVec <- as.numeric(which(array(TrueLocVec_tmp, dim = prod(DIM)) == 1, 
                               arr.ind = TRUE))

# We generate a temporary design for getting a true signal
truthdesign <- neuRosim::simprepTemporal(1, 1, onsets = 1, effectsize = 1,
                               durations = 1, TR = 1, acc = 0.1)

# Now use this to get a sphere shaped area
area <- neuRosim::simprepSpatial(regions = 1, coord = list(TrueLocations),
                       radius = ext, form = "sphere", fading = 0)
truth <- neuRosim::simVOLfmri(design = truthdesign, image = area,
                    dim = DIM, SNR = 1, noise = "none")[,,,1]
GroundTruth <- ifelse(truth > 0, 1, 0)

#######################################
#### DESIGN AND SIGNAL TIME SERIES ####
#######################################

# Generating a design matrix
X <- neuRosim::simprepTemporal(total,1,onsets = onsets,
                     effectsize = 1, durations = duration,
                     TR = TR, acc = 0.1, hrf = "double-gamma")

# Generate time series for ONE active voxel: predicted signal, this is the design
pred <- neuRosim::simTSfmri(design=X, base=100, SNR=1, noise="none", verbose=FALSE)

# Now we create the BOLD signal by converting to % BOLD signal changes
# Need to be in appropriate scale
signal_BOLDC <- BOLDC * (pred-base) + base

# Now get the smoothed (raw) signal
Rawsignal <- GroundTruth %o% signal_BOLDC

## Design parameters
# Extend the design matrix with an intercept
xIN <- cbind(1,pred)

# Contrast: not interested in intercept
CONTRAST <- matrix(c(0,1),nrow=1)

# Calculate (X'X)^(-1) with contrast
design_factor <- CONTRAST %*% (solve(t(xIN) %*% xIN )) %*% t(CONTRAST)


##################
#### GENERATE DATA
##################

# Sigma of white noise
whiteSigma <- 100

# Between study variability
tau <- 0 

# Start 10 iterations (increases efficiency since iterations run very fast)
for(ID in startIndex:endIndex){
  
  # Set starting seed
  starting.seed <- pi*ID
  set.seed(starting.seed)

  # Empty vectors
  COPE <- VARCOPE <- array(NA,dim=c(prod(DIM),nsub))
  STHEDGE <- STWEIGHTS <- STCOPE <- STVARCOPE <- STVALUE <- array(NA,dim=c(prod(DIM),nstud))
    
  # For loop over studies
  for(t in 1:nstud){
    # Create the delta: subject specific true effect, using tau as between-study
    #   heterogeneity.
    # Distributed with mean true signal and variance tau
    StudData <- Rawsignal + array(
      array(rnorm(n = prod(DIM), mean = 0, sd = tau),
            dim = DIM), dim = c(DIM, nscan))
    # For loop over nsub
    for(s in 1:nsub){
      # Make white noise
      whiteNoise <- array(rnorm(n = (prod(DIM) * nscan), mean = 0,
                                sd = whiteSigma), dim = c(DIM, nscan))
  
      # Create image for this subject
      SubjData <- StudData + whiteNoise
      # plot(SubjData[5,5,5,], type = 'l')
      
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
      rm(model.lm, b1, xIN, BETA, res, CONTRAST)
    }
    
    ####************####
    #### GROUP ANALYSIS: 2e level using OLS
    ####************####
  
    # Group COPE
    STCOPE[,t] <- apply(COPE, 1, mean)
    
    # Group Y variable
    Gr.Y <- t(COPE)
    
    # Group VARCOPE:
    # First constant: (X'X)^-1
    # X is the design matrix, column of 1's
    GrX <- matrix(1, nrow = nsub)
    GrCt <- solve(t(GrX) %*% GrX)
    # Residuals
    GrRes <- (t(Gr.Y - GrX %*% matrix(STCOPE[,t], nrow = 1))) %*% 
              (Gr.Y - GrX %*% matrix(STCOPE[,t], nrow = 1))/
              (nsub - 1)
    GrRes <- diag(GrRes)
    # Denominator: checked using t.test
    STVARCOPE[,t] <- GrRes %*% GrCt
    
    # T value
    STVALUE[,t] <- STCOPE[,t] / sqrt(STVARCOPE[,t])
    
    # Clean objects
    rm(Gr.Y, GrRes)
  }
    
  # MA on voxel 365
  voxCOPE <- STCOPE[365,]
  voxVARCOPE <- STVARCOPE[365,]
  voxTval <- STVALUE[365,]
  
  # Hedges g
  voxHedgeG <- hedgeG(t = voxTval, N = nsub)
  
  # Variance of g
  voxVarG <- varHedge(voxHedgeG, N = nsub)
  
  # Weighted average
  WA <- as.numeric(rma(yi = voxHedgeG, vi = voxVarG, 
      method = 'DL')$beta)
  
  # Bind results in data frame
  MAvec <- data.frame(sim = ID,
             Wavg = WA,
             sigma = whiteSigma,
             nstud = nstud) %>%
    bind_rows(MAvec,. )

}

# Some extra objects to save
# Cope of one voxel of subjects from latest study in latest MA
subSimDat <- data.frame(Value = 
                c(COPE[TrueLocVec,],
                  VARCOPE[TrueLocVec, ]),
                param = rep(c('COPEsub', 'VARCOPEsub'), each = nsub),
                subID = rep(1:nsub, 2),
                voxID = TrueLocVec)

# Same for latest study in the latest MA
studSimDat <- data.frame(Value = 
                  c(STCOPE[TrueLocVec,],
                    STVARCOPE[TrueLocVec, ]),
                param = rep(c('COPEstud', 'VARCOPEstud'), each = nstud),
                studID = rep(1:nstud, 2),
                voxID = TrueLocVec)


##################
#### Save
##################

# Subject data
saveRDS(object = subSimDat, file = paste0(wd, '/Results/subSimDat_', hpcID, '.rda'))

# Study data
saveRDS(object = studSimDat, file = paste0(wd, '/Results/studSimDat_', hpcID, '.rda'))

# MA data
saveRDS(object = MAvec, file = paste0(wd, '/Results/MAvec_', hpcID, '.rda'))

