####################
#### TITLE:     GLM: third level variances
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



##
###############
### Preparation
###############
##


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
  nsim <- 5
  SCEN <- 1
}
# DataWrite directory: where all temp FSL files are written to
DataWrite <- try(as.character(input)[4],silent=TRUE)

# Set starting seed: it is the product of the amount of voxels,
# the number of studies and the number of subjects!
starting.seed <- 36865
set.seed(starting.seed)

# Set WD: this is location where results are written
if(MACHINE=='HPC'){
  wd <- '/user/scratch/gent/gvo000/gvo00022/vsc40728/IBMAvsMA'
}
if(MACHINE=='MAC'){
  wd <- '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/FLAME_check_3lvl'
}

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
library(oro.nifti)
library(dplyr)
library(tibble)
library(neuRosim)
library(NeuRRoStat)
library(fMRIGI)

# Data frame with results:
flameCheck <- tibble()

##
###############
### Functions
###############
##

# Function to generate a time series (for one voxel) 
# for each subject & study.
generateTimeSeries <- function(nscan, BETA, int, X, sigma2W){
  signal <- (int + X * BETA) + rnorm(n = nscan, mean = 0, sd = sqrt(sigma2W))
  return(signal)
}

# Function to gather results into tibbles
GetTibble <-function(data, nameParam, sim, DIM, BOLDC_p, sigmaW, sigmaM, nstud, tdof_t1){
  gather_data <- data.frame('sim' = as.integer(sim),
                            'voxel' = as.vector(1:prod(DIM)),
                            'value' = matrix(data, ncol = 1),
                            'parameter' = factor(nameParam,
                                                 levels = levels(saveParam)),
                            'BOLDC' = BOLDC_p,
                            'sigmaW' = sigmaW,
                            'sigmaM' = sigmaM,
                            'nstud' = nstud,
                            'FLAMEdf_3' = tdof_t1)
  return(as.tibble(gather_data))
}

##
###############
### Simulation parameters
###############
##

###################
#### Global options
###################

# Image characteristics
DIM <- trueMCvalues(ID = 'sim_act', keyword = 'DIM')

# Signal characteristics
TR <- trueMCvalues(ID = 'sim_act', keyword = 'TR')
nscan <- trueMCvalues('sim_act', 'nscan')
total <- trueMCvalues('sim_act', 'total')
on1 <- trueMCvalues('sim_act', 'on1')
onsets <- trueMCvalues('sim_act', 'onsets')
duration <- trueMCvalues('sim_act', 'duration')

# %BOLD change => fixed quantity
#   We will change the amount of noise to change effect size, Cohen's d
# See MAvsIBMA_Act_true_values.R on how we obtained values for Cohen's d
#   and the amount of noise within subjects to achieve these ES.
BOLDC <- trueMCvalues('sim_act', 'BOLDC')[2]

# Base of signal (i.e. intercept)
intcpt <- trueMCvalues('sim_act', 'base')

# Number of subject: median sample size at 2018 = 28.5 (Poldrack et al., 2017)
nsub <- trueMCvalues('sim_act', 'nsub')

# Number of studies
nstud <- 50

# Design factors: due to first level
design_lvl1 <- trueMCvalues('sim_act', 'design_factor')

# Second level design matrix and factor
XGmat <- matrix(1, nrow = nsub)
design_lvl2 <- solve(t(XGmat) %*% XGmat)


##############################
#### Simulation parameters
##############################

# Within-subject variance --> white noise: high, medium and low amount of noise (i.e. sigma_W)
sigma2W <- 500
  #trueMCvalues('sim_act', 'TrueSigma2W')[2]

# Between-subject variability, note: sigma^2_B/(sigma^2_W * design_factor) = 0.5
sigma2B <- 100
  # trueMCvalues('sim_act', 'TrueSigma2B')[2]

# Between-study variability: based on I^2 measure.
#I2 <- 62.61 / 100
# Work this out and we get:
#sigma2M <- (I2*sigma2W + I2*sigma2B) / (1 - I2)
sigma2M <- 60
TrueSigmaM <- sqrt(sigma2M)
#TrueSigmaM <- sqrt((I2*sigma2W*design_lvl1 + I2*sigma2B*design_lvl2) / (1 - I2))

###########################
###### GROUND TRUTH #######
###########################

# NOTE: partially same code as the one from the R package fMRIGI!

# We generate a temporary design for getting a true signal
truthdesign <- simprepTemporal(1, 1, onsets = 1, effectsize = 1,
                               durations = 1, TR = 1, acc = 0.1)

#######################################
#### DESIGN AND SIGNAL TIME SERIES ####
#######################################

### FIRST LEVEL ###
# Preparing the design matrix
X_prep <- simprepTemporal(total,1,onsets = onsets,
                          effectsize = 1, durations = duration,
                          TR = TR, acc = 0.1, hrf = "double-gamma")

# Generate the design matrix for each voxel
X <- simTSfmri(design=X_prep, base = 0, SNR = 1, noise = "none", verbose = FALSE)

### SECOND LEVEL ###
XG <- rep(1, nsub)

##################
#### GENERATE DATA
##################

# Progress
progR <- floor(seq(nsim/10, nsim, by = nsim/10))

# For loop over the K simulations
for(k in 1:nsim){
  if(k %in% progR) print(paste0('AT SIMULATION ', k))
  # Empty vectors
  COPE <- VARCOPE <- array(NA,dim=c(prod(DIM),nsub))
  STCOPE <- STVARCOPE <- array(NA,dim=c(prod(DIM), nstud))
  
  # We start by generating values for each study using the model: Y_M = X_M*Beta_M + E_M
  # We need a meta-analysis level (3e level) design matrix
  XM <- rep(1, nstud)
  # The effect at population level represents the %BOLD signal change
  Beta_M <- BOLDC
  MAData <- XM * Beta_M + rnorm(n = nstud, mean = 0, sd = sqrt(sigma2M))
  
  # For loop over studies
  for(t in 1:nstud){
    # Generate the study-level (2e level) data using the model:
    # XG * Beta_G + E_G, where Beta_G comes from the meta-analysis level
    SLData <- XG * MAData[t] + rnorm(n = nsub, mean = 0, sd = sqrt(sigma2B))
    
    # For loop over nsub
    for(s in 1:nsub){
      # Generate the signal in EACH voxel. 
      # First level model: X * Beta + E
      # No smoothing of noise as we are unable to calculate the true value of
      #   the effect size if we do so.
      # Run the generateTimeSeries function for each voxel.
      # Replicate is used to run the random number generator function several times.
      # Directly in correct dimension (Y = t x V).
      Y.data <- replicate(n = prod(DIM), generateTimeSeries(nscan = nscan,
                              BETA = SLData[s],
                              int = intcpt,
                              X = X,
                              sigma2W = sigma2W),
                          simplify = "array")
      
      ####************####
      #### ANALYZE DATA: 1e level GLM
      ####************####
      # COPE (beta 1) --> fit GLM
      model.lm <- lm(Y.data ~ X)
      b1 <- coef(model.lm)['X',]
      COPE[,s] <- b1
      
      # VARCOPE --> estimate residual (we need to extend the design matrix with an intercept)
      xIN <- cbind(1,X)
      BETA <- coef(model.lm)
      res <- (t(Y.data - xIN %*% BETA) %*% (Y.data - xIN %*% BETA))/(nscan - 2)
      res <- diag(res)
      # Contrast: not interested in intercept
      CONTRAST <- matrix(c(0,1), nrow=1)
      # Calculate varcope
      VARCOPE[,s] <- CONTRAST %*% (solve(t(xIN) %*% xIN )) %*% t(CONTRAST) %*% res
      
      # Clean objects
      rm(model.lm, b1, xIN, BETA, res, CONTRAST)
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

    # Clean up objects
    rm(GRCOPE4D,GRVARCOPE4D,command)
  }
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
  
  ### Now read in the SE
  #IBMA.SE <- sqrt(matrix(readNIfTI(paste(DataWrite,"/MA_stats/varcope1.nii",sep=""), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,],ncol=1))
  IBMA.SE <- sqrt(matrix(readNIfTI(paste(DataWrite,"/MA_stats/mean_random_effects_var1.nii",sep=""), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,],ncol=1))

  ########################################################################################################################################################################
  ########################################################################################################################################################################
  
  # Remove objects in DataWrite folder
  command <- paste0('rm -r ', DataWrite, '/*')
  system(command)
  
  ########################################################################################################################################################################
  ########################################################################################################################################################################
  
  # Save the SE in data frame
  flameCheck <- data.frame('sim' = k,
                           'voxel' = 1:prod(DIM),
                           'EstSE' = IBMA.SE,
                           'TrueSE' = TrueSigmaM) %>%
    as_tibble() %>%
    bind_rows(flameCheck, .)
}


##
###############
### Save object
###############
##

# saveRDS(flameCheck, 
#   file = '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/FLAME_check_3lvl/flameCheck.rda',
#         compress = TRUE)


##
###############
### Read object
###############
##

# Also possible to read back in
# flameCheck <- readRDS(file = '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/FLAME_check_3lvl/flameCheck.rda')

# Analysis
flameCheck %>%
  summarise(AvgEstSE = mean(EstSE),
            SDEstSE = sd(EstSE),
            TrueSE = mean(TrueSE)) 















