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

# We compare the outcome of transforming the second level
# GLM to an ES and calculate a weighted average using a random effects MA model.
# This is compared with a third level GLM, mixed effects model.

# Activation and between study heterogeneity is added using random intercept + slopes.
# Signal in every voxel!


# MEASURES:
#   * CI coverage
#   * Standardized bias
#   * Average CI length


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
      SCEN <- 1
    }
  # DataWrite directory: where all temp FSL files are written to
  DataWrite <- try(as.character(input)[4],silent=TRUE)

# Set starting seed: it is the product of the amount of voxels,
  # the number of studies and the number of subjects!
starting.seed <- 36865*K
set.seed(starting.seed)

# Set WD: this is location where results are written
if(MACHINE=='HPC'){
  wd <- '/user/scratch/gent/gvo000/gvo00022/vsc40728/IBMAvsMA'
}
if(MACHINE=='MAC'){
  wd <- '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/MAvsIBMA'
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
library(AnalyzeFMRI)
library(lattice)
library(gridExtra)
library(oro.nifti)
library(ggplot2)
library(dplyr)
library(tibble)
library(reshape2)
library(RColorBrewer)
library(Hmisc)
library(MASS)
library(devtools)
library(neuRosim)
library(NeuRRoStat)
library(fMRIGI)

# Parameters that will be saved
saveParam <- factor(levels = c('CI.MA.upper.weightVar', 'CI.MA.lower.weightVar',
                 'MA.WeightedAvg',
                 'CI.IBMA.upper.t','CI.IBMA.lower.t', 'IBMA.COPE',
                 'CI.MA.weightedVariance', 'STHEDGE', 'ESTTAU',
                 'STWEIGHTS', 'STWEIGHTS_ran'))

# Data frame with results:
MAvsIBMAres <- tibble(sim = integer(),
                voxel = integer(),
                value = numeric(),
                parameter = saveParam,
                sigma = numeric(),
                eta = numeric(),
                nstud = numeric(),
                FLAMEdf_3 = numeric())

##
###############
### Functions
###############
##

# Function to generate a time series (for one voxel) with specified random intercepts and slopes
# for each subject & study.
generateTimeSeries <- function(nscan, BOLDC, int, X, 
                               d0_bsub, d1_bsub, g0_bstud, g1_bstud, sigma_wsub){
  signal_s <- (int + d0_bsub + g0_bstud) +
    ((BOLDC + d1_bsub + g1_bstud) * X) + 
    rnorm(n = nscan, mean = 0, sd = sigma_wsub)
  return(signal_s)
}

# Function to gather results into tibbles
GetTibble <-function(data, nameParam, sim, DIM, sigma, eta, nstud, tdof_t1){
  gather_data <- data.frame('sim' = as.integer(sim),
                            'voxel' = as.vector(1:prod(DIM)),
                            'value' = matrix(data, ncol = 1),
                            'parameter' = factor(nameParam,
                                                 levels = levels(saveParam)),
                            'sigma' = sigma,
                            'eta' = eta,
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
BOLDC <- trueMCvalues('sim_act', 'BOLDC')

# Base of signal (i.e. intercept)
intcpt <- trueMCvalues('sim_act', 'base')

# Number of subject: median sample size at 2018 = 28.5 (Poldrack et al., 2017)
nsub <- trueMCvalues('sim_act', 'nsub')

##############################
#### Simulation parameters
##############################

# Sigma of WITHIN-SUBJECT white noise: high, medium and low amount of noise (i.e. sigma_W)
whiteSigma_vec <- sqrt(trueMCvalues('sim_act', 'TrueSigma2W'))

# Sigma^2 of between subject noise, note: sigma^2_B/sigma^2_W = 0.5
#	--> corresponds to variance of random slopes b1_bsub
TrueSigma2B1_vec <- trueMCvalues('sim_act', 'TrueSigma2B1')

# Random intercepts between subjects
TrueSigma2B0 <- trueMCvalues('sim_act', 'TrueSigma2B0')

# It is possible to use the estimated I^2 as an indicator for between-study
# heterogeneity. 
I2 <- FALSE
if(I2){
  # We estimated the amount of between-study heterogeneity using the I^2
  #   statistic. This is the amount of between-study hereogeneity in relation to 
  #   total variability. 
  I2_vec <- trueMCvalues('sim_act', 'I2')/100
  
  # Let us denote the amount of between-study variability in the GLM notation as eta^2
  #   Hence we assume: I^2 = (eta**2)/(eta**2 + whiteSigma**2)
  eta2_vec <- (I2_vec * (whiteSigma_vec**2 + TrueSigma2B1_vec))/(1-I2_vec)
  eta_vec <- sqrt(eta2_vec)
}

# Otherwise work with TrueSigma2starB1, which is an estimate of tau^2
if(!I2){
  # Have both slopes (we call it eta in this script)
  eta_vec <- sqrt(trueMCvalues('sim_act', 'TrueSigma2starB1'))
  # as intercepts
  TrueSigma2starB0 <- trueMCvalues('sim_act', 'TrueSigma2starB0')
}

# Change number of studies in the MA.
# However, need to find more sensible values.
nstud_vec <- trueMCvalues('sim_act', 'nstud')

# At the moment, only consider one noise structure:
    # White noise only
Noise <- list(
  'S1' = c(1,0,0,0,0,0)
)

# Data frame with combinations
ParamComb <- expand.grid('whiteSigma' = whiteSigma_vec,
                         'eta' = eta_vec,
                         'nstud' = nstud_vec)
NumPar <- dim(ParamComb)[1]

###################################
#### Subject/Study specific details
###################################
# Subject parameters

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

# Preparing the design matrix
X_prep <- simprepTemporal(total,1,onsets = onsets,
                     effectsize = 1, durations = duration,
                     TR = TR, acc = 0.1, hrf = "double-gamma")

# Generate the design matrix for each voxel
X <- simTSfmri(design=X_prep, base=0, SNR=1, noise="none", verbose=FALSE)


##################
#### GENERATE DATA
##################

# For loop over the data generating parameters
for(p in 1:NumPar){
  print(paste('At parameter ', p, ' in simulation ', K, sep=''))

  # Select studies, amount of white noise, between-subject variability and between-study variability
  whiteSigma <- ParamComb[p, 'whiteSigma']
  eta <- ParamComb[p, 'eta']
  nstud <- ParamComb[p, 'nstud']
  # Between-subject variability moves with white noise
  TrueSigma2B1 <- TrueSigma2B1_vec[whiteSigma == whiteSigma_vec]

  # Empty vectors
  COPE <- VARCOPE <- array(NA,dim=c(prod(DIM),nsub))
  STHEDGE <- STWEIGHTS <- STCOPE <- STVARCOPE <- array(NA,dim=c(prod(DIM), nstud))
  
  # First generate the variance-covariance matrix for all studies (g-matrix).
  # This one has a random intercept and slope.
  # Note that we have sigma^2 for random intercepts, but we use the square root for
  # between-study variability. Hence we need to square this one.
  ETA <- rbind(c(TrueSigma2starB0, 0), c(0, eta**2))
  
  # Now generate the study specific random intercept and slope values using this matrix
  Stud_matrix <- MASS::mvrnorm(nstud, mu=c(0,0), Sigma = ETA)

  # For loop over studies
  for(t in 1:nstud){
  print(paste('At study ', t, ', parameter ', p, ' in simulation ', K, sep=''))
    # Now create the variance-covariance matrix for each subject (d-matrix).
    SIGMA <- rbind(c(TrueSigma2B0, 0), c(0, TrueSigma2B1))
    
    # Generate the subject-specific values for intercept and slope using this d-matrix
    Sub_matrix <- MASS::mvrnorm(nsub, mu=c(0,0), Sigma = SIGMA)
    
    # For loop over nsub
    for(s in 1:nsub){
      # Generate the signal in EACH voxel. 
      # No smoothing of noise as we are unable to calculate the true value of
      #   the effect size if we do so.
      # Run the generateTimeSeries function for each voxel.
      # Replicate is used to run the random number generator function several times.
      # Directly in correct dimension (Y = t x V).
      Y.data <- replicate(n = prod(DIM), generateTimeSeries(nscan, BOLDC, intcpt, X, 
                   d0_bsub = Sub_matrix[s,1], d1_bsub = Sub_matrix[s,2],
                   g0_bstud = Stud_matrix[t, 1], g1_bstud = Stud_matrix[t,2],
                   sigma_wsub = whiteSigma),
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

    ## WE WILL NEED TO HAVE THE ES WITH ITS VARIANCE FOR THE FIRST APPROACH:
    # Load in T-map
    STMAP <- readNIfTI(paste(DataWrite,"/study",t,"_stats/tstat1.nii",sep=""), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
    # Transform to an ES using hedgeG function, for each study
    HedgeG <- apply(matrix(STMAP,ncol=1), 1, FUN = hedgeG, N = nsub)
    # Calculate variance of ES
    VarianceHedgeG <- apply(matrix(HedgeG,ncol=1), 1, FUN = NeuRRoStat::varHedgeT, N = nsub)
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

  # Estimate between-study heterogeneity: DL estimator
    # Need to make lists of Hedge g and the weights for using mapply
    # Reason is that I combine matrices and per row I need Hedge g and weights
    # in a function.
  STWEIGHTSL <- as.list(as.data.frame(t(STWEIGHTS)))
  STHEDGEL <- as.list(as.data.frame(t(STHEDGE)))
  ESTTAU <- array(as.vector(mapply(NeuRRoStat::tau, Y = STHEDGEL,
              W = STWEIGHTSL, k = nstud)), dim = prod(DIM))

  # Random effect weights
  STWEIGHTS_ran <- 1/(STWEIGHTS + array(ESTTAU, dim = c(prod(DIM), nstud)))

  # Calculate weighted average.
  MA.WeightedAvg <- (apply((STHEDGE*STWEIGHTS_ran), 1, sum))/(apply(STWEIGHTS_ran, 1 ,sum))

  # CI for weighted average based on weighted variance CI
  CI.MA.weightedVariance <- (apply((STWEIGHTS_ran*(STHEDGE - MA.WeightedAvg)^2), c(1), sum))/((nstud - 1) * apply(STWEIGHTS_ran,1,sum))
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
  tdof_t1 <- readNIfTI(paste(DataWrite,"/MA_stats/tdof_t1.nii",sep=""), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[5,5,5]
  
  CI.IBMA.upper.t <- IBMA.COPE +  (qt(0.975,df=tdof_t1) * IBMA.SE)
  CI.IBMA.lower.t <- IBMA.COPE -  (qt(0.975,df=tdof_t1) * IBMA.SE)

  ########################################################################################################################################################################
  ########################################################################################################################################################################

  # Remove objects in DataWrite folder
  command <- paste0('rm -r ', DataWrite, '/*')
  system(command)

  ########################################################################################################################################################################
  ########################################################################################################################################################################

  # Create data frame with all info through looping over the factor with all
  #  parameters and bind to tibble.
  for(j in 1:length(levels(saveParam))){
    tmpObject <- get(levels(saveParam)[j])
    nameObject <- levels(saveParam)[j]
    MAvsIBMAres <- GetTibble(tmpObject, nameObject,
                             sim = K, DIM, whiteSigma, eta, nstud, tdof_t1) %>%
      bind_rows(MAvsIBMAres, .)
  }
}


##
###############
### Save object
###############
##
saveRDS(MAvsIBMAres, file = paste(wd,'/Results/ActMAvsIBMA_',K,'.rda', sep=''),
        compress = TRUE)










