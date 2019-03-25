####################
#### TITLE:     MA vs GLM using FSL's FLAME on the COPE and VARCOPES: activation.
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

# We compare two techniques to aggregate fMRI data.
# The first is the outcome of transforming the second level
# GLM to an ES and calculate a weighted average using a random effects MA model.
# This is compared with a third level GLM, mixed effects model.

# Approach is to start from third level, generate data for second level, 
# then generate data at the first level. Then we model again to third level. 
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

# Record time
t1 <- Sys.time()

# Take argument from master file
input <- commandArgs(TRUE)
# K'th simulation
K <- try(as.numeric(as.character(input)[1]),silent=TRUE)
# Which scenario: GLM, MA using DL, MA using HE or MA using REML?
SCEN <- try(as.character(input)[2], silent=TRUE)
# Which machine
MACHINE <- try(as.character(input)[3],silent=TRUE)
# If no machine is specified, then it has to be this machine!
if(is.na(MACHINE)){
  MACHINE <- 'MAC'
  K <- 1
  SCEN <- 'REML'
}
# DataWrite directory: where all temp FSL files are written to
DataWrite <- try(as.character(input)[4],silent=TRUE)

# Check if SCEN contains one of the values that we need
if(!SCEN %in% c('GLM', 'DL', 'HE')){
  error('SCEN should be one of the following: GLM, DL or HE')
}

# Set starting seed: it is the product of the amount of voxels
# and the number of subjects!
starting.seed <- 21141*K
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
library(oro.nifti)
library(dplyr)
library(tibble)
library(neuRosim)
library(metafor)
library(NeuRRoStat)
library(fMRIGI)

# Parameters that could be saved if recorded
saveParam <- factor(levels = c('CI.MA.upper.weightVar', 'CI.MA.lower.weightVar',
                               'MA.WeightedAvg',
                               'CI.IBMA.upper.t','CI.IBMA.lower.t', 'IBMA.COPE',
                               'CI.MA.weightedVariance', 'STHEDGE', 'ESTTAU',
                               'STWEIGHTS', 'STWEIGHTS_ran','IBMA.SE'))

# Data frame with results:
MAvsIBMAres <- tibble(sim = integer(),
                      voxel = integer(),
                      value = numeric(),
                      parameter = saveParam,
                      SCEN = factor(levels = c('GLM', 'DL', 'HE', 'REML')),
                      BOLDC = numeric(),
                      sigmaW = numeric(),
                      sigmaM = numeric(),
                      nstud = numeric(),
                      FLAMEdf_3 = numeric())

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
GetTibble <-function(data, nameParam, SCEN, sim, DIM, BOLDC_p, sigmaW, sigmaM, nstud, tdof_t1){
  gather_data <- data.frame('sim' = as.integer(sim),
                            'voxel' = as.vector(1:prod(DIM)),
                            'value' = matrix(data, ncol = 1),
                            'parameter' = factor(nameParam,
                                          levels = levels(saveParam)),
                            'SCEN' = factor(SCEN,
                                          levels = c('GLM', 'DL', 'HE', 'REML')),
                            'BOLDC' = BOLDC_p,
                            'sigmaW' = sigmaW,
                            'sigmaM' = sigmaM,
                            'nstud' = nstud,
                            'FLAMEdf_3' = tdof_t1)
  return(as.tibble(gather_data))
}

# Function to get HE or REML estimate out of the metafor package
getHE_REML <- function(Y, W, method_tau = 'HE'){
  # Note that we will use mapply, thus the input for THIS function is a vector
  # that comes from a list. Therefore, no lists are used in this function.
          # Check if Y and W are lists
        #  if(class(Y) != "list") stop('Y has to be of class list!')
        #  if(class(W) != "list") stop('W has to be of class list!')
  
  # Now switch again from weights to within-study variance
        #lapply(W, function(x){1/x})
  wVar <- 1/W
  
  # Fit the linear model with HE or REML as estimator
  HE_REML_est <- metafor::rma(yi = Y, vi = wVar, method = method_tau)$tau2
  
  # Return the estimate
  return(HE_REML_est)
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
# We have a condition with no activation (i.e. 0) and one with.
BOLDC <- trueMCvalues('sim_act', 'BOLDC')

# Base of signal (i.e. intercept)
intcpt <- trueMCvalues('sim_act', 'base')

# Number of subject: median sample size at 2018 = 28.5 (Poldrack et al., 2017)
nsub <- trueMCvalues('sim_act', 'nsub')

##############################
#### Simulation parameters
##############################

# Within-subject variance --> white noise: high, medium and low amount of noise (i.e. sigma_W)
whiteSigma2W_vec <- trueMCvalues('sim_act', 'TrueSigma2W')

# Between-subject variability, note: sigma^2_B/(sigma^2_W * design_factor) = 0.5
TrueSigma2B_vec <- trueMCvalues('sim_act', 'TrueSigma2B')

# Between-study variability: based on I^2 measure.
TrueSigma2M_vec <- trueMCvalues('sim_act', 'TrueSigma2M')

# Number of studies in the MA.
nstud_vec <- trueMCvalues('sim_act', 'nstud')

# Data frame with combinations
ParamComb <- expand.grid('BOLDC' = BOLDC,
                'Sigma2W' = whiteSigma2W_vec,
                'Sigma2M' = TrueSigma2M_vec,
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

# Print scenario
print(paste('RUNNING SCENARIO: ', SCEN, sep = ''))

# For loop over the data generating parameters
for(p in 1:NumPar){
  print(paste('At parameter ', p, ' in simulation ', K, sep=''))
  
  # Select studies, amount of white noise, between-subject variability and between-study variability
  sigma2W <- ParamComb[p, 'Sigma2W']
  sigma2M <- ParamComb[p, 'Sigma2M']
  nstud <- ParamComb[p, 'nstud']
  # Between-subject variability moves with white noise
  sigma2B <- TrueSigma2B_vec[sigma2W == whiteSigma2W_vec]
  # Select the effect parameter
  BOLDC_p <- ParamComb[p, 'BOLDC']

  # Empty vectors
  COPE <- VARCOPE <- array(NA,dim=c(prod(DIM),nsub))
  STHEDGE <- STWEIGHTS <- STCOPE <- STVARCOPE <- array(NA,dim=c(prod(DIM), nstud))
  
  # We start by generating values for each study using the model: Y_M = X_M*Beta_M + E_M
  # We need a meta-analysis level (3e level) design matrix
  XM <- rep(1, nstud)
  # The effect at population level represents the %BOLD signal change
  Beta_M <- BOLDC_p
  MAData <- XM * Beta_M + rnorm(n = nstud, mean = 0, sd = sqrt(sigma2M))

  # For loop over studies
  for(t in 1:nstud){
    print(paste('At study ', t, ', parameter ', p, ' in simulation ', K, sep=''))

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
    
    ## WE WILL NEED TO HAVE THE ES WITH ITS VARIANCE FOR THE FIRST APPROACH:
    # Load in T-map
    STMAP <- readNIfTI(paste(DataWrite,"/study",t,"_stats/tstat1.nii",sep=""), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
    # Transform to an ES using hedgeG function, for each study
    HedgeG <- apply(matrix(STMAP,ncol=1), 1, FUN = NeuRRoStat::hedgeG, N = nsub, type = 'exact')
    # Calculate variance of ES
    VarianceHedgeG <- apply(matrix(HedgeG, ncol = 1), 1, FUN = NeuRRoStat::varHedgeT, N = nsub)
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
  
  # Now depending on the scenario, we have different parts of the script running
  if(SCEN %in% c('DL', 'HE', 'REML')){
  
    ####************####
    #### META-ANALYSIS: classical approach
    ####************####
    
    # First I need to make lists of Hedge g and its weights for using it in mapply
    STWEIGHTSL <- as.list(as.data.frame(t(STWEIGHTS)))
    STHEDGEL <- as.list(as.data.frame(t(STHEDGE)))
    
    # Estimate the between-study heterogeneity
    # First up: DL estimator
    if(SCEN == 'DL'){
      # Reason is that I combine matrices and per row I need Hedge g and weights
      # in a function.
      ESTTAU <- array(as.vector(mapply(NeuRRoStat::tau, Y = STHEDGEL,
                                       W = STWEIGHTSL, k = nstud)), dim = prod(DIM))
    }
    if(SCEN == 'HE'){
      # Estimate tau2 in each voxel, using the individual studies
      ESTTAU <- array(as.vector(mapply(getHE_REML, Y = STHEDGEL, W = STWEIGHTSL,
                                       method_tau = 'HE')), 
                      dim = prod(DIM))
        # If you want to check:    
        #  teY <- STHEDGEL[[728]]
        #  teW <- 1/STWEIGHTSL[[728]]
        #  metafor::rma(yi = teY, vi = teW, method = "HE")$tau2
    }
    if(SCEN == 'REML'){
      # Estimate tau2 in each voxel, using the individual studies
      ESTTAU <- array(as.vector(mapply(getHE_REML, Y = STHEDGEL, W = STWEIGHTSL,
                                       method_tau = 'REML')), 
                      dim = prod(DIM))
      # If you want to check:    
      #  teY <- STHEDGEL[[728]]
      #  teW <- 1/STWEIGHTSL[[728]]
      #  metafor::rma(yi = teY, vi = teW, method = "REML")$tau2
    }
    
    # Random effect weights: inverse of sum of within- and between-study variability
    STWEIGHTS_ran <- 1/((1/STWEIGHTS) + array(ESTTAU, dim = c(prod(DIM), nstud)))
    
    # Calculate weighted average.
    MA.WeightedAvg <- (apply((STHEDGE*STWEIGHTS_ran), 1, sum))/(apply(STWEIGHTS_ran, 1 ,sum))
    
    # CI for weighted average based on weighted variance CI
    CI.MA.weightedVariance <- (apply((STWEIGHTS_ran*(STHEDGE - MA.WeightedAvg)^2), c(1), sum))/((nstud - 1) * apply(STWEIGHTS_ran,1,sum))
    CI.MA.upper.weightVar <- matrix(MA.WeightedAvg,ncol=1) + (qt(0.975,df=nstud-1) * sqrt(matrix(CI.MA.weightedVariance,ncol=1)))
    CI.MA.lower.weightVar <- matrix(MA.WeightedAvg,ncol=1) - (qt(0.975,df=nstud-1) * sqrt(matrix(CI.MA.weightedVariance,ncol=1)))

    # Degrees of freedom, not really needed for this scenario, but I add it
    # as it will be in the data frame with the results.
    tdof_t1 <- nstud - 1
  }  
  ########################################################################################################################################################################
  ########################################################################################################################################################################

  # Scenario GLM
  if(SCEN == 'GLM'){
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
    
  }
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
    tmpObject <- try(get(levels(saveParam)[j]), silent = TRUE)
      if(class(tmpObject) =='try-error') next
    nameObject <- levels(saveParam)[j]
    MAvsIBMAres <- GetTibble(data = tmpObject, nameParam = nameObject, 
                             SCEN = SCEN, sim = K,
                      DIM = DIM, BOLDC_p = BOLDC_p,
                      sigmaW = sqrt(sigma2W), sigmaM = sqrt(sigma2M),
                      nstud = nstud, tdof_t1 =  tdof_t1) %>%
      bind_rows(MAvsIBMAres, .)
  }
}


##
###############
### Save object
###############
##
saveRDS(MAvsIBMAres, file = paste(wd,'/Results/',SCEN,'/ActMAvsIBMA_',K,'.rda', sep=''),
        compress = TRUE)

# Print time
print(Sys.time() - t1)








