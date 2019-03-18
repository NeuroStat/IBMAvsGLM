####################
#### TITLE:     MA FSL's FLAME on the COPE and VARCOPES: check true value of ES
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

# Let us check the true value for the ES after three stages.

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
  K <- 250
  SCEN <- 1
}
# DataWrite directory: where all temp FSL files are written to
DataWrite <- try(as.character(input)[4],silent=TRUE)

# Set starting seed: it is the product of the amount of voxels,
starting.seed <- 36865
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
library(NeuRRoStat)
library(fMRIGI)

# Parameters that could be saved if recorded
saveParam <- factor(levels = c('CI.MA.upper.weightVar', 'CI.MA.lower.weightVar',
                               'MA.WeightedAvg',
                               'CI.IBMA.upper.t','CI.IBMA.lower.t', 'IBMA.COPE',
                               'CI.MA.weightedVariance', 'STHEDGE', 'ESTTAU',
                               'STWEIGHTS', 'STWEIGHTS_ran','IBMA.SE'))

# Data frame with results:
MACheckLVL2 <- MACheckLVL3 <- tibble()
checkLVL2 <- checkLVL3 <- c()

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

# Function to get t-value from t.test
getT <- function(x){
  tVal <- as.numeric(t.test(x)$statistic)
  return(tVal)
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
#DIM <- trueMCvalues(ID = 'sim_act', keyword = 'DIM')

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
BOLDC <- trueMCvalues('sim_act', 'BOLDC')[2]

# Base of signal (i.e. intercept)
intcpt <- trueMCvalues('sim_act', 'base')

# Number of subjects
nsub <- 100

##############################
#### Simulation parameters
##############################

# Within-subject variance
sigma2W <- 30

# Between-subject variability
sigma2B <- 10

# Between-study variability
sigma2M <- 20

# Number of studies in the MA
nstud <- 80


###################################
#### Subject/Study specific details
###################################
# Subject parameters

###########################
###### GROUND TRUTH #######
###########################


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

# solve(t(matrix(XG, ncol = 1)) %*% matrix(XG, ncol = 1))
# 1/sum(XG)
# 
# a <- sigma2B + (sigma2W * des_lvl1)
# b <- a/sum(XG)
# d <- b + sigma2M
# 3/sqrt(d)
# 
# 3/sqrt(((sigma2B + sigma2W*des_lvl1)*(1/sum(XG)) + sigma2M))
# 
# 
# CON <- matrix(c(0,1), nrow = 1)
# t(CON) %*% (t(X) %*% X)^(-1) %*% CON
# 
# des_lvl1 <- CON %*% (solve(t(Xex) %*% Xex )) %*% t(CON)

###########################
#### TRUE VALUE FOR ES ####
###########################

# First get the design factof of first level
# To do so, extend the design matrix with the intercept
Xex <- as.matrix(data.frame(intcpt,X), ncol = 2)
# Contrast: no interest in design matrix
CON <- matrix(c(0,1), nrow = 1)
des_lvl1 <- CON %*% (solve(t(Xex) %*% Xex )) %*% t(CON)

# Variance of data at second level
VARlvl2 <- sigma2B + sigma2W*des_lvl1

# Variance of data at third level
VARlvl3 <- ((sigma2B + sigma2W*des_lvl1)*(1/sum(XG)) + sigma2M)

# True value for ES at second level
TrueES_lvl2 <- BOLDC/sqrt(VARlvl2) * NeuRRoStat::corrH(Ne = nsub, type = "one.sample")

VARlvl3St <- (VARlvl3/VARlvl2) * (NeuRRoStat::corrH(nsub)**2)
BOLDC/sqrt(VARlvl3St)
BOLDC/sqrt(VARlvl3/VARlvl2)
3/(sqrt(VARlvl3)/sqrt(VARlvl2))
# True value for ES at third level
TrueES_lvl3 <- BOLDC/sqrt(VARlvl3)

sqrt(VARlvl3)/sqrt(sqrt((2*VARlvl2^2)/(nsub - 1)))

##################
#### GENERATE DATA
##################

# For loop over the simulations
for(k in 1:K){
  print(k)
  # Empty vectors
  #COPE <- VARCOPE <- array(NA,dim=c(prod(DIM),nsub))
  #STHEDGE_OLS <- STWEIGHTS_OLS <- STHEDGE <- STWEIGHTS <- array(NA,dim=c(prod(DIM), nstud))
  COPE <- VARCOPE <- array(NA, dim=nsub)
  STHEDGE <- STWEIGHTS <- array(NA, dim = nstud)
  VARHEDGE <- array(NA, dim = nstud)
  
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
      # Run the generateTimeSeries function for each voxel.
      # Replicate is used to run the random number generator function several times.
      # Directly in correct dimension (Y = t x V).
      # Y.data <- replicate(n = prod(DIM), generateTimeSeries(nscan = nscan,
      #                         BETA = SLData[s],
      #                         int = intcpt,
      #                         X = X,
      #                         sigma2W = sigma2W),
      #                     simplify = "array")
      
      Y.data <- (intcpt + X * SLData[s]) + rnorm(n = nscan, mean = 0, sd = sqrt(sigma2W))
      ####************####
      #### ANALYZE DATA: 1e level GLM
      ####************####
      # COPE (beta 1) --> fit GLM
      model.lm <- lm(Y.data ~ X)
      #b1 <- coef(model.lm)['X',]
      b1 <- coef(model.lm)['X']
      #COPE[,s] <- b1
      COPE[s] <- b1
      
      # VARCOPE --> estimate residual (we need to extend the design matrix with an intercept)
      xIN <- cbind(1,X)
      BETA <- coef(model.lm)
      res <- (t(Y.data - xIN %*% BETA) %*% (Y.data - xIN %*% BETA))/(nscan - 2)
      res <- diag(res)
      # Contrast: not interested in intercept
      CONTRAST <- matrix(c(0,1), nrow=1)
      # Calculate varcope
      #VARCOPE[,s] <- CONTRAST %*% (solve(t(xIN) %*% xIN )) %*% t(CONTRAST) %*% res
      VARCOPE[s] <- CONTRAST %*% (solve(t(xIN) %*% xIN )) %*% t(CONTRAST) %*% res
      
      # Clean objects
      rm(model.lm, b1, xIN, BETA, res, CONTRAST)
    }
    # 
    # ####************####
    # #### GROUP ANALYSIS: 2e level using FLAME
    # ####************####
    # 
    # # Write auxiliarly files to DataWrite. We need:
    # # GRCOPE in nifti
    # # GRVARCOPE in nifti
    # # 4D mask
    # # design.mat file
    # # design.grp file
    # # design.con file
    # 
    # #----- 1 ----#
    # ### Design.mat
    # fileCon <- paste(DataWrite,"/design.mat",sep="")
    # # Text to be written to the file
    # cat('/NumWaves\t1
    #     /NumPoints\t',paste(nsub,sep=''),'
    #     /PPheights\t\t1.000000e+00
    #     
    #     /Matrix
    #     ',rep("1.000000e+00\n",nsub),file=fileCon)
    # 
    # #----- 2 ----#
    # ### Design.con
    # fileCon <- file(paste(DataWrite,"/design.con", sep=""))
    # writeLines('/ContrastName1	Group Average
    #            /NumWaves	1
    #            /NumContrasts	1
    #            /PPheights		1.000000e+00
    #            /RequiredEffect		5.034
    #            
    #            /Matrix
    #            1.000000e+00
    #            ',fileCon)
    # close(fileCon)
    # 
    # #----- 3 ----#
    # ### Design.grp
    # fileCon <- paste(DataWrite,"/design.grp",sep="")
    # # Text to be written to the file
    # cat('/NumWaves\t1
    #     /NumPoints\t',paste(nsub,sep=''),'
    #     
    #     /Matrix
    #     ',rep("1\n",nsub),file=fileCon)
    # 
    # #----- 4 ----#
    # ### COPE.nii
    # GRCOPE4D <- nifti(img=array(COPE,dim=c(DIM,nsub)),dim=c(DIM,nsub),datatype = 16)
    # writeNIfTI(GRCOPE4D, filename = paste(DataWrite,'/GRCOPE',sep=''),gzipped=FALSE)
    # 
    # #----- 5 ----#
    # ### VARCOPE.nii
    # GRVARCOPE4D <- nifti(img=array(VARCOPE,dim=c(DIM,nsub)),dim=c(DIM,nsub),datatype = 16)
    # writeNIfTI(GRVARCOPE4D, filename = paste(DataWrite,'/GRVARCOPE',sep=''),gzipped=FALSE)
    # 
    # #----- 6 ----#
    # ### mask.nii
    # mask <- nifti(img=array(1, dim=c(DIM,nsub)), dim=c(DIM,nsub), datatype=2)
    # writeNIfTI(mask, filename = paste(DataWrite,'/mask',sep=''),gzipped=FALSE)
    # 
    # # FSL TIME!
    # setwd(DataWrite)
    # command <- paste(fslpath, 'flameo --cope=GRCOPE --vc=GRVARCOPE --mask=mask --ld=study',t,'_stats --dm=design.mat --cs=design.grp --tc=design.con --runmode=flame1', sep='')
    # Sys.setenv(FSLOUTPUTTYPE="NIFTI")
    # system(command)
    # 
    # ## WE WILL NEED TO HAVE THE ES WITH ITS VARIANCE FOR THE FIRST APPROACH:
    # # Load in T-map
    # STMAP <- readNIfTI(paste(DataWrite,"/study",t,"_stats/tstat1.nii",sep=""), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
    # # Transform to an ES using hedgeG function, for each study
    # HedgeG <- apply(matrix(STMAP,ncol=1), 1, FUN = NeuRRoStat::hedgeG, N = nsub, type = 'exact')
    # # Calculate variance of ES
    # VarianceHedgeG <- apply(matrix(HedgeG, ncol = 1), 1, FUN = NeuRRoStat::varHedgeT, N = nsub)
    # # Weights of this study
    # weigFix <- 1/VarianceHedgeG
    # # Now put in a vector
    # STHEDGE[,t] <- HedgeG
    # STWEIGHTS[,t] <- weigFix
    # 
    # # Clean up objects
    # rm(GRCOPE4D,GRVARCOPE4D,command,weigFix,HedgeG,STMAP)
    
    ####************####
    #### GROUP ANALYSIS: 2e level using OLS
    ####************####
    tval <- t.test(COPE)$statistic
    HedgeG <- as.numeric(NeuRRoStat::hedgeG(t = tval, N = nsub, type = 'exact'))
    VarianceHedgeG <- NeuRRoStat::varHedgeT(g = HedgeG, N = nsub)
    STHEDGE[t] <- HedgeG
    STWEIGHTS[t] <- 1/VarianceHedgeG
    VARHEDGE[t] <- VarianceHedgeG
    # Save level 2 ES
    checkLVL2 <- c(checkLVL2, HedgeG)

    # # First we get the t-value on the copes only (OLS approach)
    # TsOLS <- apply(COPE, 1, getT)
    # # Transform to an ES using hedgeG function, for each study
    # HedgeG_OLS <- apply(matrix(TsOLS,ncol=1), 1, FUN = NeuRRoStat::hedgeG, N = nsub, type = 'exact')
    # # Calculate variance of ES
    # VarianceHedgeG_OLS <- apply(matrix(HedgeG_OLS, ncol = 1), 1, FUN = NeuRRoStat::varHedgeT, N = nsub)
    # # Weights of this study
    # weigFix_OLS <- 1/VarianceHedgeG_OLS
    # # Now put in a vector
    # STHEDGE_OLS[,t] <- HedgeG_OLS
    # STWEIGHTS_OLS[,t] <- weigFix_OLS
  }
  
  # # Save the estimates at second level: FLAME
  # AvgLVL2 <- apply(STHEDGE, 2, mean)
  # # Save the estimates at second level: OLS
  # AvgLVL2_OLS <- apply(STHEDGE_OLS, 2, mean)
  
  # MACheckLVL2 <- data.frame('sim' = k,
  #     'study' = rep(1:nstud, 2),
  #     'LEVEL' = 2,
  #     'MODEL' = rep(c('FLAME', 'OLS'), each = nstud),
  #     'EstWAvg' = as.matrix(c(AvgLVL2, AvgLVL2_OLS)),
  #     'TrueES' = TrueES_lvl2) %>%
  #   as_tibble() %>%
  #   bind_rows(MACheckLVL2, .)
  
  ########################################################################################################################################################################
  ########################################################################################################################################################################
  
  ####************####
  #### META-ANALYSIS: classical approach using FLAMES model at second level
  ####************####
  
  # Estimate between-study heterogeneity: DL estimator
  # Need to make lists of Hedge g and the weights for using mapply
  # Reason is that I combine matrices and per row I need Hedge g and weights
  # in a function.
  # STWEIGHTSL <- as.list(as.data.frame(t(STWEIGHTS)))
  # STHEDGEL <- as.list(as.data.frame(t(STHEDGE)))
  # ESTTAU <- array(as.vector(mapply(NeuRRoStat::tau, Y = STHEDGEL,
  #                                  W = STWEIGHTSL, k = nstud)), dim = prod(DIM))
  # 
  # # Random effect weights: inverse of sum of within- and between-study variability
  # STWEIGHTS_ran <- 1/((1/STWEIGHTS) + array(ESTTAU, dim = c(prod(DIM), nstud)))
  # 
  # # Calculate weighted average.
  # MA.WeightedAvg <- (apply((STHEDGE*STWEIGHTS_ran), 1, sum))/(apply(STWEIGHTS_ran, 1 ,sum))
  # 
  # # Average over the voxels
  # MA.WeightedAvg_vox <- mean(MA.WeightedAvg)
  checkLVL3 <- c(checkLVL3, 
      metafor::rma(yi = c(STHEDGE), vi = c(VARHEDGE), method = "DL")$b[1])
  
  ########################################################################################################################################################################
  ########################################################################################################################################################################

  # # Remove objects in DataWrite folder
  # command <- paste0('rm -r ', DataWrite, '/*')
  # system(command)
  
  ########################################################################################################################################################################
  ########################################################################################################################################################################
  
  # # Save the weighted average in data frame
  # MACheckLVL3 <- data.frame('sim' = k,
  #                          'EstWAvg' = MA.WeightedAvg_vox,
  #                          'TrueES' = TrueES) %>%
  #   as_tibble() %>%
  #   bind_rows(MACheckLVL3, .)
}

mean(checkLVL2)
TrueES_lvl2
mean(checkLVL3)
TrueES_lvl3


##
###############
### Save object
###############
##
#saveRDS(MACheck, file = paste(wd,'/MACheck.rda', sep=''),
#        compress = TRUE)



MACheck %>%
  filter(sim >= 36) %>%
  summarise(AvgMA = mean(EstSE),
            TrueES = mean(TrueES))





NeuRRoStat::varHedgeT
NeuRRoStat::tau


gs <- rep(0.9, 30)
vargs<- NeuRRoStat::varHedgeT(g = gs, N = 20)
NeuRRoStat::tau(Y = gs, W = 1/vargs, k = 30)

test <- NeuRRoStat::varHedgeT(g = TrueES_lvl3, N = nstud)
TrueES_lvl3/test




















