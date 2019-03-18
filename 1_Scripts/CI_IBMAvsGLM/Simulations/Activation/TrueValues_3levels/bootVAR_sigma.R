####################
#### TITLE:     MA at third level: define true value through bootstrap
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
# Using bootstrap of the variance estimates at the second level.

# Could it be that ES at third level equals:
# var(GLM Y DATA AT THIRD LEVEL)/(var(sigmaHAT_scnd_level)) * h^2

# Bootstrap subjects into K studies?

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
  nBOOT <- 1000
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

# Data frame used to bootstrap
bootRES <- c()

# Data frame with results
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

# To know the true value of the ES, we will bootstrap at level 3 the squared
# estimates of the variances at the second level


##################
#### GENERATE DATA
##################

# Empty vector
bootEST <- c()

# For loop over the simulations
for(k in 1:K){
  print(k)
  # Empty vectors
  COPE <- VARCOPE <- array(NA, dim=nsub)
  STHEDGE <- STWEIGHTS <- array(NA, dim = nstud)
  VARHEDGE <- array(NA, dim = nstud)
  sigmaLVL2 <- c()
  LVL2sigmas <- c()

  # For loop over studies
  for(t in 1:nstud){
    # Generate the study-level (2e level) data using the model:
    # XG * Beta_G + E_G, where Beta_G comes from the meta-analysis level
    SLData <- XG * BOLDC + rnorm(n = nsub, mean = 0, sd = sqrt(sigma2B))
    
    # For loop over nsub
    for(s in 1:nsub){
      # signal
      Y.data <- (intcpt + X * SLData[s]) + rnorm(n = nscan, mean = 0, sd = sqrt(sigma2W))
      
      ####************####
      #### ANALYZE DATA: 1e level GLM
      ####************####
      # COPE (beta 1) --> fit GLM
      model.lm <- lm(Y.data ~ X)
      b1 <- coef(model.lm)['X']
      COPE[s] <- b1
      
      # VARCOPE --> estimate residual (we need to extend the design matrix with an intercept)
      xIN <- cbind(1,X)
      BETA <- coef(model.lm)
      res <- (t(Y.data - xIN %*% BETA) %*% (Y.data - xIN %*% BETA))/(nscan - 2)
      res <- diag(res)
      # Contrast: not interested in intercept
      CONTRAST <- matrix(c(0,1), nrow=1)
      # Calculate varcope
      VARCOPE[s] <- CONTRAST %*% (solve(t(xIN) %*% xIN )) %*% t(CONTRAST) %*% res
      
      # Clean objects
      rm(model.lm, b1, xIN, BETA, res, CONTRAST)
    }
    # Second level analysis: only save the estimate of the residual sigma
    LVL2sigmas <- c(LVL2sigmas, summary(lm(COPE ~ 1))$sigma)
  }
  # Now we run the bootstrap procedure where we generate new sets of studies
  # and then save the variance of the sigma estimates
    for(b in 1:nBOOT){
      # First sample subjects
      boot <- sample(x = COPE, size = nsub, replace = TRUE)
      # fit the OLS model and save residual sigma
      sigmaLVL2 <- c(sigmaLVL2, summary(lm(boot ~ 1))$sigma)
    }
    # Bootstrap estimate of variance of residual sigma
    # Save over all simulations
    bootEST <- c(bootEST, var(sigmaLVL2))
}

# Estimate of the variance of sigma at second level
varSIGMA <- mean(bootEST)

# So the true value of the ES at third level becomes
TrueES_lvl3 <- BOLDC/(sqrt(VARlvl3/varSIGMA) * NeuRRoStat::corrH(nsub, type = 'one.sample'))


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




















