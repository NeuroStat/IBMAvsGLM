####################
#### TITLE:   Simulation: generate and fit general linear mixed model for fMRI
#### Contents:
####
#### Source Files:
#### First Modified: 11/12/2018
#### Notes:
#################

##
###############
### Notes
###############
##

# Generate time series of T time points in N subjects (the classes)
# Fit a linear mixed model and save the parameters.
# Are we able to extract the true values of the variance-covariance matrix
# of the fixed and random effects?

# The approach is to create the design matrices for all subjects
# at the beginning. 
# The variance-covariance matrix of the fixed effects is calculated through:
# (sum over class of X'V^tX)^t


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
  wd <- '/user/scratch/gent/gvo000/gvo00022/vsc40728/MixedEffGLMfMRI/'
}
if(MACHINE=='MAC'){
  wd <- '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/Variance_2lvl_fMRI/'
}


# Load in libraries
library(dplyr)
library(tibble)
library(lme4)
library(MASS)
library(neuRosim)


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

# Number of subjects
numsubj <- 40

# Fixed parameters: 
    # beta0 = base
    # beta1 = 3% BOLD signal change (if base signal = 100)
beta0 <- 100
beta1 <- 3

# Random effects parameters: variance components
sigma_e <- 4
sigmab0 <- 0
sigmab1 <- 2

# fMRI paradigm: block design 10s ON/OFF + 100 scans
nscans <- 100
tr <- 2
total.time <- nscans * tr
dur <- 10
onsets <- seq(1, total.time, dur * 2)

# General X list with the design matrix specification
Xgen <- simprepTemporal(totaltime = total.time,
                     regions = 1,
                     onsets = onsets, 
                     effectsize = beta1,
                     durations=dur,
                     TR = tr,
                     acc=0.1, hrf="double-gamma")
# Predicted signal
pred <- simTSfmri(design = Xgen, 
                  base = 0, 
                  SNR = 1,
                  noise = "none", verbose = FALSE) 

# Generate the variance covariance matrix of the random effects
var_cov_U <- rbind(c(sigmab0**2, 0), c(0, sigmab1**2))
# Generate the values for b0 and b1
B_matrix <- MASS::mvrnorm(numsubj, mu = c(0,0), Sigma = var_cov_U)

# First we create empty vectors for X and V, as well as an empty matrix
# for the variance-covariance matrices of the fixed effects.
ComplX <- ComplV <- matrix(NA,nrow = 1, ncol = 1)
VarCovBeta_raw <- matrix(0, ncol = 2, nrow = 2)
Xlist <- Zlist <- list()

# Pre-define the true variance-covariance matrix for fixed effects parameters
for(i in 1:numsubj){
  # Predictor for this subject
  X <- cbind(1, pred)

  # Z-matrix for this subject
  Z <- X
  
  # V-matrix
  V <- Z %*% var_cov_U %*% t(Z) +
    diag(sigma_e**2, nscans)
  
  # Part of var-covar-beta matrix
  VarCovBeta_raw <- VarCovBeta_raw + t(X) %*% solve(V) %*% X
  
  # Save X and Z
  Xlist[[i]] <- X
  Zlist[[i]] <- Z
}

# Now calculate true variance-covariance matrix
VarCovBeta <- solve(VarCovBeta_raw)

# Standard error of beta = sqrt(var(beta))
SEBeta <- data.frame('term' = c('(Intercept)', 'X'),
                     'TrueSE' = sqrt(diag(VarCovBeta)), 
                     stringsAsFactors = FALSE)


##
###############
### Simulation
###############
##

# Empty data frame with simulation results
FitTotDat <- data.frame() %>% as_tibble()

# For loop over the simulations
for(r in startIndex:endIndex){
  # Set starting seed
  starting.seed <- pi*r
  set.seed(starting.seed)
  
  # Empty data frame
  TotDat <- data.frame()
  
  # Loop over the subjects
  for(i in 1:numsubj){
    # Generate data using: X*beta + Z*u + e
    dat <- Xlist[[i]] %*% matrix(c(beta0, beta1), ncol = 1) + 
      Zlist[[i]] %*% matrix(c(B_matrix[i,1], B_matrix[i,2]), ncol = 1) +
      rnorm(n = nscans, mean = 0, sd = sigma_e)
    
    # Add to data frame
    TotDat <- data.frame(Y = dat, X = Xlist[[i]][,2], subj = i) %>% as_tibble() %>%
      bind_rows(TotDat,.)
  }
  
  # Analysis
  fit <- lmer(Y ~ 1 + X + (1 + X|subj), data = TotDat, REML = TRUE)
  FitTotDat <- broom::tidy(fit) %>% 
    # Add true SE
    left_join(.,SEBeta, by = 'term') %>%
    mutate(sim = r) %>%
    bind_rows(FitTotDat,.)
}

# Save R object
saveRDS(FitTotDat, file = paste0(wd, 'Results/fMRMixEffglm_', hpcID, '.rda'))


##
###############
### Analysis
###############
##

# Fixed part
FitTotDat %>% 
  filter(term %in% c('(Intercept)', 'X')) %>%
  # Add true estimate
  left_join(.,data.frame(term = c('(Intercept)', 'X'),
                         TrueEst = c(beta0, beta1), stringsAsFactors = FALSE),
            by = 'term') %>%
  group_by(term) %>%
  summarise(AvgEst = mean(estimate),
            TrueEst = mean(TrueEst),
            AvgSE = mean(std.error),
            TrueSE = mean(TrueSE))

# Random part
FitTotDat %>% 
  filter(term %in% c('sd_Observation.Residual',
                     'sd_(Intercept).subj',
                     'sd_X.subj')) %>%
  dplyr::select(term, estimate, sim) %>%
  left_join(.,data.frame('term' = c('sd_Observation.Residual',
                                    'sd_(Intercept).subj',
                                    'sd_X.subj'),
                         'TrueSD_ran' = c(sigma_e, sigmab0, sigmab1),
                         stringsAsFactors = FALSE),
            by = 'term') %>%
  group_by(term) %>%
  summarise(AvgEst = mean(estimate),
            AvgTrueSD = mean(TrueSD_ran))


