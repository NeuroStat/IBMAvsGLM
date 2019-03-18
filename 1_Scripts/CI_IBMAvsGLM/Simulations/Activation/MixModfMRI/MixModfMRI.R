####################
#### TITLE:   Simulation: 3-level general linear mixed model for fMRI
#### Contents:
####
#### Source Files:
#### First Modified: 12/12/2018
#### Notes:
#################

##
###############
### Notes
###############
##

# Generate time series of T time points in N subjects in K studies
# ---------------------------------------------------------------------------
# ----> change comments below!

# Fit a linear mixed model and save the parameters.
# Are we able to extract the true values of the variance-covariance matrix
# of the fixed and random effects?

# The approach is to create the design matrices for all subjects
# at the beginning. 
# The variance-covariance matrix of the fixed effects is calculated through:
# (sum over class of X'V^tX)^t
# ---------------------------------------------------------------------------


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
  wd <- '/user/scratch/gent/gvo000/gvo00022/vsc40728/MixedEff3LVLfMRI/'
}
if(MACHINE=='MAC'){
  wd <- '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/MixedEff3LVLfMRI/'
}

# Set starting seed
starting.seed <- 36*K
set.seed(starting.seed)

# Load in libraries
library(oro.nifti)
library(dplyr)
library(tibble)
library(lme4)
library(MASS)
library(neuRosim)
library(NeuRRoStat)
library(fMRIGI)
# library(AnalyzeFMRI)
# library(lattice)
# library(gridExtra)
# library(ggplot2)
# library(reshape2)
# library(Hmisc)
# library(devtools)

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

##
###############
### Simulation parameters
###############
##

###################
#### SIGNAL
###################

# Image characteristics
DIM <- c(1,1,1)

# Signal characteristics
TR <- trueMCvalues(ID = 'sim_act', keyword = 'TR')
nscan <- trueMCvalues('sim_act', 'nscan')
total <- trueMCvalues('sim_act', 'total')
on1 <- trueMCvalues('sim_act', 'on1')
onsets <- trueMCvalues('sim_act', 'onsets')
duration <- trueMCvalues('sim_act', 'duration')

# Number of subject: median sample size at 2018 = 28.5 (Poldrack et al., 2017)
nsub <- trueMCvalues('sim_act', 'nsub')

###################
#### GLM PARAMETER VALUES
###################

# ------------------ #
# Fixed effects parameters: 
# 1) Base of signal (i.e. intercept)
beta0 <- trueMCvalues('sim_act', 'base')
# 2) %BOLD change => fixed quantity
#   We will change the amount of noise to change effect size, Cohen's d
# See MAvsIBMA_Act_true_values.R on how we obtained values for Cohen's d
#   and the amount of noise within subjects to achieve these ES.
beta1 <- trueMCvalues('sim_act', 'BOLDC')

# ------------------ #
# Random effects parameters
#--#: Within subject parameter: white noise: high, medium and low amount 
sigma_e <- sqrt(trueMCvalues('sim_act', 'TrueSigma2W'))
#--#: Group level parameters
sigma_b0 <- 0
sigma_b1 <- 2
#--#: Study level parameter: sigma star of b0 and b1
sigmaS_b0 <- 0
# We estimated I^2: so from here we get the variance component between studies
I2_vec <- trueMCvalues('sim_act', 'I2')/100
# We put the other variance components in a data frame with all combinations.
# Then we calculate sigma star b1 using: sigS_b1 = (sigma_2^2*I2 + sigma_b1^2*I2)/(1 - I2)
# Note: no between subject variability on the intercepts (signal is normalized)
RandEffParam <- data.frame(
  expand.grid('sig_e' = sigma_e, 'I2' = I2_vec),
  'sig_b1' = sigma_b1) %>%
  mutate(sigS_b1 = sqrt((((sig_e^2 * I2) + (sig_b1^2 * I2)))/(1 - I2))) %>%
  mutate(sig_b0 = sigma_b0, 
         sigS_b0 = sigmaS_b0) %>%
  dplyr::select(sig_e, sig_b0, sig_b1, sigS_b0, sigS_b1, I2)

###################
#### TRUE SIGNAL
###################

# General list with the design matrix specification
Xgen <- simprepTemporal(totaltime = total,
                        regions = 1,
                        onsets = onsets, 
                        effectsize = beta1,
                        durations = duration,
                        TR = TR,
                        acc=0.1, hrf="double-gamma")
# Predicted signal
pred <- simTSfmri(design = Xgen, 
                  base = 0, 
                  SNR = 1,
                  noise = "none", verbose = FALSE) 

###################
#### TRUE STUDY VARIANCE-COVARIANCE MATRIX
###################

# Generate the variance covariance matrix of the random effects
var_cov_U <- rbind(c(sigma_b0**2, 0), c(0, sigma_b1**2))

# Pre-define the true variance-covariance matrix for fixed effects parameters:
# Do this for every value of sigma_2
# Empty data frame
SEBeta <- data.frame()
for(j in 1:length(sigma_e)){
  # First we create an empty matrix
  # for the variance-covariance matrices of the fixed effects.
  VarCovBeta_raw <- matrix(0, ncol = 2, nrow = 2)
  Xlist <- Zlist <- list()
  
  # Now loop over the subjects 
  for(i in 1:nsub){
    # Predictor for this subject
    X <- cbind(1, pred)
    
    # Z-matrix for this subject
    Z <- X
    
    # V-matrix
    V <- Z %*% var_cov_U %*% t(Z) +
      diag(sigma_e[j]**2, nscan)
    
    # Part of var-covar-beta matrix
    VarCovBeta_raw <- VarCovBeta_raw + t(X) %*% solve(V) %*% X
    
    # Save X and Z
    Xlist[[i]] <- X
    Zlist[[i]] <- Z
  }
  # Now calculate true variance-covariance matrix
  VarCovBeta <- solve(VarCovBeta_raw)
  
  # Gather in data frame
  SEBeta <- data.frame('term' = c('(Intercept)', 'X'),
                       'TrueSE' = sqrt(diag(VarCovBeta)),
                       'sig_e' = sigma_e[j], 
                       stringsAsFactors = FALSE) %>%
    bind_rows(SEBeta, .)
}

###################
#### STUDY LEVEL PARAMETERS
###################

# Change number of studies in the MA.
nstud_vec <- trueMCvalues('sim_act', 'nstud')

# Data frame with all combinations
ParamComb <- expand.grid('whiteSigma' = unique(RandEffParam$sig_e),
                         'I2' = unique(RandEffParam$I2),
                         'nstud' = nstud_vec)
NumPar <- dim(ParamComb)[1]


##
###############
### Simulation
###############
##


# For loop over the data generating parameters
for(p in 1:NumPar){
  print(paste('At parameter ', p, ' in simulation ', K, sep=''))
  
  # Select studies, amount of white noise, between-subject variability and between-study variability
  whiteSigma <- ParamComb[p, 'whiteSigma']
  I2_p <- ParamComb[p, 'I2']
  nstud <- ParamComb[p, 'nstud']
  sigmaS_b1 <- dplyr::filter(RandEffParam, sig_e == whiteSigma & I2 == I2_p) %>%
    dplyr::select(sigS_b1) %>% unlist() %>% as.numeric()

  # Empty vectors
  COPE <- VARCOPE <- array(NA,dim=c(prod(DIM),nsub))
  STHEDGE <- STWEIGHTS <- STCOPE <- STVARCOPE <- array(NA,dim=c(prod(DIM), nstud))
  
  # For loop over the studies
  for(s in 1:nstud){
    # First we generate the values for beta*_0 and beta*_1
    StudEff <- c(rnorm(n = 1, mean = 0, sd = sigmaS_b0),
                 rnorm(n = 1, mean = 0, sd = sigmaS_b1))
  
    # Empty data frame
    TotDat <- data.frame()
    
    # Generate the values for b0 and b1 in this study (with mean values for the study generated earlier)
    B_matrix <- MASS::mvrnorm(nsub, mu = StudEff, Sigma = var_cov_U)
    
    # Loop over the subjects
    for(i in 1:nsub){
      # Generate data using: X*beta + Z*u + e
      dat <- Xlist[[i]] %*% matrix(c(beta0, beta1), ncol = 1) + 
        Zlist[[i]] %*% matrix(c(B_matrix[i,1], B_matrix[i,2]), ncol = 1) +
        rnorm(n = nscan, mean = 0, sd = whiteSigma)
      
      # Add to data frame
      TotDat <- data.frame(Y = dat, X = Xlist[[i]][,2], subj = i) %>% as_tibble() %>%
        bind_rows(TotDat,.)
    }
    ####************####
    #### FIRST AND SECOND LEVEL ANALYSIS USING ONE STAGE MIXED EFFECTS
    ####************####
    fit <- lmer(Y ~ 1 + X + (1 + X|subj), data = TotDat, REML = TRUE)
    FitTotDat <- broom::tidy(fit) %>% 
      # Add true SE
      left_join(.,dplyr::filter(SEBeta, sig_e == whiteSigma), by = 'term') %>%
      mutate(sim = r) %>%
      bind_rows(FitTotDat,.)
    
  }
}



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


