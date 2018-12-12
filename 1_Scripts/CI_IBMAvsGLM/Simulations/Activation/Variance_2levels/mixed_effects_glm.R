####################
#### TITLE:   Simulation: generate and fit general linear mixed model
#### Contents:
####
#### Source Files:
#### First Modified: 13/07/2018
#### Notes:
#################

##
###############
### Notes
###############
##

# Generate data of N subjects each in M classes.
# Fit a linear mixed model and save the parameters.
# Are we able to extract the true values of the variance-covariance matrix
# of the fixed effects?

# The approach is to create the predictors for all subjects in all classes
# at the beginning. Then loop over all classes.
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
  wd <- '/user/scratch/gent/gvo000/gvo00022/vsc40728/MixedEffGLM/'
}
if(MACHINE=='MAC'){
  #wd <- '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/Variance_2lvl/'
  wd <- '/Users/hanbossier/Dropbox/PhD/PhDWork/Meta Analysis/R Code/Studie_Simulation/SimulationGit/1_Scripts/CI_IBMAvsGLM/Simulations/Activation/Variance_2levels/ResultsMixEffglm/'
}


# Load in libraries
library(AnalyzeFMRI)
library(lattice)
library(gridExtra)
library(oro.nifti)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(reshape2)
library(lme4)
library(MASS)
library(RColorBrewer)
library(mvmeta)
library(metafor)
library(devtools)
library(neuRosim)
library(NeuRRoStat)
library(fMRIGI)


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

# True values
## Class and sample size
numclass <- 40
NinClass <- 50

## Fixed parameters
beta0 <- 5
beta1 <- 2

## Random parameters
sigma_e <- 1
sigmab0 <- 1
sigmab1 <- 2

# Generate the variance covariance matrix of the random effects
var_cov_U <- rbind(c(sigmab0**2, 0), c(0, sigmab1**2))
# Generate the values for b0 and b1
B_matrix <- MASS::mvrnorm(numclass, mu = c(0,0), Sigma = var_cov_U)

# First we create empty vectors for X and V, as well as an empty matrix
# for the variance-covariance matrices of the fixed effects.
ComplX <- ComplV <- matrix(NA,nrow = 1, ncol = 1)
VarCovBeta_raw <- matrix(0, ncol = 2, nrow = 2)
Xlist <- Zlist <- list()

# Pre-define the true variance-covariance matrix for fixed effects parameters
for(i in 1:numclass){
  # Predictors for this class
  X <- cbind(1, round(runif(n = NinClass, min = 1, max = 20)))
  
  # Z-matrix for this class
  Z <- X
  
  # V-matrix
  V <- Z %*% var_cov_U %*% t(Z) +
    diag(sigma_e**2, NinClass)
  
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
  
  # Loop over the classes
  for(i in 1:numclass){
      # Loop over the subjects
      for(j in 1:NinClass){
      # Generate data using: X*beta + Z*u + e
        dat <- Xlist[[i]] %*% matrix(c(beta0, beta1), ncol = 1) + 
          Zlist[[i]] %*% matrix(c(B_matrix[i,1], B_matrix[i,2]), ncol = 1) +
          rnorm(n = NinClass, mean = 0, sd = sigma_e)
      
      # Add to data frame
      TotDat <- data.frame(Y = dat, X = Xlist[[i]][,2], class = i) %>% as_tibble() %>%
        bind_rows(TotDat,.)
    }
  }
  
  # Analysis
  fit <- lmer(Y ~ 1 + X + (1 + X|class), data = TotDat, REML = TRUE)
  FitTotDat <- broom::tidy(fit) %>% 
    # Add true SE
    left_join(.,SEBeta, by = 'term') %>%
    mutate(sim = r) %>%
    bind_rows(FitTotDat,.)
}

# Save R object
saveRDS(FitTotDat, file = paste0(wd, 'Results/MixEffglm_', hpcID, '.rda'))


##
###############
### Analysis
###############
##

# FitTotDat %>% 
#   mutate(filEst = ifelse(term %in% c('(Intercept)', 'X'),
#                 estimate, NA)) %>%
#   group_by(term) %>%
#   summarise(AvgEst = mean(estimate),
#             AvgSE = mean(std.error),
#             EmpSE = sd(filEst),
#             AvgTrueSE = mean(TrueSE))

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
                     'sd_(Intercept).class',
                     'sd_X.class')) %>%
  dplyr::select(term, estimate, sim) %>%
  left_join(.,data.frame('term' = c('sd_Observation.Residual',
                                    'sd_(Intercept).class',
                                    'sd_X.class'),
                         'TrueSD_ran' = c(sigma_e, sigmab0, sigmab1),
                         stringsAsFactors = FALSE),
            by = 'term') %>%
  group_by(term) %>%
  summarise(AvgEst = mean(estimate),
            AvgTrueSD = mean(TrueSD_ran))




