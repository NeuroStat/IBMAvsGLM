####################
#### TITLE:   Analysis of ==> Simulation: generate and fit general linear mixed model for fMRI
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

# Load in libraries
library(dplyr)
library(tibble)
library(lme4)
library(MASS)
library(neuRosim)

# Number of batches run
nbatch <- 100

# Location of data
locDat <- '/Users/hanbossier/Dropbox/PhD/PhDWork/Meta Analysis/R Code/Studie_Simulation/SimulationGit/1_Scripts/CI_IBMAvsGLM/Simulations/Activation/Variance_2levels/mixed_effects_glm_fMRI/Results'

##
###############
### Read in data
###############
##

# Empty data frame with simulation results
FitTotDat <- data.frame() %>% as_tibble()

# For loop over the batches (containing multiple simulation runs)
for(i in 1:nbatch){
  FitTotDat <-
    readRDS(file = paste(locDat, '/fMRMixEffglm_', i, '.rda', sep = '')) %>%
    bind_rows(FitTotDat, .)
}

##
###############
### Analysis
###############
##

# Checks
summary(FitTotDat$sim)
unique(FitTotDat$term)

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

























