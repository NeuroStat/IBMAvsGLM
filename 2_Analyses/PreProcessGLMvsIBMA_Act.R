####################
#### TITLE:     Processing intermediate results on HPC: simulations of activation
#### Contents:
####
#### Source Files:
#### First Modified: 26/01/2018
#### Notes:
#################



##
###############
### Notes
###############
##

# This script loads in raw simulated data from the HPC, processes them into 
# summarized statistics and saves these objects. This is because objects are too large
# to transfer and process locally.

# Activation data!

# Blocked design for individual subjects.
# One condition.
# These N subjects are pooled using FLAME pooling
# The resulting images are converted to either Hedges' g and pooled using random effects meta-analysis.
# OR using 3e level GLM with again FLAME1.

# Data contains activation and between-study heterogeneity


##
###############
### Preparation
###############
##

# Take argument from master file
input <- commandArgs(TRUE)
# Which machine
MACHINE <- try(as.character(input)[1],silent=TRUE)
# If no machine is specified, then it has to be this machine!
if(is.na(MACHINE)){
  MACHINE <- 'MAC'
}

# Directories of the data for different simulations: here is raw data stored
if(MACHINE=='HPC'){
  DATAwd <- list(
    'Take[MAvsIBMA_Act]' = try(as.character(input)[2],silent=TRUE),
    'Take[MAvsIBMA_Act_RanInSl]' = try(as.character(input)[2],silent=TRUE))
}
if(MACHINE=='MAC'){
  DATAwd <- list(
    'Take[MAvsIBMA_Act]' = 
      "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/MAvsIBMA_act/Results_Parameters/Results",
    'Take[MAvsIBMA_Act_RanInSl]' = 
      "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/MAvsIBMA_act/Results_RanInSl/Results"
  )
}
NUMDATAwd <- length(DATAwd)
# Select your simulation setting here
currentWD <- 2

# Write Intermediate Results: this is location where summarized results are written
if(MACHINE=='HPC'){
  WIR <- '/user/scratch/gent/gvo000/gvo00022/vsc40728/IBMAvsMA/ProcessedResults'
}
if(MACHINE=='MAC'){
  WIR <- list(
    'Take[MAvsIBMA_Act]' = 
      '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/MAvsIBMA_act/Results_Parameters/ProcessedResults',
    'Take[MAvsIBMA_Act_RanInSl]' = 
      "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/MAvsIBMA_act/Results_RanInSl/ProcessedResults"
  )[[currentWD]]
}

# Date of today
date <- Sys.Date()

# Set starting seed
set.seed(1990)

# Print arguments
print(input)
print('-----')
print(MACHINE)

##
###############
### Functions
###############
##

# Libraries
library(ggplot2)
library(tibble)
library(tidyr)
library(dplyr)
library(AnalyzeFMRI)
library(lattice)
library(grid)
library(gridExtra)
library(oro.nifti)
library(reshape2)
library(RColorBrewer)
library(devtools)
library(neuRosim)
library(NeuRRoStat)
library(fMRIGI)


##
###############
### Simulation parameters
###############
##

# Number of conficence intervals
CIs <- c('MA-weightVar','GLM-t')
NumCI <- length(CIs)

# Data frame with number of simulations and subjects for current simulation
info <- data.frame('Sim' = c(1, 2),
                   'nsim' = c(1000, 500),
                   'nsub' = c(trueMCvalues('sim_act', 'nsub'),
                            trueMCvalues('sim_act', 'nsub')))
nsim <- info[currentWD,'nsim']
nsub <- info[currentWD,'nsub']

# However, if this machine, we select only 10 simulations (for testing code)
if(MACHINE == 'MAC') nsim <- 10

# Parameters that will be checked
saveParam <- factor(levels = c('CI.MA.upper.weightVar', 'CI.MA.lower.weightVar',
                               'MA.WeightedAvg',
                               'CI.IBMA.upper.t','CI.IBMA.lower.t', 'IBMA.COPE',
                               'CI.MA.weightedVariance', 'STHEDGE', 'ESTTAU',
                               'STWEIGHTS', 'STWEIGHTS_ran'))

# Data frame with simulation results:
MAvsIBMAres <- tibble(sim = integer(),
                      voxel = integer(),
                      value = numeric(),
                      parameter = saveParam,
                      sigma = numeric(),
                      eta = numeric(),
                      nstud = numeric(),
                      FLAMEdf_3 = numeric())

# Data frame with coverage results, averaged over all voxels:
coverage_Vox <- tibble(sim = integer(),
                      parameter = saveParam,
                      TrueD = numeric(),
                      eta = numeric(),
                      nstud = numeric(),
                      coverage = numeric(),
                      sdCov = numeric())

# Data frame with coverage results, NOT averaged over all voxels:
coverage_all <- tibble(sim = integer(),
                       voxel = integer(),
                       parameter = saveParam,
                       TrueD = numeric(),
                       eta = numeric(),
                       nstud = numeric(),
                       cov_IND = numeric())

# Data frame to calculate bias
bias_all <- tibble(sim = integer(),
                   voxel = integer(),
                    parameter = saveParam,
                    TrueD = numeric(),
                    eta = numeric(),
                    nstud = numeric(),
                    bias = numeric())

# Data frame to calculate CI length
CIlength_all <- tibble(sim = integer(),
                   voxel = integer(),
                   parameter = saveParam,
                   TrueD = numeric(),
                   eta = numeric(),
                   nstud = numeric(),
                   CIlength = numeric())

# Data frame to calculate CI length
EstVar_all <- tibble(sim = integer(),
                       voxel = integer(),
                       parameter = saveParam,
                       TrueD = numeric(),
                       eta = numeric(),
                       nstud = numeric(),
                       EstVar = numeric())

# Dimension of brain
DIM <- trueMCvalues('sim_act', 'DIM')

################
#### TRUE VALUES: load in R objects from the R package
################

# It is possible to use the estimated I^2 as an indicator for between-study
# heterogeneity. 
# --> if so: I2 = TRUE.
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
}

# Data frame with link between sigma and Cohen's d, will add this to the
# parameter combinations later on.
# NOTE: don't add eta here, as there is no link between d and eta!!
TrueParamDat <- data.frame(Nsub = trueMCvalues('sim_act', 'nsub'),
                           TrueD = trueMCvalues('sim_act', 'TrueD'),
                           TrueSigma = sqrt(trueMCvalues('sim_act', 'TrueSigma2W')),
                           TrueG = trueMCvalues('sim_act', 'TrueG'))

# Data frame with combinations of all simulations run
ParamComb <- expand.grid('TrueSigma' = sqrt(trueMCvalues('sim_act', 'TrueSigma2W')),
                         'eta' = eta_vec,
                         'nstud' = trueMCvalues('sim_act', 'nstud'))
NumPar <- dim(ParamComb)[1]

# Extend the true values with number of studies
TrueP_S <- TrueParamDat %>% 
  full_join(.,ParamComb, by = c('TrueSigma')) %>%
  # Add true COPE value
  mutate(TrueCOPE = trueMCvalues('sim_act', 'BOLDC'))

##
###############
### Data Wrangling
###############
##


# Load in all the data and summarise 
for(i in 1:nsim){
  print(paste0('In simulation ',i))
  # First load in the data of this simulation
  MAvsIBMAres <- readRDS(
    paste(DATAwd[[currentWD]],'/ActMAvsIBMA_',i, '.rda', sep='')) 
	
  #################################
  ###### PROCESSING COVERAGE ######
  #################################
  print('Calculating empirical coverage')
  
  # Data frame with estimate only (without the CI bounds)
  Estimate <- MAvsIBMAres %>%
    filter(parameter %in% c('MA.WeightedAvg', 'IBMA.COPE')) %>%
    left_join(., TrueP_S, by = c('sigma' = 'TrueSigma',
                                 'eta' = 'eta',
                                 'nstud' = 'nstud'))

  # Adding the CI bounds. Seperated per method.
  MAwide <- MAvsIBMAres %>%
    filter(parameter %in% c('CI.MA.upper.weightVar', 
                            'CI.MA.lower.weightVar')) %>%
    # Remove df column
    dplyr::select(-FLAMEdf_3) %>%
    # Need to spread the data frame with two extra columns
    tidyr::spread(., key = 'parameter', value = 'value') %>%
    # Rename to later bind in one data frame
    rename(., CI_upper = CI.MA.upper.weightVar,
           CI_lower = CI.MA.lower.weightVar) %>%
    # Add MA.WeightedAvg to dataframe
    full_join(
      filter(Estimate, parameter == 'MA.WeightedAvg') %>% 
        dplyr::select(-FLAMEdf_3),
      ., by = c("sim", "voxel", "sigma", "eta", "nstud"))
  
  # Repeat with GLM approach
  GLMwide <- MAvsIBMAres %>%
    filter(parameter %in% c('CI.IBMA.upper.t',
                            'CI.IBMA.lower.t')) %>%
    # Remove df column
    dplyr::select(-FLAMEdf_3) %>%
    tidyr::spread(., key = 'parameter', value = 'value') %>%
    # Rename
    rename(., CI_upper = CI.IBMA.upper.t,
           CI_lower = CI.IBMA.lower.t) %>%
    # Add MA.WeightedAvg to dataframe
    full_join(
      filter(Estimate, parameter == 'IBMA.COPE') %>% 
        dplyr::select(-FLAMEdf_3),
      ., by = c("sim", "voxel", "sigma", "eta", "nstud"))
  
  # Bind data frames
  ProcessedDat <- bind_rows(MAwide, GLMwide)
  
  # Add CI_coverage and summarise over all voxels
  EmpCov <- ProcessedDat %>% 
    # True value within CI limits?
    # NOTE: this is either hedges g (MA) or the true value for the BOLD signal (GLM)
    mutate(cov_IND = ifelse(parameter == "MA.WeightedAvg",
                            ifelse(TrueG >= CI_lower & TrueG <= CI_upper,1, 0),
                            ifelse(TrueCOPE >= CI_lower & TrueCOPE <= CI_upper,1, 0))) %>%
    # Drop variables that we do not need to calculate coverages
      # sigma dropped as info is in TrueD
    dplyr::select(voxel, parameter, TrueD, eta, nstud, cov_IND) %>%
    group_by(parameter, TrueD, eta, nstud) %>%
    summarise(coverage = mean(cov_IND),
              sdCov = sd(cov_IND)) %>%
    # Add ID of simulation
    mutate(sim = i)
  
  # Bind to data frame
  coverage_Vox <- bind_rows(coverage_Vox, EmpCov)
  
  # Also prepare object without summarizing over the voxels
  coverage_all <- ProcessedDat %>% 
    # True value within CI limits?
    # NOTE: this is either hedges g (MA) or the true value for the BOLD signal (GLM)
    mutate(cov_IND = ifelse(parameter == "MA.WeightedAvg",
                            ifelse(TrueG >= CI_lower & TrueG <= CI_upper,1, 0),
                            ifelse(TrueCOPE >= CI_lower & TrueCOPE <= CI_upper,1, 0))) %>%
    # Drop variables that we do not need to calculate coverages
    # sigma dropped as info is in TrueD
    dplyr::select(voxel, parameter, TrueD, eta, nstud, cov_IND) %>%
    # Add ID of simulation
    mutate(sim = i) %>%
    bind_rows(coverage_all,.)
  
  # Reset objects
  rm(MAwide, GLMwide, EmpCov)
  
  #################################
  ###### PROCESSING CI LENGTH #####
  #################################
  print('Calculating CI length')
  
  # Use ProcessedDat data frame
  CIlength_all <- ProcessedDat %>%
    mutate(CIlength = CI_upper - CI_lower) %>%
    # Select parameters
    dplyr::select(sim, voxel, parameter, TrueD, eta, nstud, CIlength) %>%
    bind_rows(CIlength_all,.)
  
  # Reset objects
  rm(ProcessedDat)
  
  #################################
  ###### PROCESSING STAN BIAS #####
  #################################
  print('Calculating standardized bias')
  
  # Use the 'Estimate' dataframe
  bias_all <- Estimate %>%
    # Mutate bias column, depending on MA or GLM
    mutate(bias = ifelse(parameter == 'MA.WeightedAvg',
                         (value - TrueG),
                         (value - TrueCOPE))) %>%
    # Select the parameters
    dplyr::select(sim, voxel, parameter, TrueD, eta, nstud, bias) %>%
    # Bind to data frame
    bind_rows(bias_all,.)
    
  # Reset objects
  rm(Estimate)
  
  ####################################
  ###### PROCESSING EST VARIANCE #####
  ####################################
  print('Calculating estimated variance')
  
  EstVar_all <- MAvsIBMAres %>%
    # Filter the IBMA approach: we need to re-calculate the variance!
    filter(parameter %in% c('IBMA.COPE', 'CI.IBMA.upper.t')) %>% 
    # Spread out CI upper bound and middle point
    spread(key = parameter, value = value) %>%
    # Now calculate the length of one side of the CI
    mutate(CIL = CI.IBMA.upper.t - IBMA.COPE) %>%
    # Use t-interval to find SE
    mutate(SE = CIL/(qt(0.975, df = FLAMEdf_3))) %>%
    # Square it
    mutate(EstVar = SE^2) %>%
    dplyr::select(sim, voxel, sigma, eta, nstud, IBMA.COPE, EstVar) %>%
    gather(key = 'parameter', value = 'value', IBMA.COPE) %>%
    dplyr::select(-value) %>%
    # Now add the data frame with MA and EstTau (estimated tau) to this
    bind_rows(.,
      MAvsIBMAres %>%
        filter(parameter %in% c('ESTTAU')) %>%
        # No need to square it as we saved squared tau already (see function 
        # ?NeuRRoStat::tau)!
        mutate(EstVar = value) %>%
        # Rename
        rename(param = parameter) %>%
        dplyr::select(-param, -value, -FLAMEdf_3) %>%
        mutate(parameter = 'MA.WeightedAvg')) %>%
    mutate(parameter = factor(parameter, levels = levels(saveParam))) %>%
    # Now have TrueD added instead of TrueSigma
    right_join(., TrueP_S, 
              by = c('sigma' = 'TrueSigma',
                    'eta' = 'eta',
                    'nstud' = 'nstud')) %>%
    dplyr::select(-sigma, -Nsub, -TrueG, -TrueCOPE) %>%
    # Bind to data frame
    bind_rows(EstVar_all, .)
}

# Save intermediate results
saveRDS(coverage_Vox, file = paste(WIR, '/coverage_Vox.rda', sep = ''))
saveRDS(coverage_all, file = paste(WIR, '/coverage_all.rda', sep = ''))
saveRDS(CIlength_all, file = paste(WIR, '/CIlength_all.rda', sep = ''))
saveRDS(bias_all, file = paste(WIR, '/bias_all.rda', sep = ''))
saveRDS(EstVar_all, file = paste(WIR, '/EstVar_all.rda', sep = ''))





