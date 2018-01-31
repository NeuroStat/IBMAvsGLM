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
    'Take[MAvsIBMA_Act]' = try(as.character(input)[2],silent=TRUE))
}
if(MACHINE=='MAC'){
  DATAwd <- list(
    'Take[MAvsIBMA_Act]' = 
      "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/MAvsIBMA_act/Results_Parameters/Results")
}
NUMDATAwd <- length(DATAwd)
# Select your simulation setting here
currentWD <- 1 

# Write Intermediate Results: this is location where summarized results are written
if(MACHINE=='HPC'){
  WIR <- '/user/scratch/gent/gvo000/gvo00022/vsc40728/IBMAvsMA/ProcessedResults'
}
if(MACHINE=='MAC'){
  WIR <- '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/MAvsIBMA/ProcessedResults'
}

# Date of today
date <- Sys.Date()

# Set starting seed
set.seed(1990)

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
library(Hmisc)
library(devtools)
library(neuRosim)
library(scatterplot3d)
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
info <- data.frame('Sim' = c(1),
                   'nsim' = c(500),
                   'nsub' = trueMCvalues('sim_act', 'nsub'))
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
                      tau = numeric(),
                      nstud = numeric(),
                      FLAMEdf_3 = numeric())

# Data frame with coverage results, averaged over all voxels:
coverage_Vox <- tibble(sim = integer(),
                      parameter = saveParam,
                      TrueD = numeric(),
                      tau = numeric(),
                      nstud = numeric(),
                      coverage = numeric(),
                      sdCov = numeric())

# Data frame with coverage results, NOT averaged over all voxels:
coverage_all <- tibble(sim = integer(),
                       voxel = integer(),
                       parameter = saveParam,
                       TrueD = numeric(),
                       tau = numeric(),
                       nstud = numeric(),
                       cov_IND = numeric())

# Data frame to calculate bias
bias_all <- tibble(sim = integer(),
                   voxel = integer(),
                    parameter = saveParam,
                    TrueD = numeric(),
                    tau = numeric(),
                    nstud = numeric(),
                    bias = numeric())

# Data frame to calculate CI length
CIlength_all <- tibble(sim = integer(),
                   voxel = integer(),
                   parameter = saveParam,
                   TrueD = numeric(),
                   tau = numeric(),
                   nstud = numeric(),
                   CIlength = numeric())

# Data frame to calculate CI length
EstVar_all <- tibble(sim = integer(),
                       voxel = integer(),
                       parameter = saveParam,
                       TrueD = numeric(),
                       tau = numeric(),
                       nstud = numeric(),
                       EstVar = numeric())

# Dimension of brain
DIM <- trueMCvalues('sim_act', 'DIM')

################
#### TRUE VALUES: load in R objects from the R package
################

# Data frame of true parameter values
TrueParamDat <- data.frame(Nsub = trueMCvalues('sim_act', 'nsub'),
                           TrueD = trueMCvalues('sim_act', 'TrueD'),
                           TrueSigma = trueMCvalues('sim_act', 'TrueSigma'),
                           TrueG = trueMCvalues('sim_act', 'TrueG'),
                           Tau = trueMCvalues('sim_act', 'Tau'))

# Smoothed area
SmGT <- trueMCvalues('sim_act', 'SmGT')

# Switch to vector dimension
SmGT_v <- data.frame(voxID = 1:prod(DIM),
             Smooth = array(SmGT, dim = prod(DIM)))

# Masked GT area
MaskGT <- trueMCvalues('sim_act', 'MaskGT')

# Data frame with combinations of all simulations run
ParamComb <- expand.grid('TrueSigma' = trueMCvalues('sim_act', 'TrueSigma'),
                         'Tau' = trueMCvalues('sim_act', 'Tau'),
                         'nstud' = trueMCvalues('sim_act', 'nstud'))
NumPar <- dim(ParamComb)[1]

# Extend the true values with number of studies
TrueP_S <- TrueParamDat %>% 
  inner_join(.,ParamComb, by = c('TrueSigma', 'Tau'))

# We will need to extend the true parameter values with the smoothed value 
#     for each voxel. We have 3 true values (Hedges' g) per voxel.
SmParam <- data.frame(voxID = rep(1:prod(DIM), 
                        each = length(TrueParamDat$TrueG))) %>% 
  # Add the voxel ID to data frame
  bind_cols(., 
      # First replicate data frame for amount of voxels (use do.call rbind to do so)
      # Then bind columns
      do.call("rbind", replicate(prod(DIM), TrueParamDat, 
                                 simplify = FALSE))) %>% as.tibble(.) %>%
  # Not using number of subjects
  select(-Nsub) %>%
  # Join smoothed area vector to dataframe
  left_join(., SmGT_v, by = 'voxID') %>%
  # Multiply TrueG with Smoothed value to obtain smoothed Hedges' g
  mutate(SmoothG = TrueG * Smooth) %>%
  # Add true COPE value as well as smoothed value
  mutate(TrueCOPE = trueMCvalues('sim_act', 'BOLDC'),
         SmoothCOPE = TrueCOPE * Smooth) %>%
  # Drop tau
  select(-Tau)


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
    left_join(., SmParam, by = c('voxel' = 'voxID',
                                 'sigma' = 'TrueSigma'))
  
  # Adding the CI bounds. Seperated per method.
  MAwide <- MAvsIBMAres %>%
    filter(parameter %in% c('CI.MA.upper.weightVar', 
                            'CI.MA.lower.weightVar')) %>%
    # Need to spread the data frame with two extra columns
    tidyr::spread(., key = 'parameter', value = 'value') %>%
    # Rename to later bind in one data frame
    rename(., CI_upper = CI.MA.upper.weightVar,
           CI_lower = CI.MA.lower.weightVar) %>%
    # Add MA.WeightedAvg to dataframe
    full_join(
      filter(Estimate, parameter == 'MA.WeightedAvg'),
      ., by = c("sim", "voxel", "sigma", "tau", "nstud"))
  
  # Repeat with GLM approach
  GLMwide <- MAvsIBMAres %>%
    filter(parameter %in% c('CI.IBMA.upper.t',
                            'CI.IBMA.lower.t')) %>%
    tidyr::spread(., key = 'parameter', value = 'value') %>%
    # Rename
    rename(., CI_upper = CI.IBMA.upper.t,
           CI_lower = CI.IBMA.lower.t) %>%
    # Add MA.WeightedAvg to dataframe
    full_join(
      filter(Estimate, parameter == 'IBMA.COPE'),
      ., by = c("sim", "voxel", "sigma", "tau", "nstud"))
  
  # Bind data frames
  ProcessedDat <- bind_rows(MAwide, GLMwide)
  
  # Add CI_coverage and summarise over all voxels
  EmpCov <- ProcessedDat %>% filter(SmoothG != 0) %>%
    # True smoothed value within CI limits?
    # NOTE: smoothed G or smoothed COPE
    mutate(cov_IND = ifelse(parameter == "MA.WeightedAvg",
                            ifelse(SmoothG >= CI_lower & SmoothG <= CI_upper,1, 0),
                            ifelse(SmoothCOPE >= CI_lower & SmoothG <= CI_upper,1, 0))) %>%
    # Drop variables that we do not need to calculate coverages
      # sigma dropped as info is in TrueD
    select(voxel, parameter, TrueD, tau, nstud, cov_IND) %>%
    group_by(parameter, TrueD, tau, nstud) %>%
    summarise(coverage = mean(cov_IND),
              sdCov = sd(cov_IND)) %>%
    # Add ID of simulation
    mutate(sim = i)
  
  # Bind to data frame
  coverage_Vox <- bind_rows(coverage_Vox, EmpCov)
  
  # Also prepare object without summarizing over the voxels
  coverage_all <- ProcessedDat %>% filter(SmoothG != 0) %>%
    # True smoothed value within CI limits?
    # NOTE: smoothed G or smoothed COPE
    mutate(cov_IND = ifelse(parameter == "MA.WeightedAvg",
                            ifelse(SmoothG >= CI_lower & SmoothG <= CI_upper,1, 0),
                            ifelse(SmoothCOPE >= CI_lower & SmoothG <= CI_upper,1, 0))) %>%
    # Drop variables that we do not need to calculate coverages
    # sigma dropped as info is in TrueD
    select(voxel, parameter, TrueD, tau, nstud, cov_IND) %>%
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
    filter(SmoothG != 0) %>%
    mutate(CIlength = CI_upper - CI_lower) %>%
    # Select parameters
    select(sim, voxel, parameter, TrueD, tau, nstud, CIlength) %>%
    bind_rows(CIlength_all,.)
  
  # Reset objects
  rm(ProcessedDat)
  
  #################################
  ###### PROCESSING STAN BIAS #####
  #################################
  print('Calculating standardized bias')
  # Use the 'Estimate' dataframe
  bias_all <- Estimate %>%
    # Only active voxels
    filter(SmoothG != 0) %>%
    # Mutate bias column, depending on MA or GLM
    mutate(bias = ifelse(parameter == 'MA.WeightedAvg',
                         (value - SmoothG),
                         (value - SmoothCOPE))) %>%
    # Select the parameters
    select(sim, voxel, parameter, TrueD, tau, nstud, bias) %>%
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
    mutate(SE = CIL/(qt(0.975,df = FLAMEdf_3))) %>%
    # Square it
    mutate(EstVar = SE^2) %>%
    select(sim, voxel, sigma, tau, nstud, IBMA.COPE, EstVar) %>%
    gather(key = 'parameter', value = 'value', IBMA.COPE) %>%
    select(-value) %>%
    # Now add the data frame with MA and EstTau (estimated tau) to this
    bind_rows(.,
      MAvsIBMAres %>%
        filter(parameter %in% c('ESTTAU')) %>%
        # No need to square it as we saved squared tau already (see function 
        # ?NeuRRoStat::tau)!
        mutate(EstVar = value) %>%
        # Rename
        rename(param = parameter) %>%
        select(-param, -value) %>%
        mutate(parameter = 'MA.WeightedAvg')) %>%
    mutate(parameter = factor(parameter)) %>%
    # Now have TrueD added instead of TrueSigma
    left_join(., select(SmParam, voxID, TrueD, TrueSigma), 
              by = c('voxel' = 'voxID',
                    'sigma' = 'TrueSigma')) %>%
    select(-sigma) %>%
    # Bind to data frame
    bind_rows(EstVar_all, .)
  
}

# Save intermediate results
saveRDS(coverage_Vox, file = paste(WIR, '/coverage_Vox.rda', sep = ''))
saveRDS(coverage_all, file = paste(WIR, '/coverage_all.rda', sep = ''))
saveRDS(CIlength_all, file = paste(WIR, '/CIlength_all.rda', sep = ''))
saveRDS(bias_all, file = paste(WIR, '/bias_all.rda', sep = ''))
saveRDS(EstVar_all, file = paste(WIR, '/EstVar_all.rda', sep = ''))





