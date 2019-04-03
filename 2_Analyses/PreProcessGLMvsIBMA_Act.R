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
# Two conditions: with 0% or 3% BOLD signal change.
# These N subjects are pooled using FLAME pooling
# The resulting images are converted to either Hedges' g and pooled using random effects meta-analysis.
# OR using 3e level GLM with again FLAME1.

# Data contains between-study heterogeneity.

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
  'Take[MAvsIBMA_Act_RanInSl]' = try(as.character(input)[2],silent=TRUE),
  'Take[GLMvsMA_wi_w_act]' = try(as.character(input)[2],silent=TRUE),
  'Take[HE]' = try(as.character(input)[2],silent=TRUE),
  'Take[Ratio]' = try(as.character(input)[2],silent=TRUE))
}
if(MACHINE=='MAC'){
  DATAwd <- list(
    'Take[MAvsIBMA_Act]' = 
      "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/MAvsIBMA_act/Results_Parameters/Results",
    'Take[MAvsIBMA_Act_RanInSl]' = 
      "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/MAvsIBMA_act/Results_RanInSl/Results",
    'Take[GLMvsMA_wi_w_act]' = "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/GLMvsMA_wi_w_act/Results",
    'Take[HE]' = "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/HE/Results",
    'Take[Ratio]' = "/Volumes/Elements/IBMAvsGLM/Ratios/"
  )
}
NUMDATAwd <- length(DATAwd)
# Select your simulation setting here
currentWD <- 5

# Which scenario: GLM, DL, HE or REML
SCEN <- try(as.character(input)[3], silent = TRUE)

# Select these values manually
if(MACHINE=='MAC'){
  SCEN <- 'HE'
}

# Write Intermediate Results: this is location where summarized results are written
if(MACHINE=='HPC'){
  WIR <- '/user/scratch/gent/gvo000/gvo00022/vsc40728/IBMAvsMA/ProcessedResults'
}
if(MACHINE=='MAC'){
  WIR <- list(
    'Take[MAvsIBMA_Act]' = 
      '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/MAvsIBMA_act/Results_Parameters/ProcessedResults',
    'Take[MAvsIBMA_Act_RanInSl]' = 
      "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/MAvsIBMA_act/Results_RanInSl/ProcessedResults",
    'Take[GLMvsMA_wi_w_act]' = 
      "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/GLMvsMA_wi_w_act/Results/ProcessedResults",
    'Take[HE]' = "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/HE/Results/ProcessedResults",
    'Take[Ratio]' = paste("/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/Ratio/", SCEN, "/Results/ProcessedResults", sep = "")
  )[[currentWD]]
}

# ---------------------------------------------------------------
## NEED TO UPDATE LOCATION OF RAW DATA AND PROCESSED RESULTS!!!
# ---------------------------------------------------------------

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
library(magrittr)
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
info <- data.frame('Sim' = c(1, 2, 3, 4, 5),
                   'nsim' = c(1000, 500, 1000, 1000, 1000),
                   'nsub' = c(trueMCvalues('sim_act', 'nsub'),
                            trueMCvalues('sim_act', 'nsub'),
                            trueMCvalues('sim_act', 'nsub'),
                            trueMCvalues('sim_act', 'nsub'),
                            trueMCvalues('sim_act', 'nsub')))
nsim <- info[currentWD,'nsim']
nsub <- info[currentWD,'nsub']

# However, if this machine, we select only 10 simulations (for testing code)
if(MACHINE == 'MAC') nsim <- 10

# Possible parameters over the various simulations that can be checked
ParamsSim <- c('CI.MA.upper.weightVar', 'CI.MA.lower.weightVar',
               'MA.WeightedAvg',
               'CI.IBMA.upper.t','CI.IBMA.lower.t', 'IBMA.COPE',
               'CI.MA.weightedVariance', 'STHEDGE', 'ESTTAU',
               'STWEIGHTS', 'STWEIGHTS_ran', 'IBMA.SE')

# Parameters that will be checked in each of the takes
if(currentWD %in% c(1,2)){
  saveParam <- factor(levels = ParamsSim[1:11])
}
if(currentWD %in% c(3,4,5)){
  saveParam <- factor(levels = ParamsSim[1:12])
}

# Dimension of brain
DIM <- trueMCvalues('sim_act', 'DIM')

# Ratio of between- over within subject variability
ratioBW_vec <- c(0.25, 0.5, 0.75)

# When saving values for each voxel (without averaging over the grid),
# we only use 4 values of the number of studies. Otherwise objects get too large
# and we cannot plot it either.
filtStud <- c(5,50)


##
###############
### Empty vectors
###############
##

# Data frame with simulation results:
MAvsIBMAres <- tibble(sim = integer(),
                      voxel = integer(),
                      value = numeric(),
                      parameter = saveParam,
                      SCEN = factor(levels = c('GLM','DL', 'HE','REML')),
                      BOLDC = numeric(),
                      sigmaW = numeric(),
                      sigmaM = numeric(),
                      nstud = numeric(),
                      FLAMEdf_3 = numeric())

# Info for each tibble
infoTibble <- tibble(sim = integer(),
                     parameter = saveParam,
                     SCEN = factor(levels = c('GLM','DL', 'HE','REML')),
                     TrueD = numeric(),
                     ratioBW = numeric(),
                     sigmaW = numeric(),
                     sigmaM = numeric(),
                     nstud = numeric())


# Data frame with coverage results, averaged over all voxels:
coverage_Vox <- infoTibble
coverage_Vox$coverage = numeric()
coverage_Vox$sdCov = numeric()

# Data frame with coverage results, NOT averaged over all voxels:
coverage_all <- infoTibble
coverage_all$voxel = integer()
coverage_all$cov_IND = numeric()

# Average CI length over all voxels
CIlength_Vox <- infoTibble
CIlength_Vox$AvCIlength <- numeric()

# Bias, averaged over all voxels
bias_Vox <- infoTibble
bias_Vox$Avbias <- numeric()

# Without averaging over voxels
bias_all <- infoTibble
bias_all$voxel <- integer()
bias_all$bias <- numeric()

# Data frame for average estimated square root of between-study variance over all voxels
EstVar_Vox <- infoTibble
EstVar_Vox$AvgEstSE <- numeric()

################
#### TRUE VALUES: load in R objects from the R package
################

# It is possible to use the estimated I^2 as an indicator for between-study
# heterogeneity. 
# --> if so: I2 = TRUE.
I2 <- TRUE
if(I2){
  # We estimated the amount of between-study heterogeneity using the I^2
  #   statistic. This is the amount of between-study hereogeneity in relation to 
  #   total variability. 
  I2_vec <- trueMCvalues('sim_act', 'I2')
}
# Otherwise work with TrueSigma2starB1, which is an estimate of tau^2
if(!I2){
  # Have both slopes (we call it eta in this script)
  eta_vec <- sqrt(trueMCvalues('sim_act', 'TrueSigma2starB1'))
}

# Data frame with link between sigma and Cohen's d, will add this to the
# parameter combinations later on.
# NOTE: don't add sigmaM here, as there is no link between d and sigmaM!!
TrueParamDat <- data.frame(Nsub = trueMCvalues('sim_act', 'nsub'),
                   TrueD = rep(trueMCvalues('sim_act', 'TrueD'), 3),
                   ratioBW = rep(c(0.25, 0.5, 0.75), each = 3),
                   TrueSigma = c(
                     sqrt(trueMCvalues('sim_act', 'TrueSigma2W', ratioBW = c(0.25))),
                     sqrt(trueMCvalues('sim_act', 'TrueSigma2W', ratioBW = c(0.50))),
                     sqrt(trueMCvalues('sim_act', 'TrueSigma2W', ratioBW = c(0.75)))),
                   TrueG = rep(trueMCvalues('sim_act', 'TrueG'), 3),
                   BOLDC = trueMCvalues('sim_act', 'BOLDC')[2])

# If preprocessing the take with null data, add extra part to TrueParamDat
if(currentWD %in% c(3,4,5)){
  tmp_TrueParamDat <- TrueParamDat
  tmp_TrueParamDat$BOLDC <- trueMCvalues('sim_act', 'BOLDC')[1]
  TrueParamDat %<>% bind_rows(
    tmp_TrueParamDat)
  # Update true values
  TrueParamDat %<>% mutate(ACT = ifelse(BOLDC > 0, 1, 0)) %>%
    mutate(TrueD = TrueD * ACT,
           TrueG = TrueG * ACT)
}

# Data frame with combinations of all simulations run: note that we expand the 
# grid with the act/no activation in the TrueParamDat data frame.
# Hence we do not need to add BOLDC here to the data frame!
ParamComb <- expand.grid('ratioBW' = ratioBW_vec,
            'sigmaW_level' = 1:length(trueMCvalues('sim_act', 'TrueSigma2W', ratioBW = 0.5)),
            'sigmaM' = sqrt(trueMCvalues('sim_act', 'TrueSigma2M')),
            'nstud' = trueMCvalues('sim_act', 'nstud')) %>%
    rowwise() %>%
    mutate(sigmaW = sqrt(trueMCvalues('sim_act', 'TrueSigma2W', 
                                      ratioBW = ratioBW)[sigmaW_level])) %>%
    ungroup() %>%
    dplyr::select(-sigmaW_level)
NumPar <- dim(ParamComb)[1]

# Extend the true values with number of studies
TrueP_S <- TrueParamDat %>% 
  full_join(.,ParamComb,
            by = c('TrueSigma' = 'sigmaW',
                   'ratioBW' = 'ratioBW')) %>% 
  # Add true COPE value (which is the same as BOLDC actually)
  mutate(TrueCOPE = ACT * trueMCvalues('sim_act', 'BOLDC')[2]) %>%
  as_tibble()


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
    paste(DATAwd[[currentWD]], '/', SCEN, '/ActMAvsIBMA_',i, '.rda', sep='')) 

  
  #################################
  ###### PROCESSING COVERAGE ######
  #################################
  print('Calculating empirical coverage')
  
  # Data frame with estimate only (without the CI bounds)
  Estimate <- MAvsIBMAres %>%
    filter(parameter %in% c('MA.WeightedAvg', 'IBMA.COPE')) %>%
    left_join(., TrueP_S, by = c('sigmaW' = 'TrueSigma',
                                 'ratioBW' = 'ratioBW',
                                 'sigmaM' = 'sigmaM',
                                 'nstud' = 'nstud',
                                 'BOLDC' = 'BOLDC'))

  # Adding the CI bounds. Seperated per method.
  if(SCEN == 'GLM'){
    MAwide <- as_tibble()
  }else{
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
        ., by = c("sim", "voxel", "SCEN", "ratioBW", "sigmaW", "sigmaB", "sigmaM", "nstud", "BOLDC"))
  }
  
  # Repeat with GLM approach
  if(SCEN != "GLM"){
    GLMwide <- as_tibble()
  }else{
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
        ., by = c("sim", "voxel", "SCEN", "ratioBW", "sigmaW", "sigmaB", "sigmaM", "nstud", "BOLDC"))
  }
  
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
      # sigmaW is not dropped anymore as we cross with 0 value TrueD
    dplyr::select(voxel, parameter, SCEN, TrueD, ratioBW, sigmaW, sigmaB, sigmaM, nstud, cov_IND) %>%
    group_by(parameter, SCEN, TrueD, ratioBW, sigmaW, sigmaB, sigmaM, nstud) %>%
    summarise(coverage = mean(cov_IND),
              sdCov = sd(cov_IND)) %>% 
    # Add ID of simulation
    mutate(sim = i)
  
  # Bind to data frame
  coverage_Vox <- bind_rows(coverage_Vox, EmpCov)
  
  # For each voxel, we reduce the amount of studies, otherwise
  # objects get too large.
  coverage_all <- ProcessedDat %>% 
    # True value within CI limits?
    # NOTE: this is either hedges g (MA) or the true value for the BOLD signal (GLM)
    mutate(cov_IND = ifelse(parameter == "MA.WeightedAvg",
              ifelse(TrueG >= CI_lower & TrueG <= CI_upper,1, 0),
              ifelse(TrueCOPE >= CI_lower & TrueCOPE <= CI_upper,1, 0))) %>%
    # Drop variables that we do not need to calculate coverages
    # sigmaW is not dropped anymore as we cross with 0 value TrueD
    dplyr::select(voxel, parameter, SCEN, TrueD, ratioBW, sigmaW, sigmaB, sigmaM, nstud, cov_IND) %>%
    # Filter out studies
    filter(nstud %in% filtStud) %>%
    # Maybe this won't work either...
    # Add ID of simulation
    mutate(sim = i) %>%
    bind_rows(coverage_all,.)

  # If we want for each voxel, we will need to run a separate script
  # Where we average over the simulations in batches.
  # Objects get too large
  # Also prepare object without summarizing over the voxels
  # ------------------------------------------------------------------------------
  # coverage_all <- ProcessedDat %>% 
  #   # True value within CI limits?
  #   # NOTE: this is either hedges g (MA) or the true value for the BOLD signal (GLM)
  #   mutate(cov_IND = ifelse(parameter == "MA.WeightedAvg",
  #                           ifelse(TrueG >= CI_lower & TrueG <= CI_upper,1, 0),
  #                           ifelse(TrueCOPE >= CI_lower & TrueCOPE <= CI_upper,1, 0))) %>%
  #   # Drop variables that we do not need to calculate coverages
  #   dplyr::select(voxel, parameter, TrueD, sigmaW, sigmaM, nstud, cov_IND) %>%
  #   # Add ID of simulation
  #   mutate(sim = i) %>%
  #   bind_rows(coverage_all,.)
  # ------------------------------------------------------------------------------
    
  # Reset objects
  rm(MAwide, GLMwide, EmpCov)
  
  #################################
  ###### PROCESSING CI LENGTH #####
  #################################
  print('Calculating CI length')
  
  # Use ProcessedDat data frame
  CIlength_Vox <- ProcessedDat %>%
    mutate(CIlength = CI_upper - CI_lower) %>%
    # Select parameters
    dplyr::select(sim, voxel, SCEN, parameter, TrueD, ratioBW, sigmaW, sigmaB, sigmaM, nstud, CIlength) %>%
    # Summarise over voxels
    group_by(sim, parameter, SCEN, TrueD, ratioBW, sigmaW, sigmaB, sigmaM, nstud) %>%
    summarise(AvCIlength = mean(CIlength)) %>%
    ungroup() %>%
    # Bind tot data frame
    bind_rows(CIlength_Vox,.)
  
  # Reset objects
  rm(ProcessedDat)
  
  #################################
  ###### PROCESSING STAN BIAS #####
  #################################
  print('Calculating standardized bias')
  
  # Use the 'Estimate' dataframe
  bias_Vox <- Estimate %>%
    # Mutate bias column, depending on MA or GLM
    mutate(bias = ifelse(parameter == 'MA.WeightedAvg',
                         (value - TrueG),
                         (value - TrueCOPE))) %>%
    # Select the parameters
    dplyr::select(sim, voxel, parameter, SCEN, TrueD, ratioBW, sigmaW, sigmaB, sigmaM, nstud, bias) %>%
    # Summarise over voxels
    group_by(sim, parameter, SCEN, TrueD, ratioBW, sigmaW, sigmaB, sigmaM, nstud) %>%
    summarise(Avbias = mean(bias)) %>%
    ungroup() %>%
    # Bind to data frame
    bind_rows(bias_Vox,.)
  
  # Without averaging over voxels
  bias_all <- Estimate %>%
    # Mutate bias column, depending on MA or GLM
    mutate(bias = ifelse(parameter == 'MA.WeightedAvg',
                         (value - TrueG),
                         (value - TrueCOPE))) %>%
    # Select the parameters
    dplyr::select(sim, voxel, parameter, SCEN, TrueD, ratioBW, sigmaW, sigmaB, sigmaM, nstud, bias) %>%
    # Filter studies
    filter(nstud %in% filtStud) %>%
    # Bind to data frame
    bind_rows(bias_all,.)
    
  # Reset objects
  rm(Estimate)
  
  ####################################
  ###### PROCESSING EST VARIANCE #####
  ####################################
  print('Calculating estimated variance')
  
  # First the average estimated between-study variance over all voxels
  EstVar_Vox <- MAvsIBMAres %>%
    # Filter estimated between-study variance components
    filter(parameter %in% c('IBMA.SE', 'ESTTAU')) %>%
    # Drop df
    dplyr::select(-FLAMEdf_3) %>%
    # Square root of tau as this is actually tau^2
    mutate(value = ifelse(parameter == 'ESTTAU', sqrt(value), value)) %>%
    # Round the sigma values
    mutate(sigmaW = round(sigmaW, 5),
           sigmaM = round(sigmaM, 5)) %>%
    # Summarise over voxels
    group_by(sim, parameter, SCEN, BOLDC, ratioBW, sigmaW, sigmaB, sigmaM, nstud) %>%
    summarise(AvgEstSE = mean(value)) %>% 
    # Ungroup
    ungroup() %>%
    #dplyr::select(-sigmaB, -SCEN) %>%
    # Now add the true value of Cohen's d
    left_join(., 
      dplyr::select(TrueP_S, -ACT, -TrueCOPE, -TrueG, -Nsub) %>%
        mutate(TrueSigma = round(TrueSigma, 5),
               sigmaM = round(sigmaM, 5)), 
       by = c('ratioBW' = 'ratioBW',
              'sigmaW' = 'TrueSigma',
              'sigmaM' = 'sigmaM',
              'nstud' = 'nstud',
              'BOLDC' = 'BOLDC')) %>% 
    # Bind to data frame
    bind_rows(EstVar_Vox,. )
    

  # # Now for each voxel the estimated between-study variance
  # EstVar_all <- MAvsIBMAres %>%
  #   # Filter estimated between-study variance components
  #   filter(parameter %in% c('IBMA.SE', 'ESTTAU')) %>%
  #   # Drop df
  #   dplyr::select(-FLAMEdf_3) %>%
  #   # Now add the true value of Cohen's d
  #   right_join(., 
  #      dplyr::select(TrueP_S, -ACT, -TrueCOPE, -TrueG), 
  #      by = c('sigmaW' = 'TrueSigma',
  #             'sigmaM' = 'sigmaM',
  #             'nstud' = 'nstud',
  #             'BOLDC' = 'BOLDC')) %>%
  #   # Rename the value to EstVar
  #   mutate(EstSE = value) %>% 
  #   dplyr::select(-value) %>%
  #   # Bind to data frame
  #   bind_rows(EstVar_all,. )
}




# Save intermediate results
saveRDS(coverage_Vox, file = paste(WIR, '/coverage_Vox.rda', sep = ''))
saveRDS(coverage_all, file = paste(WIR, '/coverage_all.rda', sep = ''))
saveRDS(CIlength_Vox, file = paste(WIR, '/CIlength_Vox.rda', sep = ''))
saveRDS(bias_Vox, file = paste(WIR, '/bias_Vox.rda', sep = ''))
saveRDS(bias_all, file = paste(WIR, '/bias_all.rda', sep = ''))
saveRDS(EstVar_Vox, file = paste(WIR, '/EstVar_Vox.rda', sep = ''))
#saveRDS(EstVar_all, file = paste(WIR, '/EstVar_all.rda', sep = ''))
















