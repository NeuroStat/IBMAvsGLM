####################
#### TITLE:     Find a voxel for illustration
#### Contents:
####
#### Source Files:
#### First Modified: 09/05/2019
#### Notes:
#################



##
###############
### Notes
###############
##

# Find a voxel to plot the CI

##
###############
### Preparation
###############
##

# Take argument from master file
input <- commandArgs(TRUE)
# Which machine
MACHINE <- try(as.character(input)[1], silent=TRUE)
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
  'Take[Ratio]' = try(as.character(input)[2],silent=TRUE),
  'Take[VariableN]' = try(as.character(input)[2],silent=TRUE),
  'Take[BetOvWith]' = try(as.character(input)[2],silent=TRUE))
}
if(MACHINE=='MAC'){
  DATAwd <- list(
    'Take[MAvsIBMA_Act]' = 
      "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/MAvsIBMA_act/Results_Parameters/Results",
    'Take[MAvsIBMA_Act_RanInSl]' = 
      "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/MAvsIBMA_act/Results_RanInSl/Results",
    'Take[GLMvsMA_wi_w_act]' = "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/GLMvsMA_wi_w_act/Results",
    'Take[HE]' = "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/HE/Results",
    'Take[Ratio]' = "/Volumes/Elements/IBMAvsGLM/Ratios/",
    'Take[VariableN]' = "/Volumes/Elements/IBMAvsGLM/VariableN",
    'Take[BetOvWith]' = "/Volumes/Elements/IBMAvsGLM/BetOvWith"
  )
}

NUMDATAwd <- length(DATAwd)
# Select your simulation setting here
currentWD <- 7

# Voxels possible
checkVox <- 134:729

# Which voxel (from 134:729)
chosenVox <- try(as.numeric(as.character(input)[3]), silent=TRUE)

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
info <- data.frame('Sim' = c(1, 2, 3, 4, 5, 6,7),
                   'nsim' = c(1000, 500, 1000, 1000, 1000, 1000, 1000),
                   'nsub' = c(trueMCvalues('sim_act', 'nsub'),
                            trueMCvalues('sim_act', 'nsub'),
                            trueMCvalues('sim_act', 'nsub'),
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
if(currentWD %in% c(3,4,5,6,7)){
  saveParam <- factor(levels = ParamsSim[1:12])
}

# Dimension of brain
DIM <- trueMCvalues('sim_act', 'DIM')

# Ratio of between- over within subject variability
ratioBW_vec <- c(0.25, 0.5, 0.75)

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
TrueParamDat$TrueSigma <- round(TrueParamDat$TrueSigma, 4)

# If preprocessing the take with null data, add extra part to TrueParamDat
if(currentWD %in% c(3,4,5,6,7)){
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
            'sigmaM' = round(sqrt(trueMCvalues('sim_act', 'TrueSigma2M')),4),
            'nstud' = trueMCvalues('sim_act', 'nstud')) %>%
    rowwise() %>%
    mutate(sigmaW = round(sqrt(trueMCvalues('sim_act', 'TrueSigma2W', 
                         ratioBW = ratioBW)[sigmaW_level]),4)) %>%
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
### Find the voxel
###############
##

# Number of simulations in the illustration
nsimIllu <- 100

# Comparison between DL and HE
ScenComs <- c('DL', 'HE')

# Empty data frame
illuCIL <- data.frame() %>% as_tibble()

# For loop over DL vs HE
for(s in 1:length(ScenComs)){
  # Select the scenario
  ScenCom <- ScenComs[s] 
  
  # Load in all the data and summarise 
  for(i in 1:nsimIllu){
    # Load data
    loadDat <- readRDS(
      paste(DATAwd[[currentWD]], '/', ScenCom, '/ActMAvsIBMA_',i, '.rda', 
            sep=''))
    
    # Process the loaded data
    rawDat <- loadDat %>%
      # Check one voxel
      filter(voxel == chosenVox) %>%
      # Filter the scenario we want
      filter(ratioBW == 0.75) %>%
      filter(., sigmaM >= 8) %>%
      filter(sigmaW >= 17 & sigmaW <= 30) %>%
      filter(BOLDC == 3) %>%
      filter(nstud == 50) %>%
      mutate(sigmaM = round(sigmaM, 4)) %>%
      mutate(sigmaW = round(sigmaW, 4)) 
    
    # Estimate without CI 
    Estimate <- 
      rawDat %>%
      # Also filter on one scenario
      filter(parameter == 'MA.WeightedAvg') %>%
      # Add the true values
      left_join(., TrueP_S, by = c('sigmaW' = 'TrueSigma',
                                   'ratioBW' = 'ratioBW',
                                   'sigmaM' = 'sigmaM',
                                   'nstud' = 'nstud',
                                   'BOLDC' = 'BOLDC'))
    # MA in wide format
    MAwide <- rawDat %>%
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
      left_join(
        filter(Estimate, parameter == 'MA.WeightedAvg') %>% 
          dplyr::select(-FLAMEdf_3),
        ., by = c("sim", "voxel", "SCEN", "ratioBW", "sigmaW", "sigmaB", 
                  "sigmaM", "nstud", "BOLDC")) %>%
      # Add indicator for the coverage
      mutate(cov_IND = ifelse(TrueG >= CI_lower & TrueG <= CI_upper,1, 0))
    
    # Add to data frame
    illuCIL <- bind_rows(illuCIL, MAwide)
  }
}

valDL <- 100 - illuCIL[illuCIL$SCEN == 'DL', 'cov_IND'] %>%
  unlist(.) %>% sum(.)
valHE <- 100 - illuCIL[illuCIL$SCEN == 'HE', 'cov_IND'] %>%
  unlist(.) %>% sum(.)
print(paste('In voxel ', chosenVox, ', we observe: DL = ', valDL, ' and HE = ', valHE, sep = ''))

if(valDL == 2 && valHE == 4){
  print(paste('Found a voxel! Check: ', chosenVox, sep = ''))
  cat(paste('Found a voxel! Check: ', chosenVox, sep = ''), 
      file = paste(DATAwd[[currentWD]], '/FoundVoxel.txt', sep = ''),
      append = TRUE)
}


