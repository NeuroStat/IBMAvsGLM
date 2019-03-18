####################
#### TITLE:     Plotting results of the MAvsIBMA_Act.R data.
#### Contents:
####
#### Source Files:
#### First Modified: 24/09/2017
#### Notes:
#################

# Created new R file for the MAvsGLM.R analysis.
# However, this file is based on same PreProcessGLMvsIBMA_Act.R file (just to confuse me).

##
###############
### Notes
###############
##

# Activation data!

# Blocked design for individual subjects.
# One condition.
# These N subjects are pooled using FLAME pooling
# The resulting images are converted to either Hedges' g and pooled using random effects meta-analysis.
# OR using 3e level GLM with again FLAME1.

# Data contains activation and between-study heterogeneity

# Set wd to base Github tree (line needs to be removed when making repo public)
setwd('/Users/hanbossier/Dropbox/PhD/PhDWork/Meta Analysis/R Code/Studie_Simulation/SimulationGit')

##
###############
### Functions
###############
##

# Libraries
library(ggplot2)
library(tibble)
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
library(scatterplot3d)
library(NeuRRoStat)
library(fMRIGI)


##
###############
### Preparation
###############
##

# Date of today
date <- Sys.Date()

# Set starting seed
set.seed(1990)

# Directories of the data for different simulations
DATAwd <- list(
  'Take[MAvsIBMA_Act]' = 
    "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/MAvsIBMA_act/Results_Parameters/Results",
  'Take[MAvsIBMA_Act_RanInSl]' = 
    "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/MAvsIBMA_act/Results_RanInSl/Results",
  'Take[GLMvsMA_wi_w_act]' = 
    "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/GLMvsMA_wi_w_act/Results"
)
NUMDATAwd <- length(DATAwd)
currentWD <- 3

# If available, load in Intermediate Results: this is location where summarized results are written
LIR <- list(
  'Take[MAvsIBMA_Act]' = 
    '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/MAvsIBMA_act/Results_Parameters/ProcessedResults',
  'Take[MAvsIBMA_Act_RanInSl]' = 
    "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/MAvsIBMA_act/Results_RanInSl/ProcessedResults",
  'Take[GLMvsMA_wi_w_act]' = 
    "~/Desktop/IBMA_tmp"
    #"/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/GLMvsMA_wi_w_act/Results/ProcessedResults"
)

# Number of conficence intervals
CIs <- c('MA-weightVar','GLM-t')
NumCI <- length(CIs)

# Data frame with number of simulations and subjects for current simulation
info <- data.frame('Sim' = c(1,2,3),
                   'nsim' = c(500, 500,1000),
                   'nsub' = rep(trueMCvalues('sim_act', 'nsub'),3))
nsim <- info[currentWD,'nsim']
nsub <- info[currentWD,'nsub']

# Parameters that will be checked
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

# Dimension of brain
DIM <- trueMCvalues('sim_act', 'DIM')


################
#### TRUE VALUES: load in R objects
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

# Data frame of true parameter values
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
  inner_join(.,ParamComb, by = c('TrueSigma')) %>%
  # Add true COPE value
  mutate(TrueCOPE = trueMCvalues('sim_act', 'BOLDC'))


##
###############
### Data Wrangling
###############
##

# Loading in raw data (TRUE) or processed data (FALSE)?
RAWDATA <- FALSE

if(RAWDATA){
# Subset of simulations
subset <- 10

#### To check some results, we have a subset of the simulations
for(i in 1:subset){
  MAvsIBMAres <- readRDS(
    paste(DATAwd[[currentWD]],'/ActMAvsIBMA_',i, '.rda', sep='')) %>%
      bind_rows(MAvsIBMAres,.)
}

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

# Data processing:
# 1) add CI_coverage and summarise
CoveragePlot <- ProcessedDat %>% filter(SmoothG != 0) %>%
  # 2) true smoothed value within CI limits?
  mutate(cov_IND = ifelse(parameter == "MA.WeightedAvg",
                          ifelse(SmoothG >= CI_lower & SmoothG <= CI_upper,1, 0),
                          ifelse(SmoothCOPE >= CI_lower & SmoothG <= CI_upper,1, 0))) %>%
  # 3) drop variables that we do not need to calculate coverages
  #       --> drop sigma as info is in TrueD
  select(sim, voxel, parameter, TrueD, tau, nstud, cov_IND) %>%
  group_by(sim, parameter, TrueD, tau, nstud) %>%
  summarise(coverage = mean(cov_IND),
            sdCov = sd(cov_IND)) 

}

# Also possible to process intermediate; summarized results
if(!RAWDATA){
  ### COVERAGE AVERAGED OVER ALL ACTIVE VOXELS ###
  # Values of each simulation for the covarage: usefull for plot with smoothing method
  SimCoverage <- readRDS(file = 
             paste(LIR[[currentWD]], '/coverage_Vox.rda', sep = ''))
  
  # Processed values (summarised) for coverage: ready to plot
  CoveragePlot <- readRDS(file = 
          paste(LIR[[currentWD]], '/coverage_Vox.rda', sep = '')) %>%
    # Rename coverage to wSimCoverage
    rename(wSimCoverage = coverage, wSimSDCov = sdCov) %>%
    # Drop wSIMSDCov (not interested in)
    dplyr::select(-wSimSDCov) %>%
    # Summarise over simulations
    group_by(parameter, TrueD, tau, nstud) %>%
    summarise(coverage = mean(wSimCoverage),
              sdCoverage = sd(wSimCoverage))

  ### COVERAGE AVERAGED OVER ALL SIMULATIONS THEN VOXELS ###  
  COVdata <- readRDS(file = 
                paste(LIR[[currentWD]], '/coverage_all.rda', sep = ''))
  # Summarise over simulations first
  COV <- COVdata %>% group_by(voxel, parameter, TrueD, tau, nstud) %>%
    summarise(AvgSimCov = mean(cov_IND)) %>%
    group_by(parameter, TrueD, tau, nstud) %>%
    summarise(AvgCOV = mean(AvgSimCov))
  
  ### COVERAGE NOT AVERAGED OVER ALL VOXELS ###
  COV_SIM <- COVdata %>% 
    # Average over all simulations
    group_by(voxel, parameter, TrueD, tau, nstud) %>%
    summarise(AvgCOV_sim = mean(cov_IND))

  # Values for CI length
  CILdata <- readRDS(file = 
              paste(LIR[[currentWD]], '/CIlength_all.rda', sep = ''))
  # Summarise over simulations
  CIL <- CILdata %>% group_by(voxel, parameter, TrueD, tau, nstud) %>%
    summarise(AvgSimCIL = mean(CIlength)) %>%
    # now summarise over voxels
    group_by(parameter, TrueD, tau, nstud) %>%
    summarise(AvgCIL = mean(AvgSimCIL))
  
  # Values for standardized bias
  BiasData <- readRDS(file = 
      paste(LIR[[currentWD]], '/bias_all.rda', sep = ''))
  # Summarise over simulations
  BIAS <- BiasData %>% group_by(voxel, parameter, TrueD, tau, nstud) %>%
    summarise(AvgBias = mean(bias),
           SDBias = sd(bias)) %>%
    mutate(StBias = AvgBias/SDBias * 100) %>%
    # Summarise over voxels
    group_by(parameter, TrueD, tau, nstud) %>%
    summarise(AvgStBias = mean(StBias))
  
  # Values for estimated variance
  EstVarData <- readRDS(file = 
            paste(LIR[[currentWD]], '/EstVar_all.rda', sep = ''))
  # Summarise over simulations
  EstVar <- EstVarData %>% 
    group_by(voxel, parameter, TrueD, tau, nstud) %>%
    summarise(AvgSimEstVar = mean(EstVar)) %>%
    ungroup() %>%
    # And then over voxels
    group_by(parameter, TrueD, tau, nstud) %>%
    summarise(AvgEstVar = mean(AvgSimEstVar))
  
}

#########################################################
###################### CI COVERAGE ######################
#########################################################

#### FIRST SECTION: AVERAGED OVER ALL VOXELS AND THEN SIMULATIONS ####
# Plot with SD bars. 
# However, these sd's are calculated on the average coverage of the grid of voxels.
# Per simulation. Hence the SD will be lower due to averaging coverage within each
# simulation. 
CoveragePlot %>%
  # For IBC, I switch from d to sigma again...
  left_join(., dplyr::select(TrueP_S, TrueD, TrueSigma), by = 'TrueD') %>% 
  # create labels for facets
  mutate(d = paste('d ~ "=" ~ ', TrueD, sep = ''),
         etaL = paste('sigma[b1]^2 ~ "=" ~ ', round(tau**2, 2), sep = ''),
         sigmaL = paste('sigma[W] ~ "=" ~ ', round(TrueSigma, 0), sep = '')) %>%
  mutate(sigmaLF = factor(sigmaL, levels = 
    paste('sigma[W] ~ "=" ~ ', 
          round(sqrt(trueMCvalues('sim_act', 'TrueSigma2W')), 0), sep = ''))) %>% 
  # 4) plot the results
  ggplot(., aes(x = nstud, y = coverage)) + 
  geom_point(aes(colour = parameter, fill = parameter), size = 1.15) +
  geom_line(aes(colour = parameter), size = 1.1) +
  geom_segment(aes(x = nstud, xend = nstud, y = coverage + sdCoverage,
                    yend = coverage - sdCoverage,
               colour = parameter), alpha = 0.75, linetype = 'dotted',
               size = .8) +
  scale_x_continuous('Number of studies in the third level') +
  scale_y_continuous('Empirical coverage') +
  scale_color_manual('Model', labels = c('random effects MA', 'mixed effects GLM'),
                     values = c('#fdb462', '#bc80bd')) +
  scale_fill_manual('Model', labels = c('random effects MA', 'mixed effects GLM'),
                     values = c('#fdb462', '#bc80bd')) +
  geom_hline(aes(yintercept = 0.95), colour = 'black') +
  ggtitle('empirical coverages of the 95% CIs') +
  facet_grid(sigmaLF ~ etaL, labeller = label_parsed) +
  theme_bw() +
  theme(legend.position = 'bottom',
        axis.text = element_text(face = 'bold', size = 9),
        strip.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(face = 'bold', size = 12),
        title = element_text(face = 'bold', size = 12),
        legend.text = element_text(face = 'bold'),
        plot.title = element_text(hjust = 0.5))
  

# Plot without the SD bars.
CoveragePlot %>%
  # create labels for facets
  mutate(d = paste('d ~ "=" ~ ', TrueD, sep = ''),
         etaL = paste('sigma[b1]^2 ~ "=" ~ ', round(eta**2, 2), sep = '')) %>%
  # 4) plot the results
  ggplot(., aes(x = nstud, y = coverage)) + 
  geom_point(aes(colour = parameter, fill = parameter), size = 0.8) +
  geom_line(aes(colour = parameter), size = 0.9) +
  scale_x_continuous('Number of studies in the third level') +
  scale_y_continuous('Empirical coverage') +
  scale_color_manual('Model', labels = c('Random effects MA', 'FLAME'),
                     values = c('#fdb462', '#bc80bd')) +
  scale_fill_manual('Model', labels = c('Random effects MA', 'FLAME'),
                    values = c('#fdb462', '#bc80bd')) +
  geom_hline(aes(yintercept = 0.95), colour = 'black') +
  ggtitle('Empirical coverages of 95% CI') +
  facet_grid(d ~ etaL, labeller = label_parsed) +
  theme_bw() 


# Plot without the SD bars: all eta on same panel
CoveragePlot %>%
  # First ungroup to round tau
  ungroup() %>%
  mutate(eta = round(eta,2)) %>%
  # Now group again
  group_by(parameter, TrueD, eta) %>%
  # create labels for facets
  mutate(d = paste('d ~ "=" ~ ', TrueD, sep = ''),
         etaL = paste('sigma[b1]^2 ~ "=" ~ ', round(eta**2, 2), sep = '')) %>%
  # 4) plot the results
  ggplot(., aes(x = nstud, y = coverage, group = interaction(parameter, factor(eta)))) + 
  geom_point(aes(colour = interaction(parameter, factor(eta)), 
                 fill = interaction(parameter, factor(eta))), 
             size = 0.8) +
  geom_line(aes(colour = interaction(parameter, factor(eta)), 
                group = interaction(parameter, factor(eta))), size = 0.7) +
  scale_x_continuous('Number of studies in the third level') +
  scale_y_continuous('Empirical coverage') +
  geom_hline(aes(yintercept = 0.95), colour = 'black') +
  ggtitle('Empirical coverages of 95% CI for each level of eta') +
  facet_grid(d ~ ., labeller = label_parsed) +
  theme_bw() 

# Check when we first average over simulations and then voxels
COV %>%
  # create labels for facets
  mutate(d = paste('d ~ "=" ~ ', TrueD, sep = ''),
         etaL = paste('sigma[b1]^2 ~ "=" ~ ', round(eta**2, 2), sep = '')) %>%
  # 4) plot the results
  ggplot(., aes(x = nstud, y = AvgCOV)) + 
  geom_point(aes(colour = parameter, fill = parameter), size = 1.15) +
  geom_line(aes(colour = parameter), size = 1.1) +
  scale_x_continuous('Number of studies in the third level') +
  scale_y_continuous('Empirical coverage') +
  scale_color_manual('Model', labels = c('Random effects MA', 'FLAME'),
                     values = c('#fdb462', '#bc80bd')) +
  scale_fill_manual('Model', labels = c('Random effects MA', 'FLAME'),
                    values = c('#fdb462', '#bc80bd')) +
  geom_hline(aes(yintercept = 0.95), colour = 'black') +
  ggtitle('Empirical coverages of 95% CI') +
  facet_grid(d ~ etaL, labeller = label_parsed) +
  theme_bw() +
  theme(legend.position = 'bottom')


# Violin plots without averaging over all voxels
COV_SIM %>% 
  ungroup() %>%
  filter(nstud %in% c(5,15,30,50)) %>%
  # For IBC, I switch from d to sigma again...
  left_join(., dplyr::select(TrueP_S, TrueD, TrueSigma), by = 'TrueD') %>% 
  # create labels for facets
  mutate(d = paste('d ~ "=" ~ ', TrueD, sep = ''),
         etaL = paste('sigma[b1]^2 ~ "=" ~ ', round(eta**2, 2), sep = ''),
         sigmaL = paste('sigma[W] ~ "=" ~ ', round(TrueSigma, 0), sep = '')) %>%
  mutate(sigmaLF = factor(sigmaL, levels = 
    paste('sigma[W] ~ "=" ~ ', 
          round(sqrt(trueMCvalues('sim_act', 'TrueSigma2W')), 0), sep = ''))) %>% 
  # 4) plot the results
  ggplot(., aes(x = factor(nstud), y = AvgCOV_sim, fill = parameter)) + 
  # Use violin plot
  geom_violin(position = 'dodge', alpha = .75, aes(factor(nstud))) +
  facet_grid(sigmaLF ~ etaL, labeller = label_parsed) +
  # Have mean line on top of it
  stat_summary(fun.y=mean, geom="line", aes(group = parameter), size = .7)  +
  scale_x_discrete('Number of studies in the third level') +
  scale_y_continuous('Empirical coverage') +
  scale_color_manual('Model', labels = c('random effects MA', 'mixed effects GLM'),
                     values = c('#fdb462', '#bc80bd')) +
  scale_fill_manual('Model', labels = c('random effects MA', 'mixed effects GLM'),
                    values = c('#fdb462', '#bc80bd')) +
  geom_hline(aes(yintercept = 0.95), colour = 'black', linetype = 'dashed') +
  ggtitle('empirical coverages of the 95% CIs') +
  theme_bw() +
  theme(legend.position = 'bottom',
        axis.text = element_text(face = 'bold', size = 9),
        strip.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(face = 'bold', size = 12),
        title = element_text(face = 'bold', size = 12),
        legend.text = element_text(face = 'bold'),
        plot.title = element_text(hjust = 0.5))


#########################################################
####################### CI LENGTH #######################
#########################################################


# Plot without the SD bars.
CIL %>%
  # create labels for facets
  mutate(d = paste('d ~ "=" ~ ', TrueD, sep = ''),
         etaL = paste('eta ~ "=" ~ ', round(eta, 2), sep = '')) %>%
  # 4) plot the results
  ggplot(., aes(x = nstud, y = AvgCIL)) + 
  geom_point(aes(colour = parameter, fill = parameter), size = 0.8) +
  geom_line(aes(colour = parameter), size = 0.9) +
  scale_x_continuous('Number of studies in the third level') +
  scale_y_continuous('Average 95% CI length') +
  scale_color_manual('Model', labels = c('Random effects MA', 'FLAME'),
                     values = c('#fdb462', '#bc80bd')) +
  scale_fill_manual('Model', labels = c('Random effects MA', 'FLAME'),
                    values = c('#fdb462', '#bc80bd')) +
  ggtitle('Average length of the 95% CI') +
  facet_grid(d ~ etaL, labeller = label_parsed) +
  theme_bw() 


#########################################################
################### STANDARDIZED BIAS ###################
#########################################################


# Plot without the SD bars.
BIAS %>%
  # For IBC, I switch from d to sigma again...
  left_join(., dplyr::select(TrueP_S, TrueD, TrueSigma), by = 'TrueD') %>% 
  # For OHBM, I convert tau to eta
  mutate(eta = ifelse(tau == 0, 0,
                      ifelse(tau == 10, 0.1, 0.4))) %>%
  # create labels for facets
  mutate(d = paste('d ~ "=" ~ ', TrueD, sep = ''),
         etaL = paste('sigma[b1]^2 ~ "=" ~ ', round(eta**2, 2), sep = ''),
         sigmaL = paste('sigma[W] ~ "=" ~ ', round(TrueSigma, 0), sep = '')) %>%
  mutate(sigmaLF = factor(sigmaL, levels = 
                  paste('sigma[W] ~ "=" ~ ', 
        round(sqrt(trueMCvalues('sim_act', 'TrueSigma2W')), 0), sep = ''))) %>% 
  # # remove wrong results of MA
  # mutate(AvgStBias = ifelse(parameter == 'MA.WeightedAvg' & TrueD != 0.14,
  #                           NA, AvgStBias)) %>%
  # 4) plot the results
  ggplot(., aes(x = nstud, y = AvgStBias)) + 
  geom_point(aes(colour = parameter, fill = parameter), size = 0.8) +
  geom_line(aes(colour = parameter), size = 0.9) +
  scale_x_continuous('Number of studies in the third level') +
  scale_y_continuous('Standardized bias') +
  scale_color_manual('Model', labels = c('random effects MA', 'mixed effects GLM'),
                     values = c('#fdb462', '#bc80bd')) +
  scale_fill_manual('Model', labels = c('random effects MA', 'mixed effects GLM'),
                    values = c('#fdb462', '#bc80bd')) +
  ggtitle('average standardised bias') +
  facet_grid(sigmaLF ~ etaL, labeller = label_parsed) +
  theme_bw() +
  theme(legend.position = 'bottom',
        axis.text = element_text(face = 'bold', size = 9),
        strip.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(face = 'bold', size = 12),
        title = element_text(face = 'bold', size = 12),
        legend.text = element_text(face = 'bold'),
        plot.title = element_text(hjust = 0.5))


##########################################################
################### ESTIMATED VARIANCE ###################
##########################################################


# Estimated variance components
EstVar %>%
  # create labels for facets
  mutate(d = paste('d ~ "=" ~ ', TrueD, sep = ''),
         etaL = paste('eta ~ "=" ~ ', round(eta, 2), sep = '')) %>%
  ungroup() %>%
  group_by(parameter, TrueD, eta) %>% 
  # 4) plot the results
  ggplot(., aes(x = nstud, y = AvgEstVar)) + 
  geom_point(aes(colour = eta, fill = eta), size = 0.8) +
  geom_line(aes(colour = eta), size = 0.9) +
  scale_x_continuous('Number of studies in the third level') +
  scale_y_continuous('Estimated Variance') +
  # scale_color_manual('Model', labels = c('Random effects MA', 'FLAME'),
  #                    values = c('#fdb462', '#bc80bd')) +
  # scale_fill_manual('Model', labels = c('Random effects MA', 'FLAME'),
  #                   values = c('#fdb462', '#bc80bd')) +
  ggtitle('Average Estimated Variance') +
  facet_wrap(d ~ parameter, labeller = label_parsed, scales = 'free_y', nrow = 3) +
  theme_bw() 

##
###############
### Plotting: coverage and weighted averages
###############
##

# BETTER PLOT 3D
levelplot()

# Function for levelplot
ValuesOnLevelPlot2D <- function(x, y, z, ...) {
  panel.levelplot(x,y,z,...)
  panel.text(x, y, round(z,3),col='red')
}

# Plotting the coverages
levelplot(array(unlist(mean.coverages[['MA']]), dim = DIM))
levelplot(array(unlist(mean.coverages[['IBMA']]), dim=rep(sqrt(prod(DIM)),2)))


CCI1 <- levelplot(array(unlist(mean.coverages[['MA']]), dim= DIM),
                  col.regions = gray(0:100/100), at=seq(0,1,by=0.05), main='Meta-Analysis',xlab='',ylab='',
                  colorkey = list(space = "bottom"))
CCI2 <- levelplot(array(unlist(mean.coverages[['IBMA']]), dim= DIM),
                  col.regions = gray(0:100/100), at=seq(0,1,by=0.05), main='3 level GLM',xlab='',ylab='',
                  colorkey = list(space = "bottom"))
grid.arrange(CCI1,CCI2,nrow=1,top = textGrob('CI - Coverage of each voxel over 1000 simulations.', gp=gpar(fontsize=20,font=1)))


# Try with ggplot
plotCI <- data.frame(CI = rbind(mean.coverages[['MA']],mean.coverages[['IBMA']]),
                     Type = rep(c("MA", "mixed effects GLM"), each = 64),
                     x = factor(rep(rep(1:8,8),2)),
                     y = factor(rep(rep(1:8, each = 8),2)))
ggplot(plotCI, aes(x=x,y=y)) + geom_tile(aes(fill = CI)) +
  facet_wrap(~Type) +
  scale_fill_gradient2(name = "", low = scales::muted("blue"), mid = "white",
                       high = scales::muted("red"), midpoint = 0.95, space = "Lab",
                       na.value = "grey50", guide = "colourbar") +
  scale_x_discrete(name = "") + scale_y_discrete(name = "") + ggtitle("") +
  theme_minimal() +
  theme(strip.text = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.9,"cm"))



# Plotting the lengths
LCI1 <- levelplot(array(mean.lengths[['MA']], dim=rep(sqrt(prod(DIM)),2)),
                  col.regions = gray(100:0/100), at=seq(0,1,by=0.05), main='Meta-Analysis',xlab='',ylab='',
                  colorkey = list(space = "bottom"),
                  panel=ValuesOnLevelPlot2D)
LCI2 <- levelplot(array(mean.lengths[['IBMA']], dim=rep(sqrt(prod(DIM)),2)),
                  col.regions = gray(100:0/100), at=seq(0,1,by=0.05), main='3 level GLM',xlab='',ylab='',
                  colorkey = list(space = "bottom"),
                  panel=ValuesOnLevelPlot2D)
grid.arrange(LCI1,LCI2,nrow=1,top = textGrob('CI - Length of each voxel over 3000 simulations.', gp=gpar(fontsize=20,font=1)))


# Plotting first the distribution of either weighted average, or COPE-value
WA_density <- ggplot(data.frame('value' = MA.MEANBETA), aes(x=value)) + geom_density(fill='#1b9e77')
COPE_density <- ggplot(data.frame('value' = IBMA.MEANBETA), aes(x=value)) + geom_density(fill='#7570b3')
grid.arrange(WA_density, COPE_density, nrow=1, top = textGrob('Average density (over all simulations) of weighted average/COPE', gp=gpar(fontsize=20,font=1)))

# Plotting the standardized bias
BCI1 <- levelplot(array(mean.bias[['MA']], dim=rep(sqrt(prod(DIM)),2)),
                  col.regions = gray(100:0/100), at=seq(0,ceiling(max(unlist(mean.bias))),length.out=100), main='Meta-Analysis',xlab='',ylab='',
                  colorkey = list(space = "bottom"),
                  panel=ValuesOnLevelPlot2D)
BCI2 <- levelplot(array(mean.bias[['IBMA']], dim=rep(sqrt(prod(DIM)),2)),
                  col.regions = gray(100:0/100), at=seq(0,ceiling(max(unlist(mean.bias))),length.out=100), main='3 level GLM',xlab='',ylab='',
                  colorkey = list(space = "bottom"),
                  panel=ValuesOnLevelPlot2D)
grid.arrange(BCI1,BCI2,nrow=1, top = textGrob('Standardized bias (%) of each voxel over 3000 simulations.', gp=gpar(fontsize=20,font=1)))

# Try with ggplot
PlotBias <- data.frame(bias = rbind(mean.bias[['MA']],mean.bias[['IBMA']]),
                       Type = rep(c("MA", "mixed effects GLM"), each = 64),
                       x = factor(rep(rep(1:8,8),2)),
                       y = factor(rep(rep(1:8, each = 8),2)))
ggplot(PlotBias, aes(x=x,y=y)) + geom_tile(aes(fill = bias)) +
  facet_wrap(~Type) +
  scale_fill_gradient2(name = "", low = scales::muted("blue"), mid = "white",
                       high = scales::muted("red"), midpoint = 0, space = "Lab",
                       na.value = "grey50", guide = "colourbar") +
  scale_x_discrete(name = "") + scale_y_discrete(name = "") + ggtitle("Average standardized bias") +
  theme_minimal()

# Maybe try the difference between the two, alongside table of average and SD
DiffBias <- data.frame(Diffbias = c(mean.bias[['MA']] - mean.bias[['IBMA']]),
                       x = factor(rep(1:8,8)),
                       y = factor(rep(1:8, each = 8)))
DiffBiasPlot <- ggplot(DiffBias, aes(x=x,y=y)) + geom_tile(aes(fill = Diffbias)) +
  scale_fill_gradient2(name = "", low = scales::muted("blue"), mid = "white",
                       high = scales::muted("red"), midpoint = 0, space = "Lab",
                       na.value = "grey50", guide = "colourbar") +
  scale_x_discrete(name = "") + scale_y_discrete(name = "") + ggtitle("") +       # Standardized bias (MA) - standardized bias (GLM)
  theme_minimal() +
  theme(strip.text = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.9,"cm"))


tableBias <- data.frame(Type = CI.bias$CI, average = CI.bias$Mean, SD = CI.bias$SD)
tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
tbl <- tableGrob(tableBias, rows=NULL, theme=tt)
# Plot chart and table into one object
grid.arrange(tbl , DiffBiasPlot,
             nrow=1,
             as.table=TRUE,
             widths=c(1,2))





##
###############
### Calculate a CI around the second level cope-map
###############
##

lvl2.CI.upper.t <- lvl2.CI.lower.t <- array(NA, dim=c(prod(DIM)*nstud,30))

# This means we go into a simulation and look at the K studies.
# We can then construct a CI around each of these studies

# Start with 30 simulations
for(i in 1:30){
  for(k in 1:nstud){
    INDEX <- (k-1) * prod(DIM)
    START <- INDEX + 1
    END <- INDEX + prod(DIM)
    cope <- readNIfTI(paste(DATAwd[[currentWD]],'/',i,'/SCEN_1/study',k,'_stats/cope1.nii',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
    varcope <- readNIfTI(paste(DATAwd[[currentWD]],'/',i,'/SCEN_1/study',k,'_stats/varcope1.nii',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]

    lvl2.CI.upper.t[START:END,i] <- matrix(cope,ncol=1) + (qt(0.975, df=c((nsub/nstud)-1)) * sqrt(matrix(varcope,ncol=1)))
    lvl2.CI.lower.t[START:END,i] <- matrix(cope,ncol=1) - (qt(0.975, df=c((nsub/nstud)-1)) * sqrt(matrix(varcope,ncol=1)))}
}






##
###############
### Comparing maps on level 2
###############
##

# I want to compare COPE maps with hedges' g on second level. Then the varcope and variance of hedges's
# Start with COPE vs Hedges'g of the first study
copeStud <- hedgeStud <- array(NA, dim=c(prod(DIM),nsim))
varcopeStud <- stweights <- array(NA, dim=c(prod(DIM),nsim))
for(i in 1:nsim){
  if(i==c(nsim/2)) print('At 50%')
  copeStud[,i] <- readNIfTI(paste(DATAwd[[currentWD]],'/',i,'/SCEN_1/study1_stats/cope1.nii', sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
  hedgeAllStud <- AllData[which(OBJ.ID=='STHEDGE'),i]
  hedgeStud[,i] <- hedgeAllStud[c(1:prod(DIM))]

  varcopeStud[,i] <- readNIfTI(paste(DATAwd[[currentWD]],'/',i,'/SCEN_1/study1_stats/varcope1.nii', sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
  stweightsAllStud <- AllData[which(OBJ.ID=='STWEIGHTS'),i]
  stweights[,i] <- stweightsAllStud[c(1:prod(DIM))]
}

str(stweights)

# Compare COPE with hedges g
# Start with comparing over all values
par(mfrow=c(1,2), oma=c(0,0,2,0))
qqplot(x=copeStud, y= hedgeStud, xlab='Level 2 in GLM - all values', ylab='Meta-analysis approach')
# Averaged over all simulations.
qqplot(x=apply(copeStud,1,mean,na.rm=TRUE), y=apply(hedgeStud,1,mean,na.rm=TRUE),xlab='Level 2 in GLM - averaged', ylab='Meta-analysis approach')
title("Q-Q plot comparing 2e level COPE and hedges g", outer=TRUE)
par(mfrow=c(1,1))

# Histogram
frameLvl2 <- data.frame('value' = matrix(c(hedgeStud,copeStud),ncol=1), 'source' = rep(c('g', 'COPE'), each = prod(dim(copeStud))))
ggplot(frameLvl2, aes(x=value)) + geom_density(aes(colour=source),size=2)

# Compare VARCOPE with variance of Weighted average
# Over all values
par(mfrow=c(1,2), oma=c(0,0,2,0))
qqplot(x=varcopeStud, y= c(1/stweights), xlab='Level 2 in GLM', ylab='Meta-analysis approach')
# Averaged over all simulations
qqplot(x=apply(varcopeStud,1,mean,na.rm=TRUE), y=c(1/apply(stweights, 1, mean, na.rm=TRUE)),
       xlab='Level 2 in GLM', ylab = 'Meta-analysis approach')
title("Q-Q plot comparing 2e level VARCOPE and variance of g in the first study", outer=TRUE)
par(mfrow=c(1,1))

# Histogram
frameLvl2 <- data.frame('value' = matrix(c(as.numeric(1/stweights),as.numeric(varcopeStud)),ncol=1), 'source' = rep(c('variance g', 'VARCOPE'), each = prod(dim(copeStud))))
ggplot(frameLvl2, aes(x=value)) + geom_histogram(aes(fill=source)) +
  theme(legend.position='top')



##
###############
### Comparing maps on level 3
###############
##

# We will go to all the cope and varcopes of each study
copeMA <- array(NA, dim=c(prod(DIM),nsim))
varcopeMA <- array(NA, dim=c(prod(DIM),nsim))
for(i in 1:nsim){
  if(i==c(nsim/2)) print('At 50%')
  copeMA[,i] <- readNIfTI(paste(DATAwd[[currentWD]],'/',i,'/SCEN_1/MA_stats/cope1.nii', sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
  varcopeMA[,i] <- readNIfTI(paste(DATAwd[[currentWD]],'/',i,'/SCEN_1/MA_stats/varcope1.nii', sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
}

# Compare COPE with Weighted average
# Start with comparing over all values
par(mfrow=c(1,2), oma=c(0,0,2,0))
qqplot(x=copeMA, y= AllData[which(OBJ.ID=='MA.WeightedAvg'),], xlab='Level 3 in GLM - all values', ylab='Meta-analysis approach')
# Averaged over all simulations (which we already did in previous part).
qqplot(x=apply(copeMA,1,mean,na.rm=TRUE), y=MA.MEANBETA,xlab='Level 3 in GLM - averaged over simulations', ylab = 'Meta-analysis approach')
title("Q-Q plot comparing 3e level COPE and weighted average in the MA", outer=TRUE)

# Histogram
frameLvl3 <- data.frame('value' = matrix(c(AllData[which(OBJ.ID=='MA.WeightedAvg'),],copeMA),ncol=1), 'source' = rep(c('WAvg', 'COPE3'), each = prod(dim(copeMA))))
ggplot(frameLvl3, aes(x=value)) + geom_density(aes(colour=source),size=2) +
  theme(legend.position='top')

# Compare VARCOPE with variance of Weighted average
# Over all values
par(mfrow=c(1,2), oma=c(0,0,2,0))
qqplot(x=varcopeMA, y= AllData[which(OBJ.ID=='CI.MA.weightedVariance'),], xlab='Level 3 in GLM - all values', ylab='Meta-analysis approach')
# Averaged over all simulations
qqplot(x=apply(varcopeMA,1,mean,na.rm=TRUE), y=apply(AllData[which(OBJ.ID=='CI.MA.weightedVariance'),], 1, mean, na.rm=TRUE),
       xlab='Level 3 in GLM - averaged over simulations', ylab = 'Meta-analysis approach')
title("Q-Q plot comparing 3e level VARCOPE and weighted variance in MA", outer=TRUE)

# Histogram
frameLvl3 <- data.frame('value' = matrix(c(AllData[which(OBJ.ID=='CI.MA.weightedVariance'),],varcopeMA),ncol=1), 'source' = rep(c('WVar', 'VARCOPE3'), each = prod(dim(copeMA))))
ggplot(frameLvl3, aes(x=value)) + geom_density(aes(colour=source),size=2) +
  theme(legend.position='top')






