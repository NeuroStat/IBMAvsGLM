####################
#### TITLE:     Plotting results of the MAvsIBMA_Act.R data.
#### Contents:
####
#### Source Files:
#### First Modified: 24/09/2017
#### Notes:
#################


##
###############
### Notes
###############
##

# Plotting the analyses of MA vs GLM approach. 
# 3-level simulation approach.

# Blocked design for individual subjects.
# One condition.
# These N subjects are pooled using FLAME pooling
# The resulting images are converted to either Hedges' g and pooled using random effects meta-analysis.
# OR using 3e level GLM with again FLAME1.

# Data contains between-study heterogeneity!

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
library(cowplot)
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

# Function to try to know what the true value of between-study variability is
GetTrueSigmaM <- function(sigmaW, nsub = 29, nstud, sigmaM){
  # The true values
  sigma2W <- fMRIGI::trueMCvalues('sim_act', 'TrueSigma2W')
  sigma2B <- fMRIGI::trueMCvalues('sim_act', 'TrueSigma2B')
  sigma2M <- fMRIGI::trueMCvalues('sim_act', 'TrueSigma2M')
  I2 <- fMRIGI::trueMCvalues('sim_act', 'I2')
  design_factor_lvl1 <- fMRIGI::trueMCvalues('sim_act', 'design_factor')
  design_factor_lvl2 <- 1/nsub
  design_factor_lvl3 <- 1/nstud
  
  # The selected value of the sigma2W vector, to take the corresponding value
  # of the sigma2B vector
  ID_sigmaW <- which(round(sqrt(sigma2W),2) == round(sigmaW,2), arr.ind = TRUE)
  
  # The selected value of the sigma2M vector, to take the corresponding value
  # of the I2 vector
  ID_sigmaM <- which(round(sqrt(sigma2M),2) == round(sigmaM,2), arr.ind = TRUE)
  
  # The true value of between-study variability due to the esimtation in three stages
  VarB <- (sigma2B[ID_sigmaW] + sigma2W[ID_sigmaW]*design_factor_lvl1) * design_factor_lvl2
  VarM <- sigma2M[ID_sigmaM] + VarB
  TrueVarCope <- VarM/nstud
  
  return(TrueVarCope)
}


##
###############
### Preparation
###############
##

# Date of today
date <- Sys.Date()

# Set starting seed
set.seed(1990)

# Source paths
source('blind_MAvsGLM.R')

# Which scenario do we analyze?
# Go for 4: variable N and true value for MA = G.
currentWD <- 5

# Number of conficence intervals
CIs <- c('MA-weightVar','GLM-t')
NumCI <- length(CIs)

# Data frame with number of simulations and subjects for current simulation
info <- data.frame('Sim' = c(1,2,3,4,5),
                   'nsim' = c(1000, 1000, 1000, 1000, 1000),
                   'nsub' = rep(trueMCvalues('sim_act', 'nsub'),5))
nsim <- info[currentWD,'nsim']
nsub <- info[currentWD,'nsub']

# Parameters that will be checked
saveParam <- factor(levels = c('CI.MA.upper.weightVar', 'CI.MA.lower.weightVar',
                               'MA.WeightedAvg',
                               'CI.IBMA.upper.t','CI.IBMA.lower.t', 'IBMA.COPE',
                               'CI.MA.weightedVariance', 'STHEDGE', 'ESTTAU',
                               'STWEIGHTS', 'STWEIGHTS_ran', 'IBMA.SE'))

# Vector of scenarios
SCEN <- c('GLM', 'DL', 'HE', 'REML')
  

# Data frame with simulation results:
MAvsIBMAres <- tibble(sim = integer(),
                      voxel = integer(),
                      value = numeric(),
                      parameter = saveParam,
                      SCEN = factor(levels = SCEN),
                      BOLDC = numeric(),
                      ratioBW = numeric(),
                      sigmaW = numeric(),
                      sigmaB = numeric(),
                      sigmaM = numeric(),
                      nstud = numeric(),
                      FLAMEdf_3 = numeric())

# Empty data frames
CoveragePlot <- CoverageViolinPlot <- CIL <- BIAS <-
  BIASviolin <- EstVar <- as_tibble()

# Dimension of brain
DIM <- trueMCvalues('sim_act', 'DIM')

# Ratio of between- over within subject variability
ratioBW_vec <- c(0.25, 0.5, 0.75)

# Want results of each voxel?
AllVox <- FALSE

# Colours
colours <- data.frame(scen = SCEN,
                      cols = c("#66c2a5",
                               "#984ea3",
                               "#8da0cb",
                               "#e78ac3"),
                      stringsAsFactors = FALSE)

# Comparisons in the graphs
comps <- list(c('GLM', 'DL'), c('GLM', 'HE'), c('GLM', 'REML'))

################
#### TRUE VALUES: load in R objects
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

# Data frame with link between sigma and Cohen's d, will add this to the
# parameter combinations later on.
# NOTE: don't add sigmaM here, as there is no link between d and sigmaM!!
TrueParamDat <- data.frame(Nsub = trueMCvalues('sim_act', 'nsub'),
                           TrueD = rep(trueMCvalues('sim_act', 'TrueD'), 2),
                           TrueSigma = rep(sqrt(trueMCvalues('sim_act', 'TrueSigma2W')), 2),
                           TrueG = rep(trueMCvalues('sim_act', 'TrueG'), 2),
                           BOLDC = rep(trueMCvalues('sim_act', 'BOLDC'), each = 3)) %>%
  mutate(ACT = ifelse(BOLDC > 0, 1, 0)) %>%
  mutate(TrueD = TrueD * ACT,
         TrueG = TrueG * ACT)

# Data frame with combinations of all simulations run: note that we expand the 
# grid with the act/no activation in the TrueParamDat data frame.
# Hence we do not need to add BOLDC here to the data frame!
ParamComb <- expand.grid(
  'sigmaW' = sqrt(trueMCvalues('sim_act', 'TrueSigma2W')),
  'sigmaM' = sqrt(trueMCvalues('sim_act', 'TrueSigma2M')),
  'nstud' = trueMCvalues('sim_act', 'nstud'))
NumPar <- dim(ParamComb)[1]

# Extend the true values with number of studies
TrueP_S <- TrueParamDat %>% 
  full_join(.,ParamComb,
            by = c('TrueSigma' = 'sigmaW')) %>% 
  # Add true COPE value (which is the same as BOLDC actually)
  mutate(TrueCOPE = ACT * trueMCvalues('sim_act', 'BOLDC')[2]) %>%
  as_tibble()


##
###############
### Data Wrangling
###############
##

# Load in intermediate; summarized results
### COVERAGE AVERAGED OVER ALL VOXELS ###
# Processed values (averaged over voxels) of coverage
for(s in 1:length(SCEN)){
  CoveragePlot <- readRDS(file = 
      paste(LIR[[currentWD]], '/', SCEN[s], '/coverage_Vox.rda', sep = '')) %>%
    # Rename coverage to wSimCoverage
    rename(wSimCoverage = coverage, wSimSDCov = sdCov) %>%
    # Drop wSIMSDCov (not interested in)
    dplyr::select(-wSimSDCov) %>%
    # Summarise over simulations
    group_by(parameter, SCEN, TrueD, ratioBW, sigmaW, sigmaM, nstud) %>%
    summarise(coverage = mean(wSimCoverage, na.rm = TRUE),
              sdCoverage = sd(wSimCoverage, na.rm = TRUE)) %>%
    bind_rows(CoveragePlot, .)

  ### COVERAGE AVERAGED OVER ALL SIMULATIONS ###
  if(AllVox){
    CoverageViolinPlot <- readRDS(file = 
      paste(LIR[[currentWD]], '/', SCEN[s], '/coverage_all.rda', sep = '')) %>%
      group_by(voxel, parameter, SCEN, TrueD, ratioBW, sigmaW, sigmaM, nstud) %>%
      summarise(coverage = mean(cov_IND, na.rm = TRUE)) %>%
      bind_rows(CoverageViolinPlot, .)
  }

  ### CI LENGHTS AVERAGED OVER ALL VOXELS ###
  CILdata <- readRDS(file = 
    paste(LIR[[currentWD]], '/', SCEN[s], '/CIlength_Vox.rda', sep = ''))
  # Summarise over simulations
  CIL <- CILdata %>% group_by(parameter, SCEN, TrueD, ratioBW, sigmaW, sigmaM, nstud) %>%
    summarise(AvgSimCIL = mean(AvCIlength, na.rm = TRUE)) %>%
    bind_rows(CIL, .)
  
  ### STANDARDIZED BIAS AVERAGED OVER ONE VOXELS ###
  BiasData <- readRDS(file = 
    paste(LIR[[currentWD]], '/', SCEN[s], '/bias_Vox.rda', sep = ''))
  # Summarise over simulations
  BIAS <- BiasData %>% 
    dplyr::select(-sigmaB) %>%
    group_by(parameter, SCEN, TrueD, ratioBW, sigmaW, sigmaM, nstud) %>%
    summarise(AvgBias = mean(Avbias, na.rm = TRUE),
              SDBias = sd(Avbias, na.rm = TRUE)) %>%
    mutate(StBias = (AvgBias/SDBias)) %>%
    bind_rows(BIAS, .)
  
  ### STANDARDIZED BIAS AVERAGED OVER ALL SIMULATIONS ###
  if(AllVox){
    BiasData_all <- readRDS(file = 
        paste(LIR[[currentWD]], '/', SCEN[s], '/bias_all.rda', sep = ''))
    # Summarise over simulations
    BIASviolin <- BiasData_all %>% 
      group_by(voxel, parameter, SCEN, TrueD, ratioBW, sigmaW, sigmaM, nstud) %>%
      summarise(AvgBias = mean(bias, na.rm = TRUE),
                SDBias = sd(bias, na.rm = TRUE)) %>%
      mutate(StBias = AvgBias/SDBias * 100) %>%
      bind_rows(BIASviolin, .)
  }

  ### ESTIMATED BETWEEN-STUDY VARIANCE AVERAGED OVER ALL VOXELS ###
  EstVarData <- readRDS(file = 
    paste(LIR[[currentWD]], '/', SCEN[s], '/EstVar_Vox.rda', sep = ''))
  # Summarise over simulations
  EstVar <- EstVarData %>% 
    group_by(parameter, SCEN, TrueD, ratioBW, sigmaW, sigmaM, nstud) %>%
    summarise(AvgSE = mean(AvgEstSE, na.rm = TRUE)) %>%
    bind_rows(EstVar, .)
}


#########################################################
###################### CI COVERAGE ######################
#########################################################
CoveragePlot %>%
  filter(SCEN=='DL') %>%
  filter(ratioBW==0.75) %>%
  filter(nstud == 50) %>%
  filter(sigmaW <= 14)

#### FIRST SECTION: AVERAGED OVER ALL VOXELS AND THEN SIMULATIONS ####
CoveragePlot <- CoveragePlot %>%
  # Column for activation YES/NO and turn into factor
  mutate(signal = ifelse(TrueD == 0, 'null', 'activation')) %>%
  mutate(signalF = factor(signal, levels = c('null', 'activation'),
                          labels = c('null data', '3% BOLD signal change'))) %>%
  # create labels for facets
  mutate(d = paste('d ~ "=" ~ ', TrueD, sep = ''),
         sigmaWL = paste('sigma[S] ~ "=" ~ ', round(sigmaW, 0), sep = ''),
         sigmaML = paste('sigma[M] ~ "=" ~ ', round(sigmaM, 0), sep = ''),
         ratioL = paste('sigma[G]/sigma[S] ~ "=" ~ ', ratioBW))

# Find the minimal Y-axis
Y_min <- min(CoveragePlot$coverage)

# Loop over the comparisons 
for(s in 1:length(comps)){
  # Loop over the ratioBWs
  for(r in 1:length(ratioBW_vec)){
    # First create variable
    LoopPlot <- 
      CoveragePlot %>%
      mutate(sigmaWLF = factor(sigmaWL, levels = unique(sigmaWL))) %>%
      mutate(sigmaMLF = factor(sigmaML, levels = unique(sigmaML))) %>%
      # Select the ratio
      filter(ratioBW == ratioBW_vec[r]) %>% 
      # Select the comparison
      filter(SCEN %in% comps[[s]]) %>% 
      # 4) plot the results
      ggplot(., aes(x = nstud, y = coverage)) +
      geom_line(aes(colour = parameter, linetype = signalF), size = 0.9) +
      geom_point(aes(colour = parameter, fill = parameter), size = 0.95, show.legend = FALSE) +
      scale_x_continuous(ifelse(r == 2,
        'Number of studies in the third level', '')) +
      scale_y_continuous(ifelse(r == 1, 'Empirical coverage', ''),
                         # I'm truncating the Y-axis otherwise you cannot see
                         # differences between the methods...
                         limits = c(0.8, 1)) +
                         #limits = c(Y_min, 1)) +
      scale_color_manual('Model',
        labels = c(paste0('random effects MA: ', comps[[s]][2]),
                                             'mixed effects GLM'),
                         values = c(colours[colours$scen == comps[[s]][2], 'cols'],
                                    colours[colours$scen == 'GLM', 'cols'])) +
      scale_fill_manual('Model',
                        labels = c(paste0('random effects MA: ', comps[[s]][2]),
                                            'mixed effects GLM'),
                        values = c(colours[colours$scen == comps[[s]][2], 'cols'],
                                   colours[colours$scen == 'GLM', 'cols'])) +
      scale_linetype_manual('', values = c('dashed', 'solid')) +
      geom_hline(aes(yintercept = 0.95), colour = 'black') +
      # ggtitle(bquote(sigma[B]/sigma[W] == .(ratioBW_vec[r]))) +
      ggtitle(paste('R = ', ratioBW_vec[r], sep = '')) +
      # ggtitle(paste0('Ratio: ', ratioBW_vec[r])) +
      facet_grid(sigmaWLF ~ sigmaMLF, labeller = label_parsed) +
      theme_bw() +
      theme(legend.position = 'none',
            axis.text = element_text(face = 'bold', size = 9, vjust = -1),
            strip.background = element_rect(fill = 'white', colour = 'white'),
            strip.text = element_text(face = 'bold', size = 12),
            title = element_text(face = 'bold', size = 12),
            legend.text = element_text(face = 'bold', size = 12),
            plot.title = element_text(hjust = 0.5))
  
    # Print the plot (otherwise it does not store the colours!)
    print(LoopPlot)
    
    # Wait for a second (otherwise objects are not created)
    Sys.sleep(2)
    
    # Now assign to variable
    assign(x = paste0('plot_', comps[[s]][2], '_ratio_', ratioBW_vec[r]), LoopPlot)
    
    # Reset
    rm(LoopPlot)
  }
}

# Use cow package to get them into one plane: GLM vs DL
legend_b <- get_legend(plot_DL_ratio_0.5 + theme(legend.position="bottom"))
prow <- plot_grid(plot_DL_ratio_0.25, 
          plot_DL_ratio_0.5, 
          plot_DL_ratio_0.75, 
          labels = c("A", "B", "C"), ncol = 3, align = "hv",
          axis = 'tblr', hjust = -1, nrow = 1)
GLMDLci <- 
  plot_grid(prow, legend_b, ncol = 1, 
            rel_heights = c(1, 0.1));rm(legend_b, prow)

# Same for GLM vs HE
legend_b <- get_legend(plot_HE_ratio_0.5 + theme(legend.position="bottom"))
prow <- plot_grid(plot_HE_ratio_0.25, 
                  plot_HE_ratio_0.5, 
                  plot_HE_ratio_0.75, 
                  labels = c("A", "B", "C"), ncol = 3, align = "hv",
                  axis = 'tblr', hjust = -1, nrow = 1)
GLMHEci <- 
  plot_grid(prow, legend_b, ncol = 1, 
            rel_heights = c(1, 0.1));rm(legend_b, prow)

# Finally for GLM vs REML
legend_b <- get_legend(plot_REML_ratio_0.5 + theme(legend.position="bottom"))
prow <- plot_grid(plot_REML_ratio_0.25, 
                  plot_REML_ratio_0.5, 
                  plot_REML_ratio_0.75, 
                  labels = c("A", "B", "C"), ncol = 3, align = "hv",
                  axis = 'tblr', hjust = -1, nrow = 1)
GLMREMLci <- 
  plot_grid(prow, legend_b, ncol = 1, 
            rel_heights = c(1, 0.1));rm(legend_b, prow)

# If we have data for each voxel
if(AllVox){
  # Violin plots without averaging over all voxels: only for ratio 0.5 and GLM vs DL
  CoverageViolinPlot %>% 
    ungroup() %>%
    # Select the ratio
    filter(ratioBW == 0.50) %>%
    # Select the comparison
    filter(SCEN %in% comps[[1]]) %>%
    # Column for activation YES/NO and turn into factor
    mutate(signal = ifelse(TrueD == 0, 'null', 'activation')) %>%
    mutate(signalF = factor(signal, levels = c('null', 'activation'),
                          labels = c('null data', '3% BOLD signal change'))) %>%
    filter(signal == 'null') %>%
    # create labels for facets
    mutate(d = paste('d ~ "=" ~ ', TrueD, sep = ''),
           sigmaWL = paste('sigma[W] ~ "=" ~ ', round(sigmaW, 0), sep = ''),
           sigmaML = paste('sigma[M] ~ "=" ~ ', round(sigmaM, 0), sep = '')) %>%
    mutate(sigmaWLF = factor(sigmaWL, levels =
     paste('sigma[W] ~ "=" ~ ',
           round(sqrt(trueMCvalues('sim_act', 'TrueSigma2W')), 0), sep = ''))) %>%
    mutate(sigmaMLF = factor(sigmaML, levels =
     paste('sigma[M] ~ "=" ~ ',
           round(sqrt(trueMCvalues('sim_act', 'TrueSigma2M')), 0), sep = ''))) %>%
    # 4) plot the results
    ggplot(., aes(x = factor(nstud), y = coverage, fill = parameter)) + 
    # Use violin plot
    geom_violin(position = 'dodge', alpha = .75, aes(factor(nstud))) +
    facet_grid(sigmaWLF ~ sigmaMLF, labeller = label_parsed) +
    # Have mean line on top of it
    stat_summary(fun.y = mean, geom="line", aes(group = parameter), size = .7)  +
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
}

#########################################################
################### STANDARDIZED BIAS ###################
#########################################################

# Maximum bias for stand. effect sizes
BIAS %>%
  filter(SCEN != 'GLM') %>%
  mutate(AbsBias = abs(StBias)) %>%
  group_by(SCEN) %>%
  summarise(maxAbsBias = max(AbsBias),
            minBias = min(StBias))

# Preparation for the bias plot
BiasPlot <- BIAS %>%
  # Column for activation YES/NO and turn into factor
  mutate(signal = ifelse(TrueD == 0, 'null', 'activation')) %>%
  mutate(signalF = factor(signal, levels = c('null', 'activation'))) %>%
  # create labels for facets
  mutate(d = paste('d ~ "=" ~ ', TrueD, sep = ''),
         sigmaWL = paste('sigma[S] ~ "=" ~ ', round(sigmaW, 0), sep = ''),
         sigmaML = paste('sigma[M] ~ "=" ~ ', round(sigmaM, 0), sep = ''))

# Loop over the comparisons 
for(s in 1:length(comps)){
  # Loop over the ratioBWs
  for(r in 1:length(ratioBW_vec)){
    # First create variable
    LoopPlot <- 
      # Plot without the SD bars.
      BiasPlot %>%
      mutate(sigmaWLF = factor(sigmaWL, levels = unique(sigmaWL))) %>%
      mutate(sigmaMLF = factor(sigmaML, levels = unique(sigmaML))) %>%    
      # Select the ratio
      filter(ratioBW == ratioBW_vec[r]) %>% 
      # Select the comparison
      filter(SCEN %in% comps[[s]]) %>% 
      # 4) plot the results
      ggplot(., aes(x = nstud, y = StBias)) + 
      geom_line(aes(colour = parameter, linetype = signalF), size = 0.95) +
      geom_point(aes(colour = parameter, fill = parameter), size = 1.05, 
                 show.legend = FALSE) +
      scale_x_continuous(ifelse(r == 2,
                         'Number of studies in the third level', '')) +
      scale_y_continuous(ifelse(r == 1, 'Standardized bias', '')) + 
                         # I'm truncating the Y-axis otherwise you cannot see
                         # differences between the methods...
                         #limits = c(0.8, 1)) +
      #limits = c(Y_min, 1)) +
      scale_color_manual('Model',
                         labels = c(paste0('random effects MA: ', comps[[s]][2]),
                                    'mixed effects GLM'),
                         values = c(colours[colours$scen == comps[[s]][2], 'cols'],
                                    colours[colours$scen == 'GLM', 'cols'])) +
      scale_fill_manual('Model',
                        labels = c(paste0('random effects MA: ', comps[[s]][2]),
                                   'mixed effects GLM'),
                        values = c(colours[colours$scen == comps[[s]][2], 'cols'],
                                   colours[colours$scen == 'GLM', 'cols'])) +
      scale_linetype_manual('', values = c('dashed', 'solid')) +
      #ggtitle(bquote(sigma[B]/sigma[W] == .(ratioBW_vec[r]))) +
      ggtitle(paste('R = ', ratioBW_vec[r], sep = '')) +
      facet_grid(sigmaWLF ~ sigmaMLF, labeller = label_parsed) +
      theme_bw() +
      theme(legend.position = 'none',
            axis.text = element_text(face = 'bold', size = 9),
            strip.background = element_rect(fill = 'white', colour = 'white'),
            strip.text = element_text(face = 'bold', size = 12),
            title = element_text(face = 'bold', size = 12),
            legend.text = element_text(face = 'bold'),
            plot.title = element_text(hjust = 0.5))

    # Print the plot (otherwise it does not store the colours!)
    print(LoopPlot)
    
    # Wait for a second (otherwise objects are not created)
    Sys.sleep(2)
    
    # Now assign to variable
    assign(x = paste0('plot_bias_', comps[[s]][2], '_ratio_', ratioBW_vec[r]), LoopPlot)
    
    # Reset
    rm(LoopPlot)
  }
}

# Use cow package to get them into one plane: GLM vs DL
legend_b <- get_legend(plot_bias_DL_ratio_0.5 + theme(legend.position="bottom"))
prow <- plot_grid(plot_bias_DL_ratio_0.25, 
                  plot_bias_DL_ratio_0.5, 
                  plot_bias_DL_ratio_0.75, 
                  labels = c("A", "B", "C"), ncol = 3, align = "hv",
                  axis = 'tblr', hjust = -1, nrow = 1)
GLMDLbias <- 
  plot_grid(prow, legend_b, ncol = 1, 
            rel_heights = c(1, 0.1));rm(legend_b, prow)

# Same for GLM vs HE
legend_b <- get_legend(plot_bias_HE_ratio_0.5 + theme(legend.position="bottom"))
prow <- plot_grid(plot_bias_HE_ratio_0.25, 
                  plot_bias_HE_ratio_0.5, 
                  plot_bias_HE_ratio_0.75, 
                  labels = c("A", "B", "C"), ncol = 3, align = "hv",
                  axis = 'tblr', hjust = -1, nrow = 1)
GLMHEbias <- 
  plot_grid(prow, legend_b, ncol = 1, 
            rel_heights = c(1, 0.1));rm(legend_b, prow)

# Finally for GLM vs REML
legend_b <- get_legend(plot_bias_REML_ratio_0.5 + theme(legend.position="bottom"))
prow <- plot_grid(plot_bias_REML_ratio_0.25, 
                  plot_bias_REML_ratio_0.5, 
                  plot_bias_REML_ratio_0.75, 
                  labels = c("A", "B", "C"), ncol = 3, align = "hv",
                  axis = 'tblr', hjust = -1, nrow = 1)
GLMREMLbias <- 
  plot_grid(prow, legend_b, ncol = 1, 
            rel_heights = c(1, 0.1));rm(legend_b, prow)

# If we have data for each voxel
if(AllVox){
  # Violin plots without averaging over all voxels
  BIASviolin %>% 
    ungroup() %>%
    # Column for activation YES/NO and turn into factor
    mutate(signal = ifelse(TrueD == 0, 'null', 'activation')) %>%
    mutate(signalF = factor(signal, levels = c('null', 'activation'),
                          labels = c('null data', '3% BOLD signal change'))) %>%
    filter(signal == 'null') %>%
    # create labels for facets
    mutate(d = paste('d ~ "=" ~ ', TrueD, sep = ''),
           sigmaWL = paste('sigma[W] ~ "=" ~ ', round(sigmaW, 0), sep = ''),
           sigmaML = paste('sigma[M] ~ "=" ~ ', round(sigmaM, 0), sep = '')) %>%
    mutate(sigmaWLF = factor(sigmaWL, levels =
     paste('sigma[W] ~ "=" ~ ',
           round(sqrt(trueMCvalues('sim_act', 'TrueSigma2W')), 0), sep = ''))) %>%
    mutate(sigmaMLF = factor(sigmaML, levels =
     paste('sigma[M] ~ "=" ~ ',
           round(sqrt(trueMCvalues('sim_act', 'TrueSigma2M')), 0), sep = ''))) %>%
    # 4) plot the results
    ggplot(., aes(x = factor(nstud), y = StBias, fill = parameter)) + 
    # Use violin plot
    geom_violin(position = 'dodge', alpha = .75, aes(factor(nstud))) +
    facet_grid(sigmaWLF ~ sigmaMLF, labeller = label_parsed) +
    # Have mean line on top of it
    stat_summary(fun.y = mean, geom="line", aes(group = parameter), size = .7)  +
    scale_x_discrete('Number of studies in the third level') +
    scale_y_continuous('Standiardized bias') +
    scale_color_manual('Model', labels = c('random effects MA', 'mixed effects GLM'),
                       values = c('#fdb462', '#bc80bd')) +
    scale_fill_manual('Model', labels = c('random effects MA', 'mixed effects GLM'),
                      values = c('#fdb462', '#bc80bd')) +
    geom_hline(aes(yintercept = 0.95), colour = 'black', linetype = 'dashed') +
    ggtitle('ST BIAS') +
    theme_bw() +
    theme(legend.position = 'bottom',
          axis.text = element_text(face = 'bold', size = 9),
          strip.background = element_rect(fill = 'white', colour = 'white'),
          strip.text = element_text(face = 'bold', size = 12),
          title = element_text(face = 'bold', size = 12),
          legend.text = element_text(face = 'bold'),
          plot.title = element_text(hjust = 0.5))
}

# Some numbers
BIAS %>%
  filter(TrueD == 1.02) %>%
  filter(ratioBW == 0.75) %>%
  filter(sigmaW < 15) %>%
  filter(nstud == 50) %>%
  filter(sigmaM > 8)

##########################################################
################### ESTIMATED VARIANCE ###################
##########################################################

# Need to update to ratio scenario!
ratio <- FALSE
if(ratio){
  # Estimated variance components in the MA model
  EstVar %>%
    # Column for activation YES/NO and turn into factor
    mutate(signal = ifelse(TrueD == 0, 'null', 'activation')) %>%
    mutate(signalF = factor(signal, levels = c('null', 'activation'))) %>%
    filter(parameter == 'ESTTAU') %>%
    # Quick hack, need to solve this, but Tau is still squared and average of squares is not the same as square of averages....
    mutate(SqAvgSE = sqrt(AvgSE)) %>%
     ggplot(., aes(x = nstud, y = SqAvgSE, group = signalF)) +
       geom_line(aes(colour = signalF), size = 1) +
       facet_grid(sigmaW~sigmaM) +
    theme_bw()
  
  # Estimated variance components in the GLM model
  EstVar %>%
    # Column for activation YES/NO and turn into factor
    mutate(signal = ifelse(TrueD == 0, 'null', 'activation')) %>%
    mutate(signalF = factor(signal, levels = c('null', 'activation'))) %>%
    filter(parameter == 'IBMA.SE') %>%
    ggplot(., aes(x = nstud, y = AvgSE, group = signalF)) +
    geom_line(aes(colour = signalF), size = 1) +
    facet_grid(sigmaW~sigmaM) +
    theme_bw()
  
  
  # Estimated variance components and trying to add true values
  EstVar %>%
    filter(parameter == 'IBMA.SE') %>%
    ungroup() %>%
    rowwise() %>%
    mutate(TrueVAR_M = GetTrueSigmaM(sigmaW = sigmaW, nsub = 29, nstud = nstud, 
                                     sigmaM = sigmaM)) %>%
    mutate(TrueSE = sqrt(TrueVAR_M)) %>% 
    tidyr::gather(key = 'type', value = 'value', AvgSE, TrueSE, factor_key = TRUE) %>%
    dplyr::select(-parameter, -TrueVAR_M) %>%
    # No difference in activation/null data, therefore only select the activation
    filter(TrueD != 0) %>%
    # create labels for facets
    mutate(d = paste('d ~ "=" ~ ', TrueD, sep = ''),
           sigmaWL = paste('sigma[W] ~ "=" ~ ', round(sigmaW, 0), sep = ''),
           sigmaML = paste('sigma[M] ~ "=" ~ ', round(sigmaM, 0), sep = '')) %>%
    mutate(sigmaWLF = factor(sigmaWL, levels =
                               paste('sigma[W] ~ "=" ~ ',
         round(sqrt(trueMCvalues('sim_act', 'TrueSigma2W')), 0), sep = ''))) %>%
    mutate(sigmaMLF = factor(sigmaML, levels =
                               paste('sigma[M] ~ "=" ~ ',
         round(sqrt(trueMCvalues('sim_act', 'TrueSigma2M')), 0), sep = ''))) %>%
    ggplot(., aes(x = nstud, y = value)) +
    geom_line(aes(colour = type), size = 1) +
    scale_color_brewer('type', labels = c('Average estimated SE', 'True (?) SE'),
                       type = 'qual', palette = 2) +
    facet_grid(sigmaWLF~sigmaMLF, labeller = label_parsed) +
    ggtitle('Estimated versus true SE of parameter estimates',
            subtitle = 'GLM only (FLAME1)') +
    theme_bw() +
    theme(legend.position = 'bottom',
          axis.text = element_text(face = 'bold', size = 9),
          strip.background = element_rect(fill = 'white', colour = 'white'),
          strip.text = element_text(face = 'bold', size = 12),
          title = element_text(face = 'bold', size = 12),
          legend.text = element_text(face = 'bold'),
          plot.title = element_text(hjust = 0))
}


#########################################################
####################### CI LENGTH #######################
#########################################################

# Variable for plotting
CILplot <- CIL %>%
  mutate(signal = ifelse(TrueD == 0, 'null', 'activation')) %>%
  mutate(signalF = factor(signal, levels = c('null', 'activation'),
                          labels = c('null data', '3% BOLD signal change'))) %>%
  # create labels for facets
  mutate(d = paste('d ~ "=" ~ ', TrueD, sep = ''),
         sigmaWL = paste('sigma[S] ~ "=" ~ ', round(sigmaW, 0), sep = ''),
         sigmaML = paste('sigma[M] ~ "=" ~ ', round(sigmaM, 0), sep = ''),
         ratioL = paste('sigma[G]/sigma[S] ~ "=" ~ ', ratioBW))

# Averaged over voxels: only for ratio 0.5 and GLM vs DL
# For loop over the ratios
for(s in 1:length(ratioBW_vec)){
  LoopPlot <- 
    CILplot %>%
    mutate(sigmaWLF = factor(sigmaWL, levels = unique(sigmaWL))) %>%
    mutate(sigmaMLF = factor(sigmaML, levels = unique(sigmaML))) %>%
    # Select the ratio
    filter(ratioBW == ratioBW_vec[s]) %>%
    # Deselect the GLM approach
    filter(parameter != 'IBMA.COPE') %>%
    # Create the plot
    ggplot(., aes(x = nstud, y = AvgSimCIL)) + 
    geom_line(aes(colour = SCEN, linetype = signalF), size = 0.5) +
    geom_point(aes(colour = SCEN, fill = SCEN), size = 0.7, 
               show.legend = FALSE) +
    scale_x_continuous('Number of studies in the third level') +
    scale_y_continuous(ifelse(s == 1, 'Average 95% CI length', '')) +
    scale_color_manual('random effects MA: ', 
                       values = colours$cols) +
    scale_linetype_manual('', values = c('dashed', 'solid')) +
    facet_grid(sigmaWLF ~ sigmaMLF, labeller = label_parsed) +
    #ggtitle(bquote(sigma[G]/sigma[S] == .(ratioBW_vec[s]))) +
    ggtitle(paste('R = ', ratioBW_vec[r], sep = '')) +
    theme_bw() +
    theme(legend.position = 'none',
          axis.text = element_text(face = 'bold', size = 9, vjust = -1),
          strip.background = element_rect(fill = 'white', colour = 'white'),
          strip.text = element_text(face = 'bold', size = 12),
          title = element_text(face = 'bold', size = 12),
          legend.text = element_text(face = 'bold', size = 12),
          plot.title = element_text(hjust = 0.5))
  
  # Print the plot (otherwise it does not store the colours!)
  print(LoopPlot)
  
  # Wait for a second (otherwise objects are not created)
  Sys.sleep(2)
  
  # Now assign to variable
  assign(x = paste0('plot_CIL_ratio_', ratioBW_vec[s]), LoopPlot)
  
  
  # Reset
  rm(LoopPlot)
}

# Use cow package to get them into one plane: GLM vs DL
legend_b <- get_legend(plot_CIL_ratio_0.5 + theme(legend.position="bottom"))
prow <- plot_grid(plot_CIL_ratio_0.25, 
                  plot_CIL_ratio_0.5, 
                  plot_CIL_ratio_0.75, 
                  labels = c("A", "B", "C"), ncol = 3, align = "hv",
                  axis = 'tblr', hjust = -1, nrow = 1)
CILplots <- 
  plot_grid(prow, legend_b, ncol = 1, 
            rel_heights = c(1, 0.1));rm(legend_b, prow)

##
###############
### Save Figures
###############
##

# Bias
ggsave(filename = paste(saveFig, 'bias_GLMDL.png', sep = ''), 
       plot = GLMDLbias,
       width = 11.2, height = 9.18, units = 'in', scale = 1)
ggsave(filename = paste(saveFig, 'bias_GLMHE.png', sep = ''), 
       plot = GLMHEbias,
       width = 11.2, height = 9.18, units = 'in', scale = 1)
ggsave(filename = paste(saveFig, 'bias_GLMREML.png', sep = ''), 
       plot = GLMREMLbias,
       width = 11.2, height = 9.18, units = 'in', scale = 1)

# Coverage
ggsave(filename = paste(saveFig, 'ci_GLMDL.png', sep = ''), 
       plot = GLMDLci,
       width = 11.2, height = 9.18, units = 'in', scale = 1)
ggsave(filename = paste(saveFig, 'ci_GLMHE.png', sep = ''), 
       plot = GLMHEci,
       width = 11.2, height = 9.18, units = 'in', scale = 1)
ggsave(filename = paste(saveFig, 'ci_GLMREML.png', sep = ''), 
       plot = GLMREMLci,
       width = 11.2, height = 9.18, units = 'in', scale = 1)

# CI length
ggsave(filename = paste(saveFig, 'ci_lengths.png', sep = ''), 
       plot = CILplots,
       width = 11.2, height = 9.18, units = 'in', scale = 1)


##
###############
### Extra figures for discussion
###############
##

###############################################################
####################### SD COMPARISON #########################
###############################################################

# Read in the files
SDcalcDat <-
  readRDS(file = 
      paste(LIR[[currentWD]], '/SDcalcDat.rda', sep = '')) %>%
  # Remove obsolete column
  dplyr::select(-sigmaB) %>%
  # Select scenario
  filter(BOLDC %in% c(0,3)) %>%
# filter(ratioBW != 0.5) %>%
  filter(nstud %in% c(5,50)) %>%
  # Summarise over simulations
  group_by(SCEN, ratioBW, sigmaW, sigmaM, nstud, BOLDC) %>%
  summarise(AvgEst = mean(value, na.rm = TRUE),
            SDEST = sd(value, na.rm = TRUE))

# Prepare data frame
SDPlot <- SDcalcDat %>%
  # Create SD limits
  mutate(sdLow = AvgEst - SDEST,
         sdUp = AvgEst + SDEST) %>%
  # create labels for facets
  mutate(sigmaWL = paste('sigma[S] ~ "=" ~ ', round(sigmaW, 0), sep = ''),
         sigmaML = paste('sigma[M] ~ "=" ~ ', round(sigmaM, 0), sep = ''),
         ratioL = paste('sigma[G]/sigma[S] ~ "=" ~ ', ratioBW))
SDPlot$sigmaWLF <- factor(SDPlot$sigmaWL, levels = unique(SDPlot$sigmaWL))
SDPlot$sigmaMLF <- factor(SDPlot$sigmaML, levels = unique(SDPlot$sigmaML))

## --- 
# Activation condition
## --- 
SDPlotAct <- filter(SDPlot, BOLDC == 3)

# Set limits
yMin <- min(SDPlotAct$sdLow)
yMax <- max(SDPlotAct$sdUp)

# Loop over the ratioBWs
for(r in 1:length(ratioBW_vec)){
  # First create variable
  LoopPlot <- 
    SDPlotAct %>%
    mutate(nstudF = factor(nstud, levels = unique(nstud))) %>%
    # Select the ratio
    filter(ratioBW == ratioBW_vec[r]) %>% 
    ggplot(., aes(x = nstudF, y = AvgEst, group = SCEN)) +
    geom_point(aes(colour = SCEN), 
               size = 0.7,
               position = position_dodge(width = 0.5)) + 
    geom_errorbar(aes(ymin = sdLow, ymax = sdUp, colour = SCEN),
                  size = 0.8,
                  width = 0.1,
                  position = position_dodge(width = 0.5)) +
    scale_color_manual('random effects MA (data with activation): ', 
                       values = colours$cols) +
    facet_grid(sigmaWLF ~ sigmaMLF, labeller = label_parsed,
               scales = 'free_y') +
    scale_x_discrete(ifelse(r == 2,
                'Number of studies in the third level', '')) +
    scale_y_continuous(ifelse(r == 1, 'estimated population effect size',
                              ''),
                       limits = c(yMin, yMax)) +
    #ggtitle(bquote(sigma[G]/sigma[S] == .(ratioBW_vec[r]))) +
    ggtitle(paste('R = ', ratioBW_vec[r], sep = '')) +
    theme_bw() +
    theme(legend.position = 'none',
          axis.text = element_text(face = 'bold', size = 9),
          strip.background = element_rect(fill = 'white', colour = 'white'),
          strip.text = element_text(face = 'bold', size = 12),
          title = element_text(face = 'bold', size = 12),
          legend.text = element_text(face = 'bold'),
          plot.title = element_text(hjust = 0.5))

  # Print the plot (otherwise it does not store the colours!)
  print(LoopPlot)
  
  # Wait for a second (otherwise objects are not created)
  Sys.sleep(2)
  
  # Now assign to variable
  assign(x = paste0('plot_CIlength_ratio_', ratioBW_vec[r]), LoopPlot)
  
  # Reset
  rm(LoopPlot)
}
  
# Combine plots using cow package
legend_b <- get_legend(plot_CIlength_ratio_0.5 + theme(legend.position="bottom"))
prow <- plot_grid(plot_CIlength_ratio_0.25, 
                  plot_CIlength_ratio_0.5, 
                  plot_CIlength_ratio_0.75, 
                  labels = c("A", "B", "C"), ncol = 3, align = "hv",
                  axis = 'tblr', hjust = -1, nrow = 1)
SDlengths <- plot_grid(prow, legend_b, ncol = 1, 
            rel_heights = c(1, 0.1));rm(legend_b, prow)

# CI length
ggsave(filename = paste(saveFig, 'SD_estimates.png', sep = ''), 
       plot = SDlengths,
       width = 10, height = 9, units = 'in', scale = 1)


## --- 
# Null condition
## --- 
SDPlotNull <- filter(SDPlot, BOLDC == 0)

# Set limits
yMin <- min(SDPlotNull$sdLow)
yMax <- max(SDPlotNull$sdUp)

# Loop over the ratioBWs
for(r in 1:length(ratioBW_vec)){
  # First create variable
  LoopPlot <- 
    SDPlotNull %>%
    mutate(nstudF = factor(nstud, levels = unique(nstud))) %>%
    # Select the ratio
    filter(ratioBW == ratioBW_vec[r]) %>% 
    ggplot(., aes(x = nstudF, y = AvgEst, group = SCEN)) +
    geom_point(aes(colour = SCEN), 
               size = 0.7,
               position = position_dodge(width = 0.5)) + 
    geom_errorbar(aes(ymin = sdLow, ymax = sdUp, colour = SCEN),
                  size = 0.8,
                  width = 0.1,
                  position = position_dodge(width = 0.5)) +
    scale_color_manual('random effects MA (null data): ', 
                       values = colours$cols) +
    facet_grid(sigmaWLF ~ sigmaMLF, labeller = label_parsed,
               scales = 'free_y') +
    scale_x_discrete(ifelse(r == 2,
                            'Number of studies in the third level', '')) +
    scale_y_continuous(ifelse(r == 1, 'estimated population effect size',
                              ''),
                       limits = c(yMin, yMax)) +
    #ggtitle(bquote(sigma[G]/sigma[S] == .(ratioBW_vec[r]))) +
    ggtitle(paste('R = ', ratioBW_vec[r], sep = '')) +
    theme_bw() +
    theme(legend.position = 'none',
          axis.text = element_text(face = 'bold', size = 9),
          strip.background = element_rect(fill = 'white', colour = 'white'),
          strip.text = element_text(face = 'bold', size = 12),
          title = element_text(face = 'bold', size = 12),
          legend.text = element_text(face = 'bold'),
          plot.title = element_text(hjust = 0.5))
  
  # Print the plot (otherwise it does not store the colours!)
  print(LoopPlot)
  
  # Wait for a second (otherwise objects are not created)
  Sys.sleep(2)
  
  # Now assign to variable
  assign(x = paste0('plot_CIlength_ratio_', ratioBW_vec[r]), LoopPlot)
  
  # Reset
  rm(LoopPlot)
}

# Combine plots using cow package
legend_b <- get_legend(plot_CIlength_ratio_0.5 + theme(legend.position="bottom"))
prow <- plot_grid(plot_CIlength_ratio_0.25, 
                  plot_CIlength_ratio_0.5, 
                  plot_CIlength_ratio_0.75, 
                  labels = c("A", "B", "C"), ncol = 3, align = "hv",
                  axis = 'tblr', hjust = -1, nrow = 1)
SDlengths_null <- plot_grid(prow, legend_b, ncol = 1, 
               rel_heights = c(1, 0.1));rm(legend_b, prow)

# CI length
ggsave(filename = paste(saveFig, 'SD_estimates_null.png', sep = ''), 
       plot = SDlengths_null,
       width = 10, height = 9, units = 'in', scale = 1)

###############################################################
####################### CI ILLUSTRATION #######################
###############################################################

# Read in the data
illuCIL <- 
  readRDS(file = 
    paste(LIR[[currentWD]], '/illuCIL.rda', sep = '')) 

# Comparison between two ratio's
SelRatios <- c(0.25, 0.75)

# SD over simulations for DL vs HE scenarios and ratio = 0.25 vs 0.75
IlluValue <- illuCIL %>%
  group_by(SCEN, ratioBW, sigmaW, nstud, TrueG) %>%
  summarise(SDest = sd(value),
            AvgEst = mean(value)) %>%
  # Add simulation at position - 2
  mutate(sim = -2)
    # --> HE is more variable!

# Not possible to set the limits: results not visible for ratio = 0.25!
#xMin <- min(illuCIL$CI_lower)
#xMax <- max(illuCIL$CI_upper)

# For loop over both ratios
for(r in 1:length(SelRatios)){
  # Select the ratio
  SelRatio <- SelRatios[r]
  
  # Create the plot
  LoopPlot <- 
  illuCIL %>%
    filter(ratioBW == SelRatio) %>%
    ggplot(aes(x = value, y = sim)) +
    geom_point(size = 0.7) + 
    geom_segment(aes(x = CI_lower, xend = CI_upper, 
                     y = factor(sim), yend = factor(sim), 
                     colour = factor(cov_IND)),
                 size = 0.6) +
    geom_vline(aes(xintercept = TrueG)) +
    scale_x_continuous('') +
    scale_y_discrete(ifelse(r == 1, 'simulation', ''),
      breaks = c('1', '25', '50', '75', '100')) +
    facet_grid(.~SCEN) +
    scale_color_manual('Contains true parameter:', labels = c('NO',  'YES'),
                       values = c('#d95f02', '#1b9e77')) +
    geom_point(data = IlluValue %>%
                 filter(ratioBW == SelRatio), 
               aes(x = AvgEst, y = sim), 
               shape = 23, size = 3, fill = 'black') +
    #ggtitle(bquote(sigma[G]/sigma[S] == .(SelRatio))) +
    ggtitle(paste('R = ', SelRatio, sep = '')) +
    theme_bw() +
    theme(legend.position = 'none',
          axis.text = element_text(face = 'bold', size = 9),
          strip.background = element_rect(fill = 'white', colour = 'white'),
          strip.text = element_text(face = 'bold', size = 12),
          title = element_text(face = 'bold', size = 12),
          legend.text = element_text(face = 'bold', size = 12),
          plot.title = element_text(hjust = 0.5))
    
  # Print the plot (otherwise it does not store the colours!)
  print(LoopPlot)
  
  # Wait for a second (otherwise objects are not created)
  Sys.sleep(2)
  
  # Now assign to variable
  assign(x = paste0('plot_illuCIL_ratio_', SelRatio), LoopPlot)
}

# Combine plots using cow package
legend_b <- get_legend(plot_illuCIL_ratio_0.25 + theme(legend.position="bottom"))
prow <- plot_grid(plot_illuCIL_ratio_0.75,
                  plot_illuCIL_ratio_0.25,
                  labels = c("A", "B"), ncol = 2, align = "hv",
                  axis = 'tblr', hjust = -1, nrow = 1)
CI_illus <- plot_grid(prow, legend_b, ncol = 1, 
                  rel_heights = c(1, 0.1));rm(legend_b, prow)


# CI illustration
ggsave(filename = paste(saveFig, 'CI_illusts.png', sep = ''), 
       plot = CI_illus,
       width = 7, height = 7, units = 'in', scale = 1)



##
###############
### Figures for OHBM 2019
###############
##

OHBMCov <- CoveragePlot %>%
  # Column for activation YES/NO and turn into factor
  mutate(signal = ifelse(TrueD == 0, 'null', 'activation')) %>%
  mutate(signalF = factor(signal, levels = c('null', 'activation'),
                          labels = c('null data', '3% BOLD signal change'))) %>%
  # create labels for facets
  mutate(d = paste('d ~ "=" ~ ', TrueD, sep = ''),
         sigmaWL = paste('sigma[S] ~ "=" ~ ', round(sigmaW, 0), sep = ''),
         sigmaML = paste('sigma[M] ~ "=" ~ ', round(sigmaM, 0), sep = ''),
         ratioL = paste('sigma[G]/sigma[S] ~ "=" ~ ', ratioBW)) %>%
  # Only have ratio = 0.5 and HE estimator
  filter(SCEN %in% c('GLM', 'HE')) %>%
  filter(ratioBW == 0.5)

OHBMcovPlot <- OHBMCov %>%
  mutate(sigmaWLF = factor(sigmaWL, levels = unique(sigmaWL))) %>%
  mutate(sigmaMLF = factor(sigmaML, levels = unique(sigmaML))) %>%
  # 4) plot the results
  ggplot(., aes(x = nstud, y = coverage)) +
  geom_line(aes(colour = parameter, linetype = signalF), size = 0.9) +
  geom_point(aes(colour = parameter, fill = parameter), size = 0.95, show.legend = FALSE) +
  scale_x_continuous('Number of studies in the third level') +
  scale_y_continuous('Empirical coverage')+
  scale_color_manual('Model',
                     labels = c('random effects MA', 'mixed effects GLM'),
                     values = c('#8da0cb', '#66c2a5')) +
  scale_fill_manual('Model',
                     labels = c('random effects MA', 'mixed effects GLM'),
                     values = c('#8da0cb', '#66c2a5')) +
  scale_linetype_manual('', values = c('dashed', 'solid')) +
  geom_hline(aes(yintercept = 0.95), colour = 'black') +
  ggtitle(bquote(sigma[G]^2/sigma[S]^2 == 0.5)) +
  facet_grid(sigmaWLF ~ sigmaMLF, labeller = label_parsed) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_text(face = 'bold', size = 9, vjust = -1),
        strip.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(face = 'bold', size = 12),
        title = element_text(face = 'bold', size = 12),
        legend.text = element_text(face = 'bold', size = 12),
        plot.title = element_text(hjust = 0.5))
OHBMcovPlot


# Preparation for the bias plot (here we use average bias instead of standardized!)
OHBMbias <- BIAS %>%
  # Column for activation YES/NO and turn into factor
  mutate(signal = ifelse(TrueD == 0, 'null', 'activation')) %>%
  mutate(signalF = factor(signal, levels = c('null', 'activation'))) %>%
  # create labels for facets
  mutate(d = paste('d ~ "=" ~ ', TrueD, sep = ''),
         sigmaWL = paste('sigma[S] ~ "=" ~ ', round(sigmaW, 0), sep = ''),
         sigmaML = paste('sigma[M] ~ "=" ~ ', round(sigmaM, 0), sep = ''),
         ratioL = paste('sigma[G]/sigma[S] ~ "=" ~ ', ratioBW)) %>%
  # Only have ratio = 0.5 and HE estimator
  filter(SCEN %in% c('GLM', 'HE')) %>%
  filter(ratioBW == 0.5)

OHBMbiasPlot <- OHBMbias %>%
  mutate(sigmaWLF = factor(sigmaWL, levels = unique(sigmaWL))) %>%
  mutate(sigmaMLF = factor(sigmaML, levels = unique(sigmaML))) %>%    
  # 4) plot the results
  ggplot(., aes(x = nstud, y = AvgBias)) + 
  geom_line(aes(colour = parameter, linetype = signalF), size = 0.95) +
  geom_point(aes(colour = parameter, fill = parameter), size = 1.05, 
             show.legend = FALSE) +
  scale_x_continuous('Number of studies in the third level') +
  scale_y_continuous('Standardized bias') + 
  scale_color_manual('Model',
                     labels = c('random effects MA', 'mixed effects GLM'),
                     values = c('#8da0cb', '#66c2a5')) +
  scale_fill_manual('Model',
                    labels = c('random effects MA', 'mixed effects GLM'),
                    values = c('#8da0cb', '#66c2a5')) +
  scale_linetype_manual('', values = c('dashed', 'solid')) +
  ggtitle(bquote(sigma[G]^2/sigma[S]^2 == 0.5)) +
  facet_grid(sigmaWLF ~ sigmaMLF, labeller = label_parsed) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_text(face = 'bold', size = 9, vjust = -1),
        strip.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(face = 'bold', size = 12),
        title = element_text(face = 'bold', size = 12),
        legend.text = element_text(face = 'bold', size = 12),
        plot.title = element_text(hjust = 0.5))
OHBMbiasPlot

# Create figure with both plots
legend_b <- get_legend(OHBMcovPlot + theme(legend.position="bottom"))
prow <- plot_grid(OHBMbiasPlot,
                  OHBMcovPlot,
                  labels = c("A", "B"), ncol = 2, align = "hv",
                  axis = 'tblr', hjust = -1, nrow = 1)
OHBMplots <- plot_grid(prow, legend_b, ncol = 1, 
              rel_heights = c(1, 0.1))

# Save figures
ggsave(filename = 
         '/Users/hanbossier/Dropbox/PhD/Events/OHBM2019/Poster/Figures/cov_OHBM.png', 
       plot = OHBMcovPlot,
       width = 7, height = 9, units = 'in', scale = 1)

ggsave(filename = 
  '/Users/hanbossier/Dropbox/PhD/Events/OHBM2019/Poster/Figures/bias_OHBM.png', 
       plot = OHBMbiasPlot,
       width = 7, height = 9, units = 'in', scale = 1)


# Both in one figure (shared legend)
ggsave(filename = 
         '/Users/hanbossier/Dropbox/PhD/Events/OHBM2019/Poster/Figures/OHBM_bias_cov.png', 
       plot = OHBMplots,
       width = 11, height = 9, units = 'in', scale = 1)








