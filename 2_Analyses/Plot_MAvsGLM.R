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

# Plotting the intermediate analyses of MA vs GLM approach. 
# Contains null and activation data.
# 3-level simulation approach.

# Blocked design for individual subjects.
# One condition.
# These N subjects are pooled using FLAME pooling
# The resulting images are converted to either Hedges' g and pooled using random effects meta-analysis.
# OR using 3e level GLM with again FLAME1.

# Data contains between-study heterogeneity

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

# Intermediate Results: this is location where summarized results are written
LIR <- list(
  'Take[GLMvsMA_wi_w_act]' = 
    "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/GLMvsMA_wi_w_act/ProcessedResults",
  'Take[Estimators]' = "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/Estimators/D_TrueVal/"
)
currentWD <- 2

# Number of conficence intervals
CIs <- c('MA-weightVar','GLM-t')
NumCI <- length(CIs)

# Data frame with number of simulations and subjects for current simulation
info <- data.frame('Sim' = c(1,2),
                   'nsim' = c(1000, 1000),
                   'nsub' = rep(trueMCvalues('sim_act', 'nsub'),2))
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
    summarise(coverage = mean(wSimCoverage),
              sdCoverage = sd(wSimCoverage)) %>%
    bind_rows(CoveragePlot, .)

  ### COVERAGE AVERAGED OVER ALL SIMULATIONS ###
  if(AllVox){
    CoverageViolinPlot <- readRDS(file = 
      paste(LIR[[currentWD]], '/', SCEN[s], '/coverage_all.rda', sep = '')) %>%
      group_by(voxel, parameter, SCEN, TrueD, ratioBW, sigmaW, sigmaM, nstud) %>%
      summarise(coverage = mean(cov_IND)) %>%
      bind_rows(CoverageViolinPlot, .)
  }

  ### CI LENGHTS AVERAGED OVER ALL VOXELS ###
  CILdata <- readRDS(file = 
    paste(LIR[[currentWD]], '/', SCEN[s], '/CIlength_Vox.rda', sep = ''))
  # Summarise over simulations
  CIL <- CILdata %>% group_by(parameter, SCEN, TrueD, ratioBW, sigmaW, sigmaM, nstud) %>%
    summarise(AvgSimCIL = mean(AvCIlength)) %>%
    bind_rows(CIL, .)
  
  ### STANDARDIZED BIAS AVERAGED OVER ALL VOXELS ###
  BiasData <- readRDS(file = 
    paste(LIR[[currentWD]], '/', SCEN[s], '/bias_Vox.rda', sep = ''))
  # Summarise over simulations
  BIAS <- BiasData %>% group_by(parameter, SCEN, TrueD, ratioBW, sigmaW, sigmaM, nstud) %>%
    summarise(AvgBias = mean(Avbias),
              SDBias = sd(Avbias)) %>%
    mutate(StBias = AvgBias/SDBias * 100) %>%
    bind_rows(BIAS, .)
  
  ### STANDARDIZED BIAS AVERAGED OVER ALL SIMULATIONS ###
  if(AllVox){
    BiasData_all <- readRDS(file = 
        paste(LIR[[currentWD]], '/', SCEN[s], '/bias_all.rda', sep = ''))
    # Summarise over simulations
    BIASviolin <- BiasData_all %>% 
      group_by(voxel, parameter, SCEN, TrueD, ratioBW, sigmaW, sigmaM, nstud) %>%
      summarise(AvgBias = mean(bias),
                SDBias = sd(bias)) %>%
      mutate(StBias = AvgBias/SDBias * 100) %>%
      bind_rows(BIASviolin, .)
  }

  ### ESTIMATED BETWEEN-STUDY VARIANCE AVERAGED OVER ALL VOXELS ###
  EstVarData <- readRDS(file = 
    paste(LIR[[currentWD]], '/', SCEN[s], '/EstVar_Vox.rda', sep = ''))
  # Summarise over simulations
  EstVar <- EstVarData %>% 
    group_by(parameter, SCEN, TrueD, ratioBW, sigmaW, sigmaM, nstud) %>%
    summarise(AvgSE = mean(AvgEstSE)) %>%
    bind_rows(EstVar, .)
}


#########################################################
###################### CI COVERAGE ######################
#########################################################


#### FIRST SECTION: AVERAGED OVER ALL VOXELS AND THEN SIMULATIONS ####
CoveragePlot <- CoveragePlot %>%
  # Column for activation YES/NO and turn into factor
  mutate(signal = ifelse(TrueD == 0, 'null', 'activation')) %>%
  mutate(signalF = factor(signal, levels = c('null', 'activation'),
                          labels = c('null data', '3% BOLD signal change'))) %>%
  # create labels for facets
  mutate(d = paste('d ~ "=" ~ ', TrueD, sep = ''),
         sigmaWL = paste('sigma[W] ~ "=" ~ ', round(sigmaW, 0), sep = ''),
         sigmaML = paste('sigma[M] ~ "=" ~ ', round(sigmaM, 0), sep = ''),
         ratioL = paste('sigma[B]/sigma[W] ~ "=" ~ ', ratioBW))

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
      ggtitle(bquote(sigma[B]/sigma[W] == .(ratioBW_vec[r]))) +
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
    
    # Print the plot
    print(LoopPlot)
    
    # Wait for a second
    Sys.sleep(1)
    
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
plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, 0.1))

# Same for GLM vs HE
legend_b <- get_legend(plot_HE_ratio_0.5 + theme(legend.position="bottom"))
prow <- plot_grid(plot_HE_ratio_0.25, 
                  plot_HE_ratio_0.5, 
                  plot_HE_ratio_0.75, 
                  labels = c("A", "B", "C"), ncol = 3, align = "hv",
                  axis = 'tblr', hjust = -1, nrow = 1)
plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, 0.1))

# Finally for GLM vs REML
legend_b <- get_legend(plot_REML_ratio_0.5 + theme(legend.position="bottom"))
prow <- plot_grid(plot_REML_ratio_0.25, 
                  plot_REML_ratio_0.5, 
                  plot_REML_ratio_0.75, 
                  labels = c("A", "B", "C"), ncol = 3, align = "hv",
                  axis = 'tblr', hjust = -1, nrow = 1)
plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, 0.1))

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
####################### CI LENGTH #######################
#########################################################


# Averaged over voxels: only for ratio 0.5 and GLM vs DL
CIL %>%
  # Select the ratio
  filter(ratioBW == 0.50) %>%
  # Select the comparison
  filter(SCEN %in% comps[[1]]) %>%
  # Column for activation YES/NO and turn into factor
  mutate(signal = ifelse(TrueD == 0, 'null', 'activation')) %>%
  mutate(signalF = factor(signal, levels = c('null', 'activation'))) %>%
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
  ggplot(., aes(x = nstud, y = AvgSimCIL)) + 
  geom_line(aes(colour = parameter, linetype = signalF), size = 0.5) +
  geom_point(aes(colour = parameter, fill = parameter), size = 1.05, 
             show.legend = FALSE) +
  scale_x_continuous('Number of studies in the third level') +
  scale_y_continuous('Average 95% CI length') +
  scale_color_manual('Model', labels = c('Random effects MA', 'FLAME'),
                     values = c('#fdb462', '#bc80bd')) +
  scale_fill_manual('Model', labels = c('Random effects MA', 'FLAME'),
                    values = c('#fdb462', '#bc80bd')) +
  scale_linetype_manual('', values = c('dashed', 'solid')) +
  facet_grid(sigmaWLF ~ sigmaMLF, labeller = label_parsed) +
  theme_bw() +
  theme(legend.position = 'bottom',
        axis.text = element_text(face = 'bold', size = 9),
        strip.background = element_rect(fill = 'white', colour = 'white'),
        strip.text = element_text(face = 'bold', size = 12),
        title = element_text(face = 'bold', size = 12),
        legend.text = element_text(face = 'bold'),
        plot.title = element_text(hjust = 0))



#########################################################
################### STANDARDIZED BIAS ###################
#########################################################

# Preparation for the bias plot
BiasPlot <- BIAS %>%
  # Column for activation YES/NO and turn into factor
  mutate(signal = ifelse(TrueD == 0, 'null', 'activation')) %>%
  mutate(signalF = factor(signal, levels = c('null', 'activation'))) %>%
  # create labels for facets
  mutate(d = paste('d ~ "=" ~ ', TrueD, sep = ''),
         sigmaWL = paste('sigma[W] ~ "=" ~ ', round(sigmaW, 0), sep = ''),
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
      ggplot(., aes(x = nstud, y = AvgBias)) + 
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
      #scale_x_continuous('Number of studies in the third level') +
      #scale_y_continuous('Standardized bias') +
      #scale_color_manual('Model', labels = c('random effects MA', 'mixed effects GLM'),
      #                   values = c('#fdb462', '#bc80bd')) +
      #scale_fill_manual('Model', labels = c('random effects MA', 'mixed effects GLM'),
      #                  values = c('#fdb462', '#bc80bd')) +
      scale_linetype_manual('', values = c('dashed', 'solid')) +
#      ggtitle('st bias = mean(estimate - true value) / sd(estimate - true value)',
#              subtitle = 'calculated over simulations') +
      ggtitle(bquote(sigma[B]/sigma[W] == .(ratioBW_vec[r]))) +
      facet_grid(sigmaWLF ~ sigmaMLF, labeller = label_parsed) +
      theme_bw() +
      theme(legend.position = 'none',
            axis.text = element_text(face = 'bold', size = 9),
            strip.background = element_rect(fill = 'white', colour = 'white'),
            strip.text = element_text(face = 'bold', size = 12),
            title = element_text(face = 'bold', size = 12),
            legend.text = element_text(face = 'bold'),
            plot.title = element_text(hjust = 0))
    
      # CoveragePlot %>%
      # mutate(sigmaWLF = factor(sigmaWL, levels = unique(sigmaWL))) %>%
      # mutate(sigmaMLF = factor(sigmaML, levels = unique(sigmaML))) %>%
      # # Select the ratio
      # filter(ratioBW == ratioBW_vec[r]) %>% 
      # # Select the comparison
      # filter(SCEN %in% comps[[s]]) %>% 
      # # 4) plot the results
      # ggplot(., aes(x = nstud, y = coverage)) +
      # geom_line(aes(colour = parameter, linetype = signalF), size = 0.9) +
      # geom_point(aes(colour = parameter, fill = parameter), size = 0.95, show.legend = FALSE) +
      # scale_x_continuous(ifelse(r == 2,
      #                           'Number of studies in the third level', '')) +
      # scale_y_continuous(ifelse(r == 1, 'Empirical coverage', ''),
      #                    # I'm truncating the Y-axis otherwise you cannot see
      #                    # differences between the methods...
      #                    limits = c(0.8, 1)) +
      # #limits = c(Y_min, 1)) +
      # scale_color_manual('Model',
      #                    labels = c(paste0('random effects MA: ', comps[[s]][2]),
      #                               'mixed effects GLM'),
      #                    values = c(colours[colours$scen == comps[[s]][2], 'cols'],
      #                               colours[colours$scen == 'GLM', 'cols'])) +
      # scale_fill_manual('Model',
      #                   labels = c(paste0('random effects MA: ', comps[[s]][2]),
      #                              'mixed effects GLM'),
      #                   values = c(colours[colours$scen == comps[[s]][2], 'cols'],
      #                              colours[colours$scen == 'GLM', 'cols'])) +
      # scale_linetype_manual('', values = c('dashed', 'solid')) +
      # geom_hline(aes(yintercept = 0.95), colour = 'black') +
      # ggtitle(bquote(sigma[B]/sigma[W] == .(ratioBW_vec[r]))) +
      # # ggtitle(paste0('Ratio: ', ratioBW_vec[r])) +
      # facet_grid(sigmaWLF ~ sigmaMLF, labeller = label_parsed) +
      # theme_bw() +
      # theme(legend.position = 'none',
      #       axis.text = element_text(face = 'bold', size = 9, vjust = -1),
      #       strip.background = element_rect(fill = 'white', colour = 'white'),
      #       strip.text = element_text(face = 'bold', size = 12),
      #       title = element_text(face = 'bold', size = 12),
      #       legend.text = element_text(face = 'bold', size = 12),
      #       plot.title = element_text(hjust = 0.5))
    
    # Print the plot
    print(LoopPlot)
    
    # Wait for a second
    Sys.sleep(1)
    
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
plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, 0.1))

# Same for GLM vs HE
legend_b <- get_legend(plot_bias_HE_ratio_0.5 + theme(legend.position="bottom"))
prow <- plot_grid(plot_bias_HE_ratio_0.25, 
                  plot_bias_HE_ratio_0.5, 
                  plot_bias_HE_ratio_0.75, 
                  labels = c("A", "B", "C"), ncol = 3, align = "hv",
                  axis = 'tblr', hjust = -1, nrow = 1)
plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, 0.1))

# Finally for GLM vs REML
legend_b <- get_legend(plot_bias_REML_ratio_0.5 + theme(legend.position="bottom"))
prow <- plot_grid(plot_bias_REML_ratio_0.25, 
                  plot_bias_REML_ratio_0.5, 
                  plot_bias_REML_ratio_0.75, 
                  labels = c("A", "B", "C"), ncol = 3, align = "hv",
                  axis = 'tblr', hjust = -1, nrow = 1)
plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, 0.1))


# # Plot without the SD bars.
#   BIAS %>%
#   # Column for activation YES/NO and turn into factor
#   mutate(signal = ifelse(TrueD == 0, 'null', 'activation')) %>%
#   mutate(signalF = factor(signal, levels = c('null', 'activation'))) %>%
#   # create labels for facets
#   mutate(d = paste('d ~ "=" ~ ', TrueD, sep = ''),
#          sigmaWL = paste('sigma[W] ~ "=" ~ ', round(sigmaW, 0), sep = ''),
#          sigmaML = paste('sigma[M] ~ "=" ~ ', round(sigmaM, 0), sep = '')) %>%
#   mutate(sigmaWLF = factor(sigmaWL, levels =
#      paste('sigma[W] ~ "=" ~ ',
#        round(sqrt(trueMCvalues('sim_act', 'TrueSigma2W')), 0), sep = ''))) %>%
#   mutate(sigmaMLF = factor(sigmaML, levels =
#      paste('sigma[M] ~ "=" ~ ',
#        round(sqrt(trueMCvalues('sim_act', 'TrueSigma2M')), 0), sep = ''))) %>%
#   # 4) plot the results
#   ggplot(., aes(x = nstud, y = AvgBias)) + 
#   geom_line(aes(colour = parameter, linetype = signalF), size = 0.95) +
#   geom_point(aes(colour = parameter, fill = parameter), size = 1.05, 
#              show.legend = FALSE) +
#   scale_x_continuous('Number of studies in the third level') +
#   scale_y_continuous('Standardized bias') +
#   scale_color_manual('Model', labels = c('random effects MA', 'mixed effects GLM'),
#                      values = c('#fdb462', '#bc80bd')) +
#   scale_fill_manual('Model', labels = c('random effects MA', 'mixed effects GLM'),
#                     values = c('#fdb462', '#bc80bd')) +
#   scale_linetype_manual('', values = c('dashed', 'solid')) +
#   ggtitle('st bias = mean(estimate - true value) / sd(estimate - true value)',
#           subtitle = 'calculated over simulations') +
#   facet_grid(sigmaWLF ~ sigmaMLF, labeller = label_parsed) +
#   theme_bw() +
#   theme(legend.position = 'bottom',
#         axis.text = element_text(face = 'bold', size = 9),
#         strip.background = element_rect(fill = 'white', colour = 'white'),
#         strip.text = element_text(face = 'bold', size = 12),
#         title = element_text(face = 'bold', size = 12),
#         legend.text = element_text(face = 'bold'),
#         plot.title = element_text(hjust = 0))


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

##########################################################
################### ESTIMATED VARIANCE ###################
##########################################################


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






