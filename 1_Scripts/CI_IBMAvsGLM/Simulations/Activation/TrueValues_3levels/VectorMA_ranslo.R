####################
#### TITLE:   Vector simulation: GLM with random slope
#### Contents:
####
#### Source Files:
#### First Modified: 06/02/2018
#### Notes:
#################

##
###############
### Notes
###############
##

# Simulate fMRI time series per subject using a GLM with a random slope for the study.
# I will simulate on grid, but only save middle voxel (no smoothing).

##
###############
### Preparation
###############
##
t1 <- Sys.time()

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
  wd <- '/user/scratch/gent/gvo000/gvo00022/vsc40728/VectorMA_RanSlope'
}
if(MACHINE=='MAC'){
  wd <- '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/VectorMA_RanSlope'
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
library(Hmisc)
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


###################
#### Global variables
###################

# Number of subject
nsub <- 50

# Number of studies
nstud <- 80

# Value for sigma in the model
sigma_m <- 100

# Value for eta in the model
eta_m <- 25

###################
#### Data characteristics
###################

# Signal characteristics
TR <- 2
nscan <- 200
total <- TR*nscan
on1 <- seq(1,total,40)
onsets <- list(on1)
duration <- list(20)

###################
#### Generate a design: GROUND TRUTH DESIGN
###################

# true %BOLD change
BOLDC <- 3

# Base/intercept of signal
base <- 100

#######################################
#### DESIGN AND SIGNAL TIME SERIES ####
#######################################

# Generating a design matrix: convolution of block design with double-gamma HRF
X <- neuRosim::simprepTemporal(total,1,onsets = onsets,
                               effectsize = 1, durations = duration,
                               TR = TR, acc = 0.1, hrf = "double-gamma")

# X vector for one subject = predicted signal
X_s <- neuRosim::simTSfmri(design=X, base=0, SNR=1, noise="none", verbose=FALSE)

# Now the model will be: base + (BOLDC + b) * pred

# Now we create the BOLD/beta 1 signal 
# signal_BOLDC <- BOLDC * (pred-base) + base

## Design parameters
# Extend the design matrix with the intercept
xIN <- cbind(base,X_s)

# Contrast: not interested in intercept
CONTRAST <- matrix(c(0,1),nrow=1)

# Calculate (X'X)^(-1) with contrast
design_factor <- CONTRAST %*% (solve(t(xIN) %*% xIN )) %*% t(CONTRAST)


##################
#### GENERATE DATA
##################

# Empty lmer results data frame
LMER_res <- MA_res <- comb_res <- data.frame() %>% as_tibble()

# Generate D matrix: variance-covariance matrix of random intercept + slope
# Variance of slope = eta_m**2
var_cov_D <- rbind(c(1.0**2, 0), c(0, eta_m**2))
# Generate values using this D-matrix for intercept and slope per study
B_matrix <- MASS::mvrnorm(nstud, mu=c(0,0), Sigma=var_cov_D)

# Start 10 iterations (increases efficiency since iterations run very fast)
for(ID in startIndex:endIndex){
  # Set starting seed
  starting.seed <- pi*ID
  set.seed(starting.seed)
  
  # Empty vector
  Y <- data.frame() %>% as_tibble()
  
  # For loop over studies
  for(t in 1:nstud){
    # Take values for b_slope and b_int
    b_int <- B_matrix[t,1]
    b_sl <- B_matrix[t,2]
    
    # For loop over all subjects
    for(i in 1:nsub){
    # Generate nscan values, corresponding to time series of one subject 
      # within a study
    Y_s <- base + b_int + ((BOLDC + b_sl) * X_s) + rnorm(n = nscan, mean = 0, sd = sigma_m)
    
    # Add to data frame
    Y <- data.frame(Y = Y_s, X = X_s, sub = i, stud = t) %>% as_tibble() %>%
      bind_rows(Y, .)
    }
  }
  #####################################
  #### LINEAR MIXED MODEL APPROACH ####
  # Fit model with random intercept per subject and study
    # and random slope for study. Get coefficients using tidy.
  LMER_res <- broom::tidy(lmer(Y ~ 1 + X + (1 | sub) + (1 + X|stud), data = Y)) %>%
    as_tibble() %>% mutate(sim = ID) 
  # %>% bind_rows(LMER_res, .)
  
  ##################################################
  # Meta-analysis using standardized effect sizes ##
  # First calculate per study the beta estimate using OLS
  StudyData <- 
    Y %>% group_by(stud) %>%
    nest() %>%
    mutate(T_study = purrr::map_dbl(data, ~ 
      summary(lm(Y ~ 1 + X + sub, data = .x))$coefficients['X', 't value'])) %>%
    dplyr::select(-data) %>%
    mutate(nsub = nsub) %>%
    mutate(HedgeG = hedgeG(t = T_study, N = nsub)) %>%
    mutate(varHedge = varHedge(g = HedgeG, N = nsub))
  
    # Now do meta-analysis
    MA_res <- data.frame(
      estimate = metafor::rma(yi = HedgeG, vi = varHedge, method = 'DL', 
                     data = StudyData)$beta['intrcpt',1],
      tau = sqrt(metafor::rma(yi = HedgeG, vi = varHedge, method = 'DL', 
                data = StudyData)$tau2)) %>%
      mutate(sim = ID) %>% as_tibble() 
    # %>% bind_rows(MA_res,.)

    # Combine LMER and MA
    comb_res <- data.frame(estimate = c(filter(LMER_res, term == 'X') %>%
                              dplyr::select(estimate) %>% as.numeric(.),
                            MA_res$estimate),
               parameter = c('BETA', 'STAN_ES'),
               EstTau = c(filter(LMER_res, term == 'sd_X.stud') %>%
                            dplyr::select(estimate) %>% as.numeric(.),
                          MA_res$tau),
               model = c('LMER', 'MA'),
               sim = ID) %>%
      bind_rows(comb_res,.)
}

Sys.time() - t1 

# Save R object
saveRDS(comb_res, file = paste0(wd, '/Results/LMER_', hpcID, '.rda'))
     

  
  
  
  
  
  
  
  
  
  











    
    
    
    
    
    
    
    
    
    

























