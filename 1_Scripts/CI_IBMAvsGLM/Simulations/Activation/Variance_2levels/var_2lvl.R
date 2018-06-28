####################
#### TITLE:   Variance estimation at group level + CI coverage
#### Contents:
####
#### Source Files:
#### First Modified: 14/06/2018
#### Notes:
#################

##
###############
### Notes
###############
##

# Simulate fMRI time series per subject using a GLM with between-subject variability.
# I will simulate on grid, but only save middle voxel (no smoothing).

# Then estimate the variance using FSL mixed effects + R lmer.
# Then construct CI and calculate empirical coverage (EC) of the CI 
# around the true parameter. 


##
###############
### Preparation
###############
##

# Let us just run it locally for now
for(l in 1:200){
print(l)
# Take argument from master file
input <- commandArgs(TRUE)
# K'th simulation
hpcID <- try(as.numeric(as.character(input)[1]),silent=TRUE)
# Which machine
MACHINE <- try(as.character(input)[2],silent=TRUE)

# DataWrite directory: where all temporary files are written to
DataWrite <- try(as.character(input)[3],silent=TRUE)

# If no machine is specified, then it has to be this machine!
if(is.na(MACHINE)){
  MACHINE <- 'MAC'
  hpcID <- l
  DataWrite <- '~/Desktop/VAR2LVL'
}

# Give path to FSL
if(MACHINE=='HPC'){
  fslpath <- ''
}
if(MACHINE=='MAC'){
  fslpath <- '/usr/local/fsl/bin/'
}

# Implement for loop over r iterations here: hpcID goes from 1 to 100 in master file
rIter <- 10
startIndex <- try(1 + rIter * (hpcID - 1), silent = TRUE)
endIndex <- try(startIndex + (rIter - 1), silent = TRUE)


# Set WD: this is location where results are written
if(MACHINE=='HPC'){
  wd <- '/user/scratch/gent/gvo000/gvo00022/vsc40728/Variance_2lvl/'
}
if(MACHINE=='MAC'){
  wd <- '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/Variance_2lvl/'
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


###################
#### Global variables
###################

# Number of subjects
nsub <- 50

# Value for sigma in the model
sigma_eps <- 100

# Between subject variability (variance of random slope)
sigma_b2 <- c(0, 5, 10)[3]

# Variance of random intercept
sigma_b1 <- c(0, 1, 1)[3]


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
intcpt <- 100

#######################################
#### DESIGN AND SIGNAL TIME SERIES ####
#######################################

# Generating a design matrix: convolution of block design with double-gamma HRF
X <- neuRosim::simprepTemporal(total,1,onsets = onsets,
                               effectsize = 1, durations = duration,
                               TR = TR, acc = 0.1, hrf = "double-gamma")

# X vector for one subject = predicted signal
X_s <- neuRosim::simTSfmri(design=X, base=0, SNR=1, noise="none", verbose=FALSE)

# Now the model will be: (intcpt + b1) + (BOLDC + b2) * pred + epsilon

## Design parameters
# Extend the design matrix with the intercept
xIN <- cbind(intcpt, X_s)

# Contrast: not interested in intercept
CONTRAST <- matrix(c(0,1),nrow=1)

# Calculate (X'X)^(-1) with contrast
design_factor <- CONTRAST %*% (solve(t(xIN) %*% xIN )) %*% t(CONTRAST)


##################
#### GENERATE DATA
##################

# Empty lmer results data frame
LMER_res <- FLAME_res <- comb_res <- data.frame() %>% as_tibble()

# Start some iterations (increases efficiency since iterations run very fast)
for(ID in startIndex:endIndex){
  # Set starting seed
  starting.seed <- pi*ID
  set.seed(starting.seed)

  # Generate D matrix: variance-covariance matrix of random intercept + slope
  # Variance of slope = sigma_b**2
  var_cov_D <- rbind(c(sigma_b1**2, 0), c(0, sigma_b2**2))
  
  # Generate the subject-specific values for intercept and slope using this D-matrix
  B_matrix <- MASS::mvrnorm(nsub, mu=c(0,0), Sigma = var_cov_D)
  
  # Empty vector
  Y <- data.frame() %>% as_tibble()
  
  # For loop over all subjects
  for(i in 1:nsub){
    # Generate nscan values, corresponding to time series of one subject 
      # note: random intercept and random slope generated earlier
    Y_s <- (intcpt + B_matrix[i,1]) + ((BOLDC + B_matrix[i,2]) * X_s) + 
      rnorm(n = nscan, mean = 0, sd = sigma_eps)
    
    # Add to data frame
    Y <- data.frame(Y = Y_s, X = X_s, sub = as.integer(i)) %>% as_tibble() %>%
      bind_rows(Y, .)
  }
  
  #############################################
  #### LINEAR MIXED MODEL APPROACH USING R ####
  #############################################
  
  # Fit model with random intercept and random slope for subject.
  # Get coefficients using tidy.
  LMER_results <- broom::tidy(lmer(Y ~ 1 + X + (1 + X|sub), data = Y)) %>%
    as_tibble() %>% mutate(sim = ID) 
  
  #################################################
  #### LINEAR MIXED MODEL APPROACH USING WATER ####
  #################################################
  
  # For this, we need to first analyze each subject individually, save COPE and VARCOPE
  # and then proceed.
  # We call this object secLevel
  secLevel <- Y %>% 
    group_by(sub) %>%
    do(., 
       # For each subject, fit linear model with an intercept and X as predictors
       broom::tidy( 
         lm(Y ~ 1 + X, data = .))) %>%
    # Filter on predictor
    filter(term == 'X') %>%
    # Now select the estimate and standard error
    dplyr::select(sub, estimate, std.error) %>%
    # Create variance
    mutate(varCope = std.error^2)
  
  # Create 4D images (all voxels in first 3 dimensions are the same), otherwise FSL crashes!
  # Then convert the estimates and variance to nifti images
  COPE4D <- nifti(img=array(rep(as.numeric(secLevel$estimate), each = 8), 
                            dim=c(2,2,2,nsub)),
                  dim=c(2,2,2,nsub), datatype = 16)
  VARCOPE4D <- nifti(img=array(rep(as.numeric(secLevel$varCope), each = 8),
                               dim=c(2,2,2,nsub)), 
                     dim=c(2,2,2,nsub), datatype = 16)
  
  # Write them to DataWrite
  writeNIfTI(COPE4D, filename = paste(DataWrite,'/COPE',sep=''), gzipped=FALSE)
  writeNIfTI(VARCOPE4D, filename = paste(DataWrite,'/VARCOPE',sep=''), gzipped=FALSE)

  # Write auxiliarly files to DataWrite. We need:
  # GRCOPE in nifti
  # GRVARCOPE in nifti
  # 4D mask
  # design.mat file
  # design.grp file
  # design.con file
  
  #----- 1 ----#
  ### Design.mat
  fileCon <- paste(DataWrite,"/design.mat",sep="")
  # Text to be written to the file
  cat('/NumWaves\t1
      /NumPoints\t',paste(nsub,sep=''),'
      /PPheights\t\t1.000000e+00
      
      /Matrix
      ',rep("1.000000e+00\n",nsub),file=fileCon)
  
  #----- 2 ----#
  ### Design.con
  fileCon <- file(paste(DataWrite,"/design.con", sep=""))
  writeLines('/ContrastName1	Group Average
             /NumWaves	1
             /NumContrasts	1
             /PPheights		1.000000e+00
             /RequiredEffect		5.034
             
             /Matrix
             1.000000e+00
             ',fileCon)
  close(fileCon)
  
  #----- 3 ----#
  ### Design.grp
  fileCon <- paste(DataWrite,"/design.grp",sep="")
  # Text to be written to the file
  cat('/NumWaves\t1
      /NumPoints\t',paste(nsub,sep=''),'
      
      /Matrix
      ',rep("1\n",nsub),file=fileCon)

  #----- 4 ----#
  ### mask.nii
  mask <- nifti(img=array(1, dim=c(2,2,2,nsub)), dim=c(2,2,2,nsub), datatype=2)
  writeNIfTI(mask, filename = paste(DataWrite,'/mask',sep=''), gzipped=FALSE)
  
  # FSL TIME!
  setwd(DataWrite)
  command <- paste(fslpath, 'flameo --cope=COPE --vc=VARCOPE --mask=mask --ld=FSL_stats --dm=design.mat --cs=design.grp --tc=design.con --runmode=flame1', sep='')
  Sys.setenv(FSLOUTPUTTYPE="NIFTI")
  system(command)
  
  # Read back results
  FLAME_results <- data.frame(value = c(
      readNIfTI(paste(DataWrite,"/FSL_stats/cope1.nii",sep=""), 
                         verbose=FALSE, warn=-1, reorient=TRUE, 
                         call=NULL)[1,1,1],
              readNIfTI(paste(DataWrite,"/FSL_stats/varcope1.nii",sep=""), 
                        verbose=FALSE, warn=-1, reorient=TRUE, 
                        call=NULL)[1,1,1])) %>%
    mutate(parameter = c('estimate', 'variance'))
  
  # Degrees of freedom:
  tdof_t1 <- readNIfTI(paste(DataWrite,"/FSL_stats/tdof_t1.nii",sep=""), 
                       verbose=FALSE, warn=-1, reorient=TRUE, 
                       call=NULL)[1,1,1]
  
  # The estimated between-subject variability
  var_bsub <- readNIfTI(paste(DataWrite,"/FSL_stats/mean_random_effects_var1.nii",sep=""), 
              verbose=FALSE, warn=-1, reorient=TRUE, 
              call=NULL)[1,1,1]
  
  ############################################################
  #### CONSRUCT 95% CONFIDENCE INTERVALS AND CALCULATE EC ####
  ############################################################
  
  LMER_res <-
    LMER_results %>% filter(term == 'X') %>%
    dplyr::select(term, estimate, std.error) %>%
    # CI around beta: using std.error of parameter!
    mutate(CIlow = estimate - qt(0.975, df = tdof_t1) * std.error,
           CIup = estimate + qt(0.975, df = tdof_t1) * std.error) %>%
    mutate(EC = ifelse(BOLDC >= CIlow & BOLDC <= CIup, 1, 0)) %>%
    # Now select the estimate of between-subject variability (SD)
    mutate(sd_X.sub = unlist(LMER_results %>% filter(term == 'sd_X.sub') %>%
             dplyr::select(estimate))) %>%
    # Add variance of parameter estimate (VARCOPE)
    mutate(variance = std.error^2) %>%
    # re-arrange
    dplyr::select(estimate, std.error, variance, sd_X.sub, CIlow, CIup, EC) %>%
    # Rename
    rename(estimate = estimate, SE_beta = std.error, variance = variance,
           SD_bsub = sd_X.sub, CIlow = CIlow, CIup = CIup, EC = EC) %>%
    mutate(type = 'LMER', simID = ID)
  

  FLAME_res <- FLAME_results %>% 
    tidyr::spread(key = parameter, value = value) %>% 
    mutate(CIlow = estimate - qt(0.975, df = tdof_t1) * sqrt(variance),
           CIup = estimate + qt(0.975, df = tdof_t1) * sqrt(variance)) %>%
    mutate(EC = ifelse(BOLDC >= CIlow & BOLDC <= CIup, 1, 0)) %>%
    # Add info and rename data object
    mutate(type = 'FLAME', simID = ID, SE_beta = sqrt(variance),
           SD_bsub = sqrt(var_bsub)) %>%
    # Re-order
    dplyr::select(estimate, SE_beta, variance, SD_bsub, CIlow, CIup,
                  EC, type, simID) %>% as_tibble()
  
  #########################################
  #### COMBINE DATA AND WRITE TO FILES ####
  #########################################
  
  comb_res <- bind_rows(comb_res, LMER_res, FLAME_res)
  
  # Remove objects in DataWrite folder
  command <- paste0('rm -r ', DataWrite, '/*')
  system(command)
}


# Save R object
saveRDS(comb_res, file = paste0(wd, 'Results_bsub_',sigma_b2,'/VAR2LVL_', hpcID, '.rda'))

# Reset
rm(list = ls())

} 
    
    

  


  
  
  