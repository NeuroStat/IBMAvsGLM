####################
#### TITLE:   Variance estimation at third (studies) level + CI coverage
#### Contents:
####
#### Source Files:
#### First Modified: 21/06/2018
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
# Then pass these estimates to third level.
# Then construct CI and calculate empirical coverage (EC) of the CI 
# around the true parameter. 


##
###############
### Preparation
###############
##


# Take argument from master file
input <- commandArgs(TRUE)
# K'th simulation
hpcID <- try(as.numeric(as.character(input)[1]),silent=TRUE)
# Which machine
MACHINE <- try(as.character(input)[2],silent=TRUE)

# DataWrite directory: where all temporary files are written to
DataWrite <- try(as.character(input)[3],silent=TRUE)

# Select the amount of variability (between-subject and between-study)
inptVar <- try(as.numeric(as.character(input)[4]),silent=TRUE)

# If no machine is specified, then it has to be this machine!
if(is.na(MACHINE)){
  MACHINE <- 'MAC'
  hpcID <- 104
  DataWrite <- '~/Desktop/VAR3LVL'
  inptVar <- 3
}

# Give path to FSL
if(MACHINE=='HPC'){
  fslpath <- ''
}
if(MACHINE=='MAC'){
  fslpath <- '/usr/local/fsl/bin/'
}

# Implement for loop over r iterations here: hpcID goes from 1 to 200 in master file
rIter <- 10
startIndex <- try(1 + rIter * (hpcID - 1), silent = TRUE)
endIndex <- try(startIndex + (rIter - 1), silent = TRUE)


# Set WD: this is location where results are written
if(MACHINE=='HPC'){
  wd <- '/user/scratch/gent/gvo000/gvo00022/vsc40728/Variance_3lvl/Results/'
}
if(MACHINE=='MAC'){
  wd <- '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/Variance_3lvl/'
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

# Number of studies
nstud <- 40
  
# Number of subjects
nsub <- 50

# Value for sigma in the model
sigma_eps <- 100

# Between subject variability (variance of random slope)
sigma_b2 <- c(0, 5, 10)[inptVar]

# Variance of random intercept
sigma_b1 <- c(0, 1, 1)[inptVar]

# Between study variability: intercept
sigma_eta1 <- c(0, 1, 1)[inptVar]

# Between study variability: slope
sigma_eta2 <- c(0, 2, 4)[inptVar]

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
  
  # Generate var-covar matrix of studies: random intercept + slope
  var_cov_Study <- rbind(c(sigma_eta1**2, 0), c(0, sigma_eta2**2))
  
  # Now generate the study specific random intercept and slope values using this matrix
  Stud_matrix <- MASS::mvrnorm(nstud, mu=c(0,0), Sigma = var_cov_Study)
  
  # Empty vector: studies
  Y_stud <- thirLevel <- data.frame() %>% as_tibble()
  
  # For loop over the studies
  for(s in 1:nstud){
    
    # Empty vector: subjects
    Y_sub <- data.frame() %>% as_tibble()

    # Generate D matrix (subject): variance-covariance matrix of random intercept + slope
    # Variance of slope = sigma_b**2
    var_cov_D <- rbind(c(sigma_b1**2, 0), c(0, sigma_b2**2))
    
    # Generate the subject-specific values for intercept and slope using this D-matrix
    B_matrix <- MASS::mvrnorm(nsub, mu=c(0,0), Sigma = var_cov_D)
    
    # For loop over all subjects
    for(i in 1:nsub){
      # Generate nscan values, corresponding to time series of one subject 
      # note: random intercepts and random slopes generated earlier
      Y_s <- (intcpt + Stud_matrix[s,1] + B_matrix[i,1]) + 
              ((BOLDC + Stud_matrix[s,2] + B_matrix[i,2]) * X_s) + 
              rnorm(n = nscan, mean = 0, sd = sigma_eps)
      
      # Add to data frame
      Y_sub <- data.frame(Y = Y_s, X = X_s, sub = as.integer(i)) %>% as_tibble() %>%
        bind_rows(Y_sub, .)
    }
    
    ################################################
    #### LINEAR MIXED MODEL APPROACH USING LMER ####
    ################################################
    
    # For now, we add the data to one big data frame, which will be analyzed later
    Y_stud <- data.frame(Y_sub, stud = as.integer(s)) %>% as_tibble() %>%
      bind_rows(Y_stud, .)
    
    
    #################################################
    #### LINEAR MIXED MODEL APPROACH USING WATER ####
    #################################################
    
    # For this, we need to first analyze each subject individually, save COPE and VARCOPE
    # and then proceed.
    # We call this object secLevel
    secLevel <- Y_sub %>% 
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
    FLAME_results_2lvl <- data.frame(value = c(
      readNIfTI(paste(DataWrite,"/FSL_stats/cope1.nii",sep=""), 
                verbose=FALSE, warn=-1, reorient=TRUE, 
                call=NULL)[1,1,1],
      readNIfTI(paste(DataWrite,"/FSL_stats/varcope1.nii",sep=""), 
                verbose=FALSE, warn=-1, reorient=TRUE, 
                call=NULL)[1,1,1])) %>%
      mutate(parameter = c('estimate', 'variance'))

    # Add to data frame
    thirLevel <- tidyr::spread(FLAME_results_2lvl, key = parameter, value = value) %>%
      mutate(study = s) %>%
      dplyr::select(study, estimate, variance) %>%
      rename(study = study, estimate = estimate, varCope = variance) %>% 
      as.tibble() %>%
      bind_rows(thirLevel, .)
    
    # Remove objects in DataWrite folder
    command <- paste0('rm -r ', DataWrite, '/*')
    system(command)
  }
  
  #############################################
  #### LINEAR MIXED MODEL APPROACH USING R ####
  #############################################
  
  # Fit model with random intercept and random slope for subject and study.
    # Subjects are nested within studies!
  # Get coefficients using tidy.
    # Note: no need to make explicit factors for subject and study.
  LMER_results <- broom::tidy(lmer(Y ~ 1 + X + (1 + X|stud/sub), data = Y_stud)) %>%
    as_tibble() %>% mutate(sim = ID) 
      
  ##########################################################
  #### LEVEL 3: LINEAR MIXED MODEL APPROACH USING FLAME ####
  ##########################################################
  
  # Create 4D images (all voxels in first 3 dimensions are the same), otherwise FSL crashes!
  # Then convert the estimates and variance to nifti images
  COPE4D <- nifti(img=array(rep(as.numeric(thirLevel$estimate), each = 8), 
                            dim=c(2,2,2, nstud)),
                  dim=c(2,2,2,nstud), datatype = 16)
  VARCOPE4D <- nifti(img=array(rep(as.numeric(thirLevel$varCope), each = 8),
                               dim=c(2,2,2, nstud)), 
                     dim=c(2,2,2,nstud), datatype = 16)
  
  # Write them to DataWrite
  writeNIfTI(COPE4D, filename = paste(DataWrite,'/COPE_S',sep=''), gzipped=FALSE)
  writeNIfTI(VARCOPE4D, filename = paste(DataWrite,'/VARCOPE_S',sep=''), gzipped=FALSE)
  
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
      /NumPoints\t',paste(nstud,sep=''),'
      /PPheights\t\t1.000000e+00
      
      /Matrix
      ',rep("1.000000e+00\n",nstud),file=fileCon)
  
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
      /NumPoints\t',paste(nstud,sep=''),'
      
      /Matrix
      ',rep("1\n",nstud),file=fileCon)
  
  #----- 4 ----#
  ### mask.nii
  mask <- nifti(img=array(1, dim=c(2,2,2,nstud)), dim=c(2,2,2,nstud), datatype=2)
  writeNIfTI(mask, filename = paste(DataWrite,'/mask',sep=''), gzipped=FALSE)
  
  # FSL TIME!
  setwd(DataWrite)
  command <- paste(fslpath, 'flameo --cope=COPE_S --vc=VARCOPE_S --mask=mask --ld=FSL_stats --dm=design.mat --cs=design.grp --tc=design.con --runmode=flame1', sep='')
  Sys.setenv(FSLOUTPUTTYPE="NIFTI")
  system(command)
  
  # Read back results
  FLAME_results_3lvl <- data.frame(value = c(
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
  
  # The estimated between-study variability
  var_bstud <- readNIfTI(paste(DataWrite,"/FSL_stats/mean_random_effects_var1.nii",sep=""), 
                        verbose=FALSE, warn=-1, reorient=TRUE, 
                        call=NULL)[1,1,1]
  
  
  ############################################################
  #### CONSRUCT 95% CONFIDENCE INTERVALS AND CALCULATE EC ####
  ############################################################
  
  LMER_res <- LMER_results %>% filter(term == 'X') %>%
    dplyr::select(term, estimate, std.error) %>%
    # CI around beta: using std.error of parameter!
    mutate(CIlow = estimate - qt(0.975, df = tdof_t1) * std.error,
           CIup = estimate + qt(0.975, df = tdof_t1) * std.error) %>%
    mutate(EC = ifelse(BOLDC >= CIlow & BOLDC <= CIup, 1, 0)) %>%
    # Now select the estimate of between-study variability (SD)
    mutate(sd_X.stud = unlist(LMER_results %>% filter(term == 'sd_X.stud') %>%
                               dplyr::select(estimate))) %>%
    # Add variance of parameter estimate (VARCOPE)
    mutate(variance = std.error^2) %>%
    # re-arrange
    dplyr::select(estimate, std.error, variance, sd_X.stud, CIlow, CIup, EC) %>%
    # Rename
    rename(estimate = estimate, SE_beta = std.error, variance = variance,
           SD_bstud = sd_X.stud, CIlow = CIlow, CIup = CIup, EC = EC) %>%
    mutate(type = 'LMER', simID = ID)

  FLAME_res <- FLAME_results_3lvl %>% 
    tidyr::spread(key = parameter, value = value) %>% 
    mutate(CIlow = estimate - qt(0.975, df = tdof_t1) * sqrt(variance),
           CIup = estimate + qt(0.975, df = tdof_t1) * sqrt(variance)) %>%
    mutate(EC = ifelse(BOLDC >= CIlow & BOLDC <= CIup, 1, 0)) %>%
    # Add info and rename data object
    mutate(type = 'FLAME', simID = ID, SE_beta = sqrt(variance),
           SD_bstud = sqrt(var_bstud)) %>%
    # Re-order
    dplyr::select(estimate, SE_beta, variance, SD_bstud, CIlow, CIup,
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
saveRDS(comb_res, file = paste0(wd, 'Results_bsub_', sigma_b2, '_bstud_', sigma_eta2, '/VAR3LVL_', hpcID, '.rda'))



 








