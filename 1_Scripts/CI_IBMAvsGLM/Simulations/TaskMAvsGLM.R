####################
#### TITLE:     MA vs Mixed Model GLM using FSL's FLAME on activation data.
#### Contents:
####
#### Source Files:
#### First Modified: 03/05/2017
#### Notes:
#################



##
###############
### Notes
###############
##

# Now simulate activation data with realistic noise, sample sizes, etc, ...
# In this scenario, we will compare the outcome of transforming the second level GLM to an ES and execute a meta-analsysis with the scenario in which you use a third level GLM.
# MEASURES:
#   * CI coverage
#   * Standardized bias


##
###############
### Preparation
###############
##

# Reset working directory
rm(list=ls())
gc(verbose = FALSE)

# Take argument from master file
input <- commandArgs(TRUE)
  # K'th simulation
  K <- try(as.numeric(as.character(input)[1]),silent=TRUE)
  # Which scenario
  SCEN <- try(as.numeric(as.character(input)[2]),silent=TRUE)
  # Which machine
  MACHINE <- try(as.character(input)[3],silent=TRUE)
    # If no machine is specified, then it has to be this machine!
    if(is.na(MACHINE)){
      MACHINE <- 'MAC'
      K <- 1
      SCEN <- 1
    }
  # DataWrite directory: where all files are written to
  DataWrite <- try(as.character(input)[4],silent=TRUE)

# Set starting seed: it is the product of the amount of voxels, the number of studies and the number of subjects!
starting.seed <- 36865*K
set.seed(starting.seed)


# Set WD
if(MACHINE=='HPC'){
  wd <- '/user/scratch/gent/gvo000/gvo00022/vsc40728/IBMAvsMA'
}
if(MACHINE=='MAC'){
  wd <- '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/MAvsIBMA'
}
setwd(wd)

# Give path to FSL
if(MACHINE=='HPC'){
  fslpath <- ''
}
if(MACHINE=='MAC'){
  fslpath <- '/usr/local/fsl/bin/'
}

# DataWrite if machine = Mac
if(MACHINE=='MAC'){
  DataWrite <- '~/Desktop/IBMAtmp'
}


# Load in libraries
library(AnalyzeFMRI)
library(fmri)
library(lattice)
library(gridExtra)
library(oro.nifti)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(Hmisc)
library(devtools)
library(neuRosim)


# Load in functions from FixRan study: THIS HAS TO COME AFTER ALL LIBRARIES ARE LOADED AS WE SOMETIMES FIX FUNCTIONS THAT ARE BUGGED IN THE PACKAGES
if(MACHINE == 'MAC'){
  source('~/Dropbox/PhD/PhDWork/Meta\ Analysis/R\ Code/Studie_FixRan/FixRanStudyGit.git/Development/functions.R')
}
if(MACHINE == 'HPC'){
  source('/user/scratch/gent/gvo000/gvo00022/vsc40728/Simulation/functions.R')
}



##
###############
### Simulation steps
###############
##

####************####
#### Global options
####************####
TR <- 2
nscan <- 200
total <- TR*nscan
on1 <- seq(1,total,40)
onsets <- list(on1)
duration <- list(20)
effect <- list(1) 			                            ## Effect of 1 for designmatrix
DIM <- c(4,4,4)


####************####
#### Scenario specific simulation details
####************####

# Let us randomly draw sample sizes for the individual studies from the Poldrack et al. paper 
# For this, source the data from David et al., as well as the added data 
# scraping Neurosynth

"https://raw.githubusercontent.com/poldracklab/ScanningTheHorizon/master/SampleSize/neurosynth_study_data.txt"

"https://raw.githubusercontent.com/poldracklab/ScanningTheHorizon/master/SampleSize/david_sampsizedata.txt"


