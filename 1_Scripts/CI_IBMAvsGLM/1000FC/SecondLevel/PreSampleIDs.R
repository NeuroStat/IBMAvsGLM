####################
#### TITLE:     Pre-sample subject IDs into studies: 1000FC part of hierarchical study.
#### Contents:
####
#### Source Files:
#### First Modified: 31/01/2017
#### Notes:
#################



##
###############
### Notes
###############
##

# In this script, we sample each time 35 subjects from a scanning site from the 1000 functional connectomes project into 5 studies.
# This is done for 3000 times (amount of simulations in the simulations part of the study).

# We can then run the analyses in parallel on the HPC.


##
###############
### Preparation
###############
##


# Reset working memory
rm(list=ls())
gc(verbose = FALSE)

# Libraries
library(dplyr)

# Working directory
wd <- '/Users/hanbossier/Dropbox/PhD/PhDWork/Meta Analysis/R Code/Studie_Simulation/SimulationGit/1_Scripts/CI_IBMAvsGLM/1000FC/SecondLevel'

# Writing directory
WriteDir <- '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IBMAvsGLM/Results/Cambridge/SecondLevel/IDs'

# Set seed
seed <- 1
set.seed(seed)

# Number of simulations
nsim <- 3000

# Number of subjects per study
NSUB <- 20

# Number of studies per MA
NSTUD <- 5

##
###############
### Read, sample and write
###############
##


# First read in all subject IDs
IDs <- read.table(file = '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/FCON_1000/Cambridge/Cambridge_Buckner_subjects.txt', header = FALSE, sep = '\n')

# Now for loop over simulations: sample and write
for(i in 1:nsim){
  # Setting seed
  seed_new <- i * (NSUB * NSTUD)
  set.seed(seed_new)

  # Sample IDs: over all studies
  subjIter <- data.frame('IDs' = sample(x = IDs[,1], size = (NSUB * NSTUD)))

  # Now split them in studies and write to folder
  subjIter$study <- rep(1:NSTUD, each = NSUB)
  for(j in 1:NSTUD){
    filter(subjIter, study == j) %>% select(IDs) %>%
      write.table(., file = paste(WriteDir, '/', i, '/subjects_STUD_', j, '.txt', sep = ''),
                                      sep = '\n', row.names = FALSE,
                                      col.names = FALSE, quote = FALSE)
  }
}




































