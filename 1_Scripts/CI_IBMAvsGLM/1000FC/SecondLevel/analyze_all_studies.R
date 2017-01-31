####################
#### TITLE:     MA vs IBMA using FSL's FLAME on the COPE and VARCOPES.
#### Contents:
####
#### Source Files:
#### First Modified: 24/02/2016
#### Notes:
#################



##
###############
### Notes
###############
##

# In this scenario, we will compare the outcome of transforming the second level GLM to an ES and execute a meta-analsysis with the scenario in which you use a third level GLM.
# MEASURES:
#   * CI coverage
#   * Standardized bias
#   * Average CI length



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
  # Number of simulation
  simulID <- try(as.numeric(as.character(input)[1]),silent=TRUE)
  # Number of subjects
  NSUB <- try(as.numeric(as.character(input)[2]),silent=TRUE)
  # Number of studies
  NSTUD <- try(as.numeric(as.character(input)[3]),silent=TRUE)
  # Scanning site of the 1000FC
  SITE <- try(as.character(input)[4],silent=TRUE)
  # Which machine
  MACHINE <- try(as.character(input)[5],silent=TRUE)
    # If no machine is specified, then it has to be this machine!
    if(is.na(MACHINE)){
      MACHINE <- 'MAC'
      simulID <- 1
      NSUB <- 20
      NSTUD <- 5
      SITE <- 'Cambridge'
    }
  # DataWrite directory: where all files are written to
  #DataWrite <- try(as.character(input)[5],silent=TRUE)

# Libraries
library(dplyr)

# Set WD
if(MACHINE=='HPC'){
  wd <- paste('/user/scratch/gent/gvo000/gvo00022/vsc40728/1000FC/', SITE, '/SecondLevel', sep = '')
}
if(MACHINE=='MAC'){
  wd <- '/Users/hanbossier/Dropbox/PhD/PhDWork/Meta Analysis/R Code/Studie_Simulation/SimulationGit/1_Scripts/CI_IBMAvsGLM/1000FC/SecondLevel'
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
# if(MACHINE=='MAC'){
#   DataWrite <- '~/Desktop/IBMA2'
# }



# Read in the template design file, which we will alter
if(MACHINE=='HPC'){
  Base_Design <- readLines(con = paste(wd, '/Design_templates/GLMSecondLevel.fsf', sep = ''))
}
if(MACHINE=='MAC'){
  Base_Design <- readLines(con = paste(wd, '/Design_templates/GLMSecondLevel.fsf', sep = ''))
}

# Location of subject IDs
if(MACHINE=='HPC'){
  LocSubjID <- paste(wd, '/IDs/', simulID, sep = '')
}
if(MACHINE=='MAC'){
  LocSubjID <- paste('/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IBMAvsGLM/Results/', SITE, '/SecondLevel/IDs/', simulID, sep = '')
}


##
###############
### Analysis parameters
###############
##

#NSUB <- 20
#NSTUD <- 5

# List of smoothing parameters
SMOOTHINGs <- list(
  '4' = '4mm',
  '8' = '8mm'
  )
SMOOTHING <- SMOOTHINGs[['8']]

# List of design parameters
DESIGNs <- list(
  '10' = 'boxcar10',
  '30' = 'boxcar30'
  )
DESIGN <- DESIGNs[['10']]

# Subject IDs based on study
subjIDs <-

studyID <- 1

# ----

# New output directory
#!!! Base_Output <- paste('/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IBMAvsGLM/Results/',SITE,'/SecondLevel/', SMOOTHING, '/', DESIGN, '/', sep = '')

# TR (depending on scanning site)
TRs <- list(
  'Cambridge' = 3,
  'Beijing' = 2
  )
TR <- TRs[[SITE]]

# Number of volumes (= NSUB)
npts <- NSUB

# Number of first-level analyses (= NSUB)
multiple <- NSUB

# Directory of input data
#!!! Base_Input <- paste('/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IBMAvsGLM/Results/',SITE,'/',SMOOTHING,'/',DESIGN,'/', sep = '')

# Higher-level EV value for EV

# Group membership


##
###############
### Pooling subjects
###############
##

# Start for loop over NSTUD
for(k in 1:NSTUD){
  # Read in the subject IDs for this study
  SubjStudID <-
    read.table(file = paste(LocSubjID, '/subjects_STUD_', k,'.txt', sep = ''), header = FALSE, sep = '\n', stringsAsFactors = FALSE) %>%
      as.matrix(., ncol = 1)

  # Let's change analysis options in the base design file
    # Output directory
    output <- paste(Base_Output, simulID,'/', 'NSTUD_', studyID, sep = '')
    Base_Design[grepl('outputdir', Base_Design)] <- paste('set fmri(outputdir) "', output, '"', sep = '')

    # TR
    Base_Design[grepl('fmri(tr)', Base_Design, fixed = TRUE)] <-
        paste('set fmri(tr) ', TR, sep = '')

    # number of volumes and first-level analyses
    Base_Design[grepl('npts', Base_Design, fixed = TRUE)] <-
        paste('set fmri(npts) ', npts, sep ='')
    Base_Design[grepl('multiple', Base_Design, fixed = TRUE)] <-
        paste('set fmri(multiple) ', multiple, sep ='')

    # Directories of input data
      # Start with getting first and last line
      BeginInput <- which(grepl('4D AVW data or FEAT directory (1)', Base_Design, fixed = TRUE)) - 1
      EndInput <- which(grepl('4D AVW data or FEAT directory (5)', Base_Design, fixed = TRUE)) + 3

      # Now create new lines, based on number of subjects and their IDs
      InputLines <- c()
      for(s in 1:NSUB){
        InputLines <- rbind(InputLines, rbind(
          paste('# 4D AVW data or FEAT directory (',s,') ', sep = ''),
          paste('set feat_files(',s,') "', Base_Input, SubjStudID[s],'/results.feat"', sep = ''),
          c(""))
        )
      }

      # Add this between first and last line of Base_Design
      Base_Design <- rbind(matrix(Base_Design[c(1:BeginInput)], ncol = 1),
          InputLines,
          matrix(Base_Design[c(EndInput:length(Base_Design))], ncol = 1)
          ) %>% as.character()

    # Higher-level EV value for EV
      # Start and end of EV value
      BeginEV <- which(grepl('Higher-level EV value for EV 1 and input 1', Base_Design, fixed = TRUE)) - 1
      EndEV <- which(grepl('Higher-level EV value for EV 1 and input 5', Base_Design, fixed = TRUE)) + 3

      # Create lines
      EVLines <- c()
      for(s in 1:NSUB){
        EVLines <- rbind(EVLines, rbind(
          paste('# Higher-level EV value for EV 1 and input ',s, sep = ''),
          paste('set fmri(evg',s,'.1) 1', sep = ''),
          c(""))
        )
      }

      # Add this between first and last line of Base_Design
      Base_Design <- rbind(matrix(Base_Design[c(1:BeginEV)], ncol = 1),
          EVLines,
          matrix(Base_Design[c(EndEV:length(Base_Design))], ncol = 1)
          ) %>% as.character()

    # Group membership
      # Start and end of group membership value
      BeginGM <- which(grepl('Group membership for input 1', Base_Design, fixed = TRUE)) - 1
      EndGM <- which(grepl('Group membership for input 5', Base_Design, fixed = TRUE)) + 3

      # Create lines
      GMLines <- c()
      for(s in 1:NSUB){
        GMLines <- rbind(GMLines, rbind(
          paste('# Group membership for input ',s, sep = ''),
          paste('set fmri(groupmem.',s,') 1', sep = ''),
          c(""))
        )
      }

      # Add this between first and last line of Base_Design
      Base_Design <- rbind(matrix(Base_Design[c(1:BeginGM)], ncol = 1),
          GMLines,
          matrix(Base_Design[c(EndGM:length(Base_Design))], ncol = 1)
          ) %>% as.character()





s<-1
subjIDs <- c(1:20)




writeLines(text = test, con = '/Users/hanbossier/Dropbox/PhD/PhDWork/Meta Analysis/R Code/Studie_Simulation/SimulationGit/1_Scripts/CI_IBMAvsGLM/1000FC/FirstLevel/Design_templates/GLMCambridgeTest.fsf')







