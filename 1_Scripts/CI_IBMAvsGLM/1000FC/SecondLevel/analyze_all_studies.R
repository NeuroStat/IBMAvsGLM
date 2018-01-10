####################
#### TITLE:     IBMA versus GLM: resting state fMRI (1000FC) second level.
#### Contents:
####
#### Source Files: ~\\Studie_Simulation/SimulationGit/1_Scripts/CI_IBMAvsGLM/1000FC/SecondLevel
#### First Modified: 31/01/2017
#### Notes:
#################



##
###############
### Notes
###############
##


# After running the First Level analyze_all_subjects_fsl.sh, we end up with
# pre-processed subjects analyzed using different paradigms.

# Here we take a design template for a second level group analysis from FSL (obtained by trying out a second level analyis manually).
# Then we modify parameters in this design file, combine subjects into studies and run this .fsf analysis.


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
  # Smoothing parameter
  SMOOTHING <- try(as.character(input)[6],silent=TRUE)
  # Design parameter
  DESIGN <- try(as.character(input)[7],silent=TRUE)
  # Location of subject IDs in this simulation
  LocSubjID <- try(as.character(input)[8],silent=TRUE)
  # Level
  LEVEL <- try(as.character(input)[9],silent=TRUE)
    # If no machine is specified, then it has to be this machine in which we are testing code!
    if(is.na(MACHINE)){
      MACHINE <- 'MAC'
      simulID <- 1
      NSUB <- 20
      NSTUD <- 5
      SITE <- 'Cambridge'
      SMOOTHING <- '8mm'
      DESIGN <- 'boxcar10'
      LEVEL <- 'SecondLevel'
      # Location of IDs of subjects which we sample
      LocSubjIDBOOLEAN <- TRUE
    }

# Libraries
library(dplyr)

# Set WD and base data directory of data processed at first level
if(MACHINE=='HPC'){
  wd <- paste('/user/scratch/gent/gvo000/gvo00022/vsc40728/1000FC/',LEVEL , sep = '')
  data_firstLevel_dir <- paste('/user/data/gent/gvo000/gvo00022/vsc40728', sep = '')
}
if(MACHINE=='MAC'){
  wd <- paste('/Users/hanbossier/Dropbox/PhD/PhDWork/Meta Analysis/R Code/Studie_Simulation/SimulationGit/1_Scripts/CI_IBMAvsGLM/1000FC/', LEVEL, sep = '')
  data_firstLevel_dir <- '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IBMAvsGLM/Results'
  if(isTRUE(LocSubjIDBOOLEAN)) LocSubjID <- paste(data_firstLevel_dir, '/', SITE, '/', LEVEL ,'/IDs/', simulID, sep = '')
}
setwd(wd)

# Give path to FSL
if(MACHINE=='HPC'){
  fslpath <- ''
}
if(MACHINE=='MAC'){
  fslpath <- '/usr/local/fsl/bin/'
}

# Read in the template design file, which we will modify
if(MACHINE=='HPC'){
  Base_Design <- readLines(con = paste(wd, '/Design_templates/GLMSecondLevel.fsf', sep = ''))
}
if(MACHINE=='MAC'){
  Base_Design <- readLines(con = paste(wd, '/Design_templates/GLMSecondLevel.fsf', sep = ''))
}


##
###############
### Analysis parameters
###############
##

# Define the new output directory
if(MACHINE=='HPC'){
  Base_Output <- paste(wd, '/Results/', SITE, '/', SMOOTHING, '/', DESIGN, '/', sep = '')
}
if(MACHINE=='MAC'){
  Base_Output <- paste(data_firstLevel_dir, '/',SITE,'/', LEVEL, '/', SMOOTHING, '/', DESIGN, '/', sep = '')
}

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
  # In HPC version: first level analyses are copied next to the ouptut directory to resolve interference.
if(MACHINE=='HPC'){
  Base_Input <- paste(Base_Output, simulID,'/NSTUD_', sep = '')
}
if(MACHINE=='MAC'){
  Base_Input <- paste(data_firstLevel_dir, '/',SITE,'/',SMOOTHING,'/',DESIGN,'/', sep = '')
}


##
###############
### Pooling subjects: create design files
###############
##

# Start for loop over NSTUD
for(k in 1:NSTUD){
  # Read in the subject IDs for this study
  SubjStudID <- read.table(file = paste(LocSubjID, '/subjects_STUD_', k,'.txt', sep = ''),
                header = FALSE, sep = '\n', stringsAsFactors = FALSE) %>%
                  as.matrix(., ncol = 1)

  # Let's change analysis options in the base design file
    # Give new name
    Modified_Design <- Base_Design
    # Output directory: FSL will automatically create a new folder using the provided name + .gfeat
    output <- paste(Base_Output, simulID,'/', 'NSTUD_', k, sep = '')
    Modified_Design[grepl('outputdir', Modified_Design)] <- paste('set fmri(outputdir) "', output, '"', sep = '')

    # TR
    Modified_Design[grepl('fmri(tr)', Modified_Design, fixed = TRUE)] <-
        paste('set fmri(tr) ', TR, sep = '')

    # number of volumes and first-level analyses
    Modified_Design[grepl('npts', Modified_Design, fixed = TRUE)] <-
        paste('set fmri(npts) ', npts, sep ='')
    Modified_Design[grepl('multiple', Modified_Design, fixed = TRUE)] <-
        paste('set fmri(multiple) ', multiple, sep ='')

    # Directories of input data
      # Start with getting first and last line
      BeginInput <- which(grepl('4D AVW data or FEAT directory (1)', Modified_Design, fixed = TRUE)) - 1
      EndInput <- which(grepl('4D AVW data or FEAT directory (5)', Modified_Design, fixed = TRUE)) + 3

      # Now create new lines, based on number of subjects and their IDs
        # In HPC: create new folder NSTUD_$k_data, which will be next to NSTUD_$k.gfeat folder 
      InputLines <- c()
      for(s in 1:NSUB){
        InputLines <- rbind(InputLines, rbind(
          paste('# 4D AVW data or FEAT directory (',s,') ', sep = ''),
          if(MACHINE=='HPC') paste('set feat_files(',s,') "', Base_Input, k, '_data/', SubjStudID[s], '/results.feat"', sep = ''),
          if(MACHINE=='MAC') paste('set feat_files(',s,') "', Base_Input, SubjStudID[s],'/results.feat"', sep = ''),
          c(""))
        )
      }

      # Add this between first and last line of Modified_Design
      Modified_Design <- rbind(matrix(Modified_Design[c(1:BeginInput)], ncol = 1),
          InputLines,
          matrix(Modified_Design[c(EndInput:length(Modified_Design))], ncol = 1)
          ) %>% as.character()

    # Higher-level EV value for EV
      # Start and end of EV value
      BeginEV <- which(grepl('Higher-level EV value for EV 1 and input 1', Modified_Design, fixed = TRUE)) - 1
      EndEV <- which(grepl('Higher-level EV value for EV 1 and input 5', Modified_Design, fixed = TRUE)) + 3

      # Create lines
      EVLines <- c()
      for(s in 1:NSUB){
        EVLines <- rbind(EVLines, rbind(
          paste('# Higher-level EV value for EV 1 and input ',s, sep = ''),
          paste('set fmri(evg',s,'.1) 1', sep = ''),
          c(""))
        )
      }

      # Add this between first and last line of Modified_Design
      Modified_Design <- rbind(matrix(Modified_Design[c(1:BeginEV)], ncol = 1),
          EVLines,
          matrix(Modified_Design[c(EndEV:length(Modified_Design))], ncol = 1)
          ) %>% as.character()

    # Group membership
      # Start and end of group membership value
      BeginGM <- which(grepl('Group membership for input 1', Modified_Design, fixed = TRUE)) - 1
      EndGM <- which(grepl('Group membership for input 5', Modified_Design, fixed = TRUE)) + 3

      # Create lines
      GMLines <- c()
      for(s in 1:NSUB){
        GMLines <- rbind(GMLines, rbind(
          paste('# Group membership for input ',s, sep = ''),
          paste('set fmri(groupmem.',s,') 1', sep = ''),
          c(""))
        )
      }

      # Add this between first and last line of Modified_Design
      Modified_Design <- rbind(matrix(Modified_Design[c(1:BeginGM)], ncol = 1),
          GMLines,
          matrix(Modified_Design[c(EndGM:length(Modified_Design))], ncol = 1)
          ) %>% as.character()

  # Write design file to output directory
  writeLines(text = Modified_Design, con = paste(output, '.fsf', sep = ''))
}







