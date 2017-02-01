#!/bin/bash

####################
#### TITLE:     Run second and third level analyses.
#### Contents:
####
####
#### Source Files: ~\\Studie_Simulation/SimulationGit/1_Scripts/CI_IBMAvsGLM/1000FC/SecondLevel
#### First Modified: 01/02/2017
#### Notes:
#################




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step One: Global variables
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Simulation ID
SIMUL=$1
# Which machine: HPC or MAC
MACHINE=HPC
# Number of simulation
NSIM=3000
# Number of subjects
NSUB=20
# Number of studies
NSTUD=5
# Scanning site of the 1000FC
SITE=Cambridge
# Smooting parameter in first level analyses
SMOOTHING=8mm
# Design in the first level analyses
DESIGN=boxcar10
# Location of the scripts and data
if [ "$MACHINE" = MAC ] ; then
	SCRPT=/Users/hanbossier/Dropbox/PhD/PhDWork/Meta Analysis/R Code/Studie_Simulation/SimulationGit/1_Scripts/CI_IBMAvsGLM/1000FC/SecondLevel
  BASE_DATA=/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IBMAvsGLM/Results
fi
if [ "$MACHINE" = HPC ] ; then
	SCRPT=/user/scratch/gent/gvo000/gvo00022/vsc40728/1000FC/SecondLevel
  BASE_DATA=/user/data/gent/gvo000/gvo00022/vsc40728
fi
# Maximum number of threads (if HPC, this is 16 I believe)
if [ "$MACHINE" = MAC ] ; then
	MaximumThreads=4
fi
if [ "$MACHINE" = HPC ] ; then
	MaximumThreads=16
fi


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step Two: create second level design files
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Note: executed using R script: analyze_all_subjects.R
# Go to script locations
cd $SCRPT

# Run Rscript
Rscript "${SCRPT}"/analyze_all_studies.R $SIMUL $NSUB $NSTUD $SITE $MACHINE $SMOOTHING $DESIGN &> ROutput_design_files.txt


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step Three: Run feat on all design files
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Go to directory of this simulation
cd "${BASE_DATA}"/"${SITE}"/SecondLevel/"${SMOOTHING}"/"${DESIGN}"/"${SIMUL}"

# Put threads to 0
threads=0

# For loop over the amount of studies
for k in $(eval echo "{1..$NSTUD}")
do
  # Run feat on the design files
  feat $BASE_DATA/$SITE/SecondLevel/$SMOOTHING/$DESIGN/$SIMUL/NSTUD_$k.fsf &> feat_design_nstud_$k.txt
  ((threads++))

  if [ $threads -eq "$MaximumThreads" ]; then
      wait
          threads=0
      fi
done













