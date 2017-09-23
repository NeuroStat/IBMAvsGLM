#!/bin/bash

####################
#### TITLE:     Run second level analyses.
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
NSIM=1000
# Level
LEVEL=SecondLevel
# Number of subjects
NSUB=20
# Number of studies
NSTUD=5
# Scanning site of the 1000FC
SITE=Cambridge
# Smooting parameter in first level analyses
SMOOTHING=8mm
# Design in the first level analyses
DESIGN=event2
# Location of the scripts, data, and IDs for sampling subjects
if [ "$MACHINE" = MAC ] ; then
	SCRPT=/Users/hanbossier/Dropbox/PhD/PhDWork/Meta Analysis/R Code/Studie_Simulation/SimulationGit/1_Scripts/CI_IBMAvsGLM/1000FC/"${LEVEL}"
  BASE_DATA=/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IBMAvsGLM/Results
	LocSubjID="${BASE_DATA}"/"${SITE}"/"${LEVEL}"/IDs/"${SIMUL}"
fi
if [ "$MACHINE" = HPC ] ; then
	SCRPT=/user/scratch/gent/gvo000/gvo00022/vsc40728/1000FC/"${LEVEL}"
  BASE_DATA=/user/data/gent/gvo000/gvo00022/vsc40728
	LocSubjID="${SCRPT}"/Results/"${SITE}"/IDs/"${SIMUL}"
fi
# Maximum number of threads (if HPC, this is 16 I believe)
if [ "$MACHINE" = MAC ] ; then
	MaximumThreads=4
fi
if [ "$MACHINE" = HPC ] ; then
	MaximumThreads=16
fi

# Create directory of this simulation
if [ "$MACHINE" = MAC ] ; then
	cd "${BASE_DATA}"/"${SITE}"/"${LEVEL}"/"${SMOOTHING}"/"${DESIGN}"
	mkdir $SIMUL
	SIMULDIR="${BASE_DATA}"/"${SITE}"/"${LEVEL}"/"${SMOOTHING}"/"${DESIGN}"/"${SIMUL}"
fi
if [ "$MACHINE" = HPC ] ; then
	cd "${SCRPT}"/Results/"${SITE}"/"${SMOOTHING}"/"${DESIGN}"
	mkdir $SIMUL
	SIMULDIR="${SCRPT}"/Results/"${SITE}"/"${SMOOTHING}"/"${DESIGN}"/"${SIMUL}"
fi


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step Two: create second level design files
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Note: executed using R script: analyze_all_subjects.R
# Go to script locations
cd $SCRPT

# Run Rscript
Rscript "${SCRPT}"/analyze_all_studies.R $SIMUL $NSUB $NSTUD $SITE $MACHINE $SMOOTHING $DESIGN $LocSubjID $LEVEL &> output/R_design_$SIMUL.txt


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step Three: Run feat on all design files
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# We need to copy the subjects to a separate location, as running the simulations in parallel might interfere the analyses
# --> e.g. featregapply, which is the transformation to standard space happens in the subject specific folder.
# --> Part of the command unzips copes and varcopes.

# Go to directory of this simulation
cd "${SIMULDIR}"

# Put threads to 0
threads=0

# For loop over the amount of studies
for k in $(eval echo "{1..$NSTUD}")
do
	# Get into simulation folder
	cd "${SIMULDIR}"

	# Create folder NSTUD_$k_data.gfeat
	mkdir NSTUD_"$k"_data

	# Copy the subjects, read in from the LocSubjID file to NSTUD_$K.gfeat
	while IFS= read -r ID
	do
		cp -r "${BASE_DATA}"/"${SITE}"/"${SMOOTHING}"/"${DESIGN}"/"$ID" NSTUD_"$k"_data
	done < "$LocSubjID"/subjects_STUD_$k.txt

  # Run feat on the design files
  feat "${SIMULDIR}"/NSTUD_$k.fsf &> feat_design_nstud_$k.txt
  ((threads++))

  if [ $threads -eq "$MaximumThreads" ]; then
      wait
          threads=0
      fi

	# Delete the NSTUD_$k_data folder
	rm -r NSTUD_"$k"_data

	# Delete files in the NSTUD_$k.gfeat folder (takes too much space)
	cd NSTUD_$k.gfeat
		rm bg_image.nii.gz
		rm design.*
		rm design_*.*
		rm -r logs
		rm report*.html
		rm -r inputreg
		rm mean_func.nii.gz
		rm -r .files
		rm .ramp.gif

		cd cope1.feat
			rm -r	cluster_mask_zstat1.nii.gz
			rm -r cluster_zstat1_std.html
			rm -r cluster_zstat1_std.txt
			rm -r design.con
			rm -r	design.fsf
			rm -r design.grp
			rm -r	design.lcon
			rm -r	design.lev
			rm -r	design.mat
			rm -r	design.png
			rm -r	design.ppm
			rm -r	example_func.nii.gz
			rm -r	filtered_func_data.nii.gz
			rm -r	lmax_zstat1_std.txt
			rm -r	logs
			rm -r	mean_func.nii.gz
			rm -r	rendered_thresh_zstat1.nii.gz
			rm -r	rendered_thresh_zstat1.png
			rm -r	report_log.html
			rm -r	report_poststats.html
			rm -r	report_stats.html
			rm -r report.html
			rm -r thresh_zstat1.nii.gz
			rm -r thresh_zstat1.vol
			rm -r tsplot
			rm -r var_filtered_func_data.nii.gz

			cd stats
				rm -r logfile
				rm -r mean_random_effects_var1.nii.gz
				rm -r pe1.nii.gz
				rm -r res4d.nii.gz
				rm -r smoothness
				rm -r weights1.nii.gz
done










