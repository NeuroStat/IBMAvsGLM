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
# Location of the scripts, data, and IDs for sampling subjects
if [ "$MACHINE" = MAC ] ; then
	SCRPT=/Users/hanbossier/Dropbox/PhD/PhDWork/Meta Analysis/R Code/Studie_Simulation/SimulationGit/1_Scripts/CI_IBMAvsGLM/1000FC/SecondLevel
  BASE_DATA=/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IBMAvsGLM/Results
	LocSubjID="${BASE_DATA}"/"${SITE}"/SecondLevel/IDs/"${SIMUL}"
fi
if [ "$MACHINE" = HPC ] ; then
	SCRPT=/user/scratch/gent/gvo000/gvo00022/vsc40728/1000FC/SecondLevel
  BASE_DATA=/user/data/gent/gvo000/gvo00022/vsc40728
	LocSubjID="${BASE_DATA}"/"${SITE}"/SecondLevel/IDs/"${SIMUL}"
fi
# Maximum number of threads (if HPC, this is 16 I believe)
if [ "$MACHINE" = MAC ] ; then
	MaximumThreads=4
fi
if [ "$MACHINE" = HPC ] ; then
	MaximumThreads=16
fi

# Create directory of this simulation
cd "${BASE_DATA}"/"${SITE}"/SecondLevel/"${SMOOTHING}"/"${DESIGN}"
mkdir $SIMUL

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step Two: create second level design files
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Note: executed using R script: analyze_all_subjects.R
# Go to script locations
cd $SCRPT

# Run Rscript
Rscript "${SCRPT}"/analyze_all_studies.R $SIMUL $NSUB $NSTUD $SITE $MACHINE $SMOOTHING $DESIGN $LocSubjID &> output/R_design_$SIMUL.txt


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step Three: Run feat on all design files
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# We need to copy the subjects to a separate location, as running the simulations in parallel might interfere the analyses
# --> e.g. featregapply, which is the transformation to standard space happens in the subject specific folder.
# --> Part of the command unzips copes and varcopes.

# Go to directory of this simulation
cd "${BASE_DATA}"/"${SITE}"/SecondLevel/"${SMOOTHING}"/"${DESIGN}"/"${SIMUL}"

# Put threads to 0
threads=0

# For loop over the amount of studies
for k in $(eval echo "{1..$NSTUD}")
do
	# Get into folder
	cd "${BASE_DATA}"/"${SITE}"/SecondLevel/"${SMOOTHING}"/"${DESIGN}"/"${SIMUL}"

	# Create folder NSTUD_$K.gfeat
	mkdir NSTUD_"$k"_data

	# Copy the subjects, read in from the LocSubjID file to NSTUD_$K.gfeat
	while IFS= read -r ID
	do
		cp -r "${BASE_DATA}"/"${SITE}"/"${SMOOTHING}"/"${DESIGN}"/"$ID" NSTUD_"$k"_data
	done < "$LocSubjID"/subjects_STUD_$k.txt

  # Run feat on the design files
  feat $BASE_DATA/$SITE/SecondLevel/$SMOOTHING/$DESIGN/$SIMUL/NSTUD_$k.fsf &> feat_design_nstud_$k.txt
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



# for i in /user/data/gent/gvo000/gvo00022/vsc40728/Cambridge/8mm/boxcar10/* ; do
# 	# Go to current directory
# 	cd $i
# 	# remove multiple masks
# 	cd results.feat/reg_standard
# 	rm mask.nii
# done
#


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step Four: 3e level: GLM (mixed effects)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#



# Now create the needed files for flameo (design.mat, design.con and design.grp and a mask) through R file designCrossVal.R
cd "${Studies}"
Rscript "${SCRPT}"/designCrossVal.R "$WRITE$MET/" "$(($NumSub + 1))" &> ROutput_designCrossVal.txt












