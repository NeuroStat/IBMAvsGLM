#!/bin/bash

####################
#### TITLE:     BLOCK DESIGN 1000FC -- REST
#### Contents:
####		Analyze all subjects of the 1000 functional connectomes project using
####		BLOCK design to create resting state fMRI.
####
#### Source Files: //Meta\ Analyis/R\ Code/Studie_Simulation/SimulationGit/
#### First Modified: 20/01/2017
#### Notes:
#################


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# NOTES
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Data downloaded from https://www.nitrc.org/frs/?group_id=296
# Renamed and saved on HDD.
# Scripts adapted from https://github.com/wanderine/ParametricMultisubjectfMRI


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step 0: Global variables
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Reset working memory
clear

# Source general settings for the analyses
source general_analysis_settings.sh

# Directory of design templates
design_directory=/Users/hanbossier/Dropbox/PhD/PhDWork/Meta\ Analysis/R\ Code/Studie_Simulation/SimulationGit/1_Scripts/CI_IBMAvsGLM/1000FC/FirstLevel

# Variable to select the dataset from 1000FC
Study=Cambridge
#Study=Beijing

# Which computer: HPC or MAC
COMPUTER=MAC
# Location of the data
if [ "$COMPUTER" = MAC ] ; then
	data=/Volumes/2_TB_WD_Elements_10B8_Han/PhD/FCON_1000/${Study}
	MaximumThreads=4 # Maximum number of CPU threads to use
fi
if [ "$COMPUTER" = HPC ] ; then
	# Which is your vsc number
	vsc=${13}
	data=/user/scratch/gent/gvo000/gvo00022/vsc"$vsc"/1000FC/${Study}
fi

# Which type of design: block or event
DESIGN=EVENT

# Type of design of 'task'
#DesignNew=boxcar10
#DesignNew=boxcar30
DesignNew=event2

# Length of boxcar periods (seconds)
BoxcarOffNew="set fmri(off1) 30"
BoxcarOnNew="set fmri(on1) 30"

# Highpass filter cutoff (twice boxcar length)
HighPassNew="set fmri(paradigm_hp) 60" #60 for 30 seconds on off, 20 for 10 seconds on off


#for Smoothing in 1 2 3 4 ; do
for Smoothing in 3 ; do

	threads=0

	# Amount of smoothing (FWHM in mm)
	if [ "$Smoothing" -eq "1" ]; then
		SmoothingNew="set fmri(smooth) 4.0"
		SmoothingOutputNew=4mm
	elif [ "$Smoothing" -eq "2" ]; then
		SmoothingNew="set fmri(smooth) 6.0"
		SmoothingOutputNew=6mm
	elif [ "$Smoothing" -eq "3" ]; then
		SmoothingNew="set fmri(smooth) 8.0"
		SmoothingOutputNew=8mm
	elif [ "$Smoothing" -eq "4" ]; then
		SmoothingNew="set fmri(smooth) 10.0"
		SmoothingOutputNew=10mm
	fi

	# Loop over all subjects
	for i in /Volumes/2_TB_WD_Elements_10B8_Han/PhD/FCON_1000/${Study}/* ; do
			# Check if fMRI data exists for this directory
		if [ -e ${i}/func/rest.nii.gz ]; then

			# Go to current directory
			cd $i

			# Get subject name
		 	SubjectNew=${PWD##*/}
		  echo "Processing" $SubjectNew

			# Go back to original directory containing design templates
			cd "${design_directory}"

			# Create data variable for each subject
			data_directory=${data}/${SubjectNew}/func

      #---------------
      # Copy template design
			if [ "$DESIGN" = "BLOCK" ]; then
      	cp "${design_directory}"/Design_templates/GLM${Study}.fsf ${data_directory}/
			elif [ "$DESIGN" = "EVENT" ]; then
				cp "${design_directory}"/Design_templates/GLM${Study}${DesignNew}.fsf ${data_directory}/
			fi

      # Change smoothing output
			if [ "$DESIGN" = "BLOCK" ]; then
      	sed -i '' "s/${SmoothingOutputOld}/${SmoothingOutputNew}/g" ${data_directory}/GLM${Study}.fsf
			elif [ "$DESIGN" = "EVENT" ]; then
				sed -i '' "s/${SmoothingOutputOld}/${SmoothingOutputNew}/g" ${data_directory}/GLM${Study}${DesignNew}.fsf
			fi

      # Change design output
			if [ "$DESIGN" = "BLOCK" ]; then
      	sed -i '' "s/${DesignOld}/${DesignNew}/g" ${data_directory}/GLM${Study}.fsf
			elif [ "$DESIGN" = "EVENT" ]; then
				sed -i '' "s/${DesignOld}/${DesignNew}/g" ${data_directory}/GLM${Study}${DesignNew}.fsf
			fi

      # Change subject name
			if [ "$DESIGN" = "BLOCK" ]; then
      	sed -i '' "s/${SubjectOld}/${SubjectNew}/g" ${data_directory}/GLM${Study}.fsf
			elif [ "$DESIGN" = "EVENT" ]; then
      	sed -i '' "s/${SubjectOld}/${SubjectNew}/g" ${data_directory}/GLM${Study}${DesignNew}.fsf
			fi

      # Change smoothing
			if [ "$DESIGN" = "BLOCK" ]; then
      	sed -i '' "s/${SmoothingOld}/${SmoothingNew}/g" ${data_directory}/GLM${Study}.fsf
			elif [ "$DESIGN" = "EVENT" ]; then
				sed -i '' "s/${SmoothingOld}/${SmoothingNew}/g" ${data_directory}/GLM${Study}${DesignNew}.fsf
			fi

      # Change boxcar period time
			if [ "$DESIGN" = "BLOCK" ]; then
      	sed -i '' "s/${BoxcarOffOld}/${BoxcarOffNew}/g" ${data_directory}/GLM${Study}.fsf
      	sed -i '' "s/${BoxcarOnOld}/${BoxcarOnNew}/g" ${data_directory}/GLM${Study}.fsf
      # Change highpass filter
      	sed -i '' "s/${HighPassOld}/${HighPassNew}/g" ${data_directory}/GLM${Study}.fsf
			fi

			#---------------
			# Run analyses in parallel
			# ------------ BLOCK ---------------
			if [ "$DESIGN" = "BLOCK" ]; then
				feat ${data_directory}/GLM${Study}.fsf &
				((threads++))

				if [ $threads -eq "$MaximumThreads" ]; then
						wait
								threads=0
						fi
			fi
			# ------------ EVENT ---------------
			if [ "$DESIGN" = "EVENT" ]; then
				feat ${data_directory}/GLM${Study}${DesignNew}.fsf &
				((threads++))

				if [ $threads -eq "$MaximumThreads" ]; then
						wait
								threads=0
						fi
			fi

			else
			echo "This directory does not contain any fMRI data"
		fi

	done
done

