#!/bin/sh
#
#
#PBS -N GLMvsMA
#PBS -o output/output.file
#PBS -e error/error.file
#PBS -m a
#PBS -l walltime=11:20:00
#

#----------------------------------------------------#
# SWAP CLUSTERS
module swap cluster/delcatty
#----------------------------------------------------#

#----------------------------------------------------#
# MODULES TO LOAD IN
module load FSL/5.0.9-intel-2016a
module load R/3.2.3-intel-2016a
  . $FSLDIR/etc/fslconf/fsl.sh
#----------------------------------------------------#

#----------------------------------------------------#
# CHANGE YOUR VSC NUMBER HERE AND GOD WILL DO THE REST
vsc=40728
#----------------------------------------------------#

#----------------------------------------------------#
# LOCATION OF SCRIPT TO RUN
srcdir=/user/scratch/gent/gvo000/gvo00022/vsc"$vsc"/IBMAvsMA
cd $srcdir
#----------------------------------------------------#

#----------------------------------------------------#
# CHOOSE YOUR SCENARIO: ESTIMATOR FOR BETWEEN-STUDY
  # HETEROGENEITY IN CASE OF STNADARDIZED EFFECT SIZES
  # OPTIONS: GLM, DL, HE or REML
SCEN="HE"
#----------------------------------------------------#

#----------------------------------------------------#
# CHOOSE YOUR RATIO OF BETWEEN- OVER WITHIN-STUDY 
	# VARIABILITY
	# OPTIONS: 0.25, 0.5 or 0.75
RATIO=0.25
#----------------------------------------------------#

#----------------------------------------------------#
# CREATE THE FOLDERS IN RESULTS
cd Results
mkdir $SCEN
cd $SCEN
mkdir RatioBW_$RATIO
cd RatioBW_$RATIO
# This folder is used to write intermediate files!
mkdir ${PBS_ARRAYID}
cd $srcdir
#----------------------------------------------------#


#----------------------------------------------------#
# GO TIME: RUN THIS SCENARIO
Rscript MAvsGLM.R ${PBS_ARRAYID} "$SCEN" "$RATIO" "HPC" "$srcdir/Results/$SCEN/RatioBW_$RATIO/${PBS_ARRAYID}"
#----------------------------------------------------#
