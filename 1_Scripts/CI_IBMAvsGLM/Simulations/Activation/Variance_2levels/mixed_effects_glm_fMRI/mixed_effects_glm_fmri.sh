#!/bin/sh
#
#
#PBS -N fMRIMixedEffGLM
#PBS -o output/output.file
#PBS -e error/error.file
#PBS -m a
#PBS -l walltime=10:00:00
#PBS -l vmem=30GB
#


#----------------------------------------------------#
# MODULES TO LOAD IN
module load R/3.5.1-foss-2018b 
#----------------------------------------------------#

#----------------------------------------------------#
# CHANGE YOUR VSC NUMBER HERE AND GOD WILL DO THE REST
vsc=40728
#----------------------------------------------------#

#----------------------------------------------------#
# LOCATION OF SCRIPT TO RUN
srcdir=/user/scratch/gent/gvo000/gvo00022/vsc"$vsc"/MixedEffGLMfMRI
cd $srcdir
#----------------------------------------------------#


#----------------------------------------------------#
# GO TIME!
Rscript mixed_effects_glm_fmri.R ${PBS_ARRAYID} "HPC"
#----------------------------------------------------#



