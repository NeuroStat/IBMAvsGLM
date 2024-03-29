#!/bin/sh
#
#
#PBS -N MixedEffGLM
#PBS -o output/output.file
#PBS -e error/error.file
#PBS -m a
#PBS -l walltime=03:00:00
#PBS -l vmem=30GB
#


#----------------------------------------------------#
# MODULES TO LOAD IN
module load R/3.2.3-intel-2016a
#----------------------------------------------------#

#----------------------------------------------------#
# CHANGE YOUR VSC NUMBER HERE AND GOD WILL DO THE REST
vsc=40728
#----------------------------------------------------#

#----------------------------------------------------#
# LOCATION OF SCRIPT TO RUN
srcdir=/user/scratch/gent/gvo000/gvo00022/vsc"$vsc"/MixedEffGLM
cd $srcdir
#----------------------------------------------------#


#----------------------------------------------------#
# GO TIME!
Rscript mixed_effects_glm.R ${PBS_ARRAYID} "HPC"
#----------------------------------------------------#



