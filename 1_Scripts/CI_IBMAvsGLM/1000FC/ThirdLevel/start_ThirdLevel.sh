#!/bin/sh
#
#
#PBS -N SimulMAvsGLM
#PBS -o output/output.file
#PBS -e error/error.file
#PBS -m a
#PBS -l walltime=11:30:00
#PBS -l vmem=10GB
#


#----------------------------------------------------#
# MODULES TO LOAD IN
module load R/3.2.3-intel-2016a
module load FSL/5.0.9-intel-2016a
  . $FSLDIR/etc/fslconf/fsl.sh
#----------------------------------------------------#

#----------------------------------------------------#
# CHANGE YOUR VSC NUMBER HERE AND GOD WILL DO THE REST
vsc=40728
#----------------------------------------------------#

#----------------------------------------------------#
# LOCATION OF SCRIPT TO RUN
srcdir=/user/scratch/gent/gvo000/gvo00022/vsc40728/1000FC/ThirdLevel
cd $srcdir
#----------------------------------------------------#


#----------------------------------------------------#
# SCRIPT PARAMETERS
# Which machine: HPC or MAC
MACHINE=HPC
# Level
LEVEL=ThirdLevel
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
#----------------------------------------------------#


#----------------------------------------------------#
# Make folder
mkdir "${srcdir}"/Results/"$SITE"/"$SMOOTHING"/"$DESIGN"/"${PBS_ARRAYID}"
#----------------------------------------------------#


#----------------------------------------------------#
# GO TIME!
Rscript ThirdLevel.R ${PBS_ARRAYID} $NSUB $NSTUD $SITE $MACHINE $SMOOTHING $DESIGN $LEVEL
#----------------------------------------------------#


echo "job finished"
