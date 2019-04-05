#!/bin/sh
#
#
#PBS -N HEPreProcessAct
#PBS -o output/output.file
#PBS -e error/error.file
#PBS -m a
#PBS -l walltime=05:00:00
#PBS -l vmem=20GB
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
# GO TIME: RUN THAT THING
Rscript PreProcessGLMvsIBMA_Act.R "HPC" "/user/data/gent/gvo000/gvo00022/vsc40728/IBMAvsGLM/Estimators/" "$SCEN" 
#----------------------------------------------------#



