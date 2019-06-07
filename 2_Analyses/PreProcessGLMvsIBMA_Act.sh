#!/bin/sh
#
#
#PBS -N DLPreProcess
#PBS -o output/
#PBS -e error/
#PBS -m a
#PBS -l walltime=05:00:00
#PBS -l vmem=20GB
#

#----------------------------------------------------#
# SWAP CLUSTERS!
# module swap cluster/delcatty
module swap cluster/golett
#----------------------------------------------------#

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
SCEN="DL"
#----------------------------------------------------#

#----------------------------------------------------#
# GO TIME: RUN THAT THING
Rscript PreProcessGLMvsIBMA_Act.R "HPC" "/user/scratch/gent/gvo000/gvo00022/vsc40728/IBMAvsMA/Results/" "$SCEN" 
#----------------------------------------------------#



