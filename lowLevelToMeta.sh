#!/bin/sh
#
#
#PBS -N Simul
#PBS -o output/output.file
#PBS -e error/error.file
#PBS -m a
#PBS -l walltime=10:00:00
#PBS -l vmem=30GB
#


#----------------------------------------------------#
# MODULES TO LOAD IN
module load R/3.2.1-intel-2015a
#----------------------------------------------------#



#----------------------------------------------------#
# CHANGE YOUR VSC NUMBER HERE AND GOD WILL DO THE REST
vsc=40728
#----------------------------------------------------#



#----------------------------------------------------#
# LOCATION OF SCRIPT TO RUN
srcdir=/user/scratch/gent/gvo000/gvo00022/vsc"$vsc"/Simulation
cd $srcdir
#----------------------------------------------------#


#----------------------------------------------------#
# GO TIME: CAREFUL WITH ORDER OF ARGUMENTS!
Rscript lowLevelToMeta.R ${PBS_ARRAYID}
#----------------------------------------------------#

echo "job finished"

