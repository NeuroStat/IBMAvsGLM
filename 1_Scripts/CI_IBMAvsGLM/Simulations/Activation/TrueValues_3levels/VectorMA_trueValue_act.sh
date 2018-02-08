#!/bin/sh
#
#
#PBS -N VectorMA
#PBS -o output/output.file
#PBS -e error/error.file
#PBS -m a
#PBS -l walltime=00:45:00
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
srcdir=/user/scratch/gent/gvo000/gvo00022/vsc"$vsc"/VectorMA
cd $srcdir
#----------------------------------------------------#


#----------------------------------------------------#
# GO TIME: let us do it in chunks of 10 simulations: let ID go from 1 --> 100
Rscript VectorMA_trueValue_act.R ${PBS_ARRAYID} "HPC" 
#----------------------------------------------------#
