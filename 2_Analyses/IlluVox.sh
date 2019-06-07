#!/bin/sh
#
#
#PBS -N FindVoxel
#PBS -o output/
#PBS -e error/
#PBS -m a
#PBS -l walltime=00:15:00
#

#----------------------------------------------------#
# SWAP CLUSTERS!
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
# GO TIME: RUN THAT THING
Rscript IlluVox.R "HPC" "/user/scratch/gent/gvo000/gvo00022/vsc40728/IBMAvsMA/Results" ${PBS_ARRAYID}
#----------------------------------------------------#



