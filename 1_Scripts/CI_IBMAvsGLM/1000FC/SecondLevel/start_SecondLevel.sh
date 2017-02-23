#!/bin/sh
#
#
#PBS -N IBMAvsGLM
#PBS -o output/output.file
#PBS -e error/error.file
#PBS -m a
#PBS -l walltime=5:00:00
#

#----------------------------------------------------#
# MODULES TO LOAD IN
module load R/3.2.3-intel-2016a
module load FSL/5.0.9-intel-2016a
module swap Java Java/1.8.0_92
. $FSLDIR/etc/fslconf/fsl.sh
#----------------------------------------------------#

# Note: newer version of FSL available...

#----------------------------------------------------#
# LOCATION OF SCRIPT TO RUN
srcdir=/user/scratch/gent/gvo000/gvo00022/vsc40728/1000FC/SecondLevel
cd $srcdir
#----------------------------------------------------#


#----------------------------------------------------#
# GO TIME!
./SecondLevel.sh ${PBS_ARRAYID}
#----------------------------------------------------#



echo "job finished"
