#!/bin/sh
#
#
#PBS -N ThirdLevelVar
#PBS -o output/output.file
#PBS -e error/error.file
#PBS -m a
#PBS -l walltime=03:00:00
#PBS -l vmem=30GB
#


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
srcdir=/user/scratch/gent/gvo000/gvo00022/vsc"$vsc"/Variance_3lvl
cd $srcdir
#----------------------------------------------------#

#----------------------------------------------------#
# Level of between-subject and between-study variability
BVAR=1
#----------------------------------------------------#

#----------------------------------------------------#
# GO TIME!
Rscript Var_3lvl.R ${PBS_ARRAYID} "HPC" "$srcdir/temp/${PBS_ARRAYID}" "$BVAR"
#----------------------------------------------------#



