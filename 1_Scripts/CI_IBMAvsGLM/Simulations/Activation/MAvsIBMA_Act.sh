#!/bin/sh
#
#
#PBS -N ActSimulMAvsIBMA
#PBS -o output/output.file
#PBS -e error/error.file
#PBS -m a
#PBS -l walltime=11:30:00
#PBS -l vmem=30GB
#


#----------------------------------------------------#
# MODULES TO LOAD IN
module load R/3.2.1-intel-2015a
module load FSL/5.0.6-ictce-5.5.0
  . $FSLDIR/etc/fslconf/fsl.sh
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
# NUMBER OF SCENARIOS
NSCEN=1
#----------------------------------------------------#

#----------------------------------------------------#
# CREATE THE FOLDERS IN RESULTS
cd Results
mkdir ${PBS_ARRAYID}
  cd ${PBS_ARRAYID}
  mkdir $(printf "SCEN_%1i " $(seq 1 $NSCEN))
cd $srcdir
#----------------------------------------------------#


#----------------------------------------------------#
for i in $(eval echo "{1..$NSCEN}"); do
  # GO TIME: FOR LOOP OVER ALL SCENARIOS:
  Rscript MAvsIBMA_Act.R ${PBS_ARRAYID} "$i" "HPC" "$srcdir/Results/${PBS_ARRAYID}/SCEN_$i"
done
#----------------------------------------------------#
