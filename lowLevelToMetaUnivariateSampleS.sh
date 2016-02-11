#!/bin/sh

# Run this locally (MACHINE = MAC)
#----------------------------------------------------#
# LOCATION OF SCRIPT TO RUN
srcdir=/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/Take4
cd $srcdir
#----------------------------------------------------#

#----------------------------------------------------#
# NUMBER OF SCENARIOS
NSCEN=45
#----------------------------------------------------#


#----------------------------------------------------#
# NUMBER OF FOR LOOPS
PBS_ARRAYID=100
#----------------------------------------------------#


#----------------------------------------------------#
# CREATE THE FOLDERS IN RESULTS
for s in $(eval echo "{1..$PBS_ARRAYID}"); do
echo At simulation $s
cd Results
mkdir "$s"
	cd "$s"
  mkdir $(printf "SCEN_%1i " $(seq 1 $NSCEN))
cd $srcdir  
  # GO TIME: FOR LOOP OVER ALL SCENARIOS:
	for i in $(eval echo "{1..$NSCEN}"); do
	  Rscript lowLevelToMetaUnivariateSampleS.R "$s" "$i" "MAC" >> Text_Files/output_"$s".txt 2>> Text_Files/warning_"$s".txt
	done
done
#----------------------------------------------------#
