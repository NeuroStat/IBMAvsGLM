####################
#### TITLE:     Simulate null data to calculate CI on ES and its coverage.
#### Contents:
####
#### Source Files: /Users/hanbossier/Dropbox/PhD/PhDWork/Meta Analysis/R Code/Studie_Simulation/SimulationGit/simulNull.R
#### First Modified: 12/01/2016
#### Notes:
#################



##
###############
### Notes
###############
##

# Simulating simple null images. These contain beta values. Then we calculate the simple OLS T-value.
# The images are group maps, containing only within study error (white noise only).
# No spatial smoothing.


##
###############
### Preparation
###############
##

# Reset working directory
rm(list=ls())
gc(verbose = FALSE)






