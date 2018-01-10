####################
#### TITLE:     Third Level of the resting state fMRI
#### Contents:
####
#### Source Files:
#### First Modified: 09/02/2017
#### Notes:
#################



##
###############
### Notes
###############
##

# We will compare the outcome of transforming the second level GLM to an ES and execute a meta-analsysis with a mixed effects approach.
# MEASURES:
#   * CI coverage
#   * Standardized bias
#   * Average CI length



##
###############
### Preparation
###############
##


# Reset working directory
rm(list=ls())
gc(verbose = FALSE)


# Take argument from master file
input <- commandArgs(TRUE)
  # ID of simulation
  simulID <- try(as.numeric(as.character(input)[1]),silent=TRUE)
  # Number of subjects
  NSUB <- try(as.numeric(as.character(input)[2]),silent=TRUE)
  # Number of studies
  NSTUD <- try(as.numeric(as.character(input)[3]),silent=TRUE)
  # Scanning site of the 1000FC
  SITE <- try(as.character(input)[4],silent=TRUE)
  # Which machine
  MACHINE <- try(as.character(input)[5],silent=TRUE)
  # Smoothing parameter
  SMOOTHING <- try(as.character(input)[6],silent=TRUE)
  # Design parameter
  DESIGN <- try(as.character(input)[7],silent=TRUE)
  # Level
  LEVEL <- try(as.character(input)[8],silent=TRUE)
    # If no machine is specified, then it has to be this machine in which we are testing code!
    if(is.na(MACHINE)){
      MACHINE <- 'MAC'
      simulID <- 1
      NSUB <- 20
      NSTUD <- 5
      SITE <- 'Cambridge'
      SMOOTHING <- '8mm'
      DESIGN <- 'boxcar10'
      LEVEL <- 'ThirdLevel'
    }

# Libraries
library(dplyr)
library(oro.nifti)

# Functions of study about CBMA
if(MACHINE=='HPC'){
  source('/user/scratch/gent/gvo000/gvo00022/vsc40728/1000FC/ThirdLevel/functions.R')
}
if(MACHINE=='MAC'){
  source('/Users/hanbossier/Dropbox/PhD/PhDWork/Meta Analysis/R Code/Studie_FixRan/FixRanStudyGit.git/Development/functions.R')
}




# Set base data directory of data processed at second level
if(MACHINE=='HPC'){
  data_dir <- paste('/user/scratch/gent/gvo000/gvo00022/vsc40728/1000FC', sep = '')
}
if(MACHINE=='MAC'){
  wd <- paste('/Users/hanbossier/Dropbox/PhD/PhDWork/Meta Analysis/R Code/Studie_Simulation/SimulationGit/1_Scripts/CI_IBMAvsGLM/1000FC/', LEVEL, sep = '')
  data_dir <- '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IBMAvsGLM/Results'
  setwd(wd)
}

# Give path to FSL
if(MACHINE=='HPC'){
  fslpath <- ''
}
if(MACHINE=='MAC'){
  fslpath <- '/usr/local/fsl/bin/'
}


# Dimensions of brain
DIM <- c(91,109,91)

##
###############
### Analysis parameters
###############
##

# Define the input and output directory
if(MACHINE=='HPC'){
  Base_Output <- paste(data_dir, '/', LEVEL, '/Results/', SITE,'/', SMOOTHING, '/', DESIGN, '/', sep = '')
  Base_Input <- paste(data_dir, '/SecondLevel/Results/', SITE, '/', SMOOTHING, '/', DESIGN, '/', sep = '')
}
if(MACHINE=='MAC'){
  Base_Output <- paste(data_dir, '/',SITE,'/', LEVEL, '/', SMOOTHING, '/', DESIGN, '/', sep = '')
  Base_Input <- paste(data_dir, '/',SITE,'/SecondLevel/', SMOOTHING, '/', DESIGN, '/', sep = '')
}




##
###############
### Read in second level data
###############
##

# Vector with amount of subjects in each study
NSUBstud <- rep(NSUB, NSTUD)

# Empty mask, T-map, COPE, VARCOPE and DOF vector
MASKstud <- TMAPstud <- DOFstud <- VARCOPEstud <- COPEstud <- array(NA, dim = c(prod(DIM),NSTUD))

# Read MASK, COPE, VARCOPE, T-map and dof of each study
print("Reading in second level studies")
for(k in 1:NSTUD){
  # Read in mask first
  MASKstud[,k] <- readNIfTI(paste(Base_Input,simulID,'/NSTUD_',k,'.gfeat/mask.nii.gz',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]

  COPEstud[,k] <- readNIfTI(paste(Base_Input,simulID,'/NSTUD_',k,'.gfeat/cope1.feat/stats/cope1.nii.gz',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,] %>%
    matrix(.,ncol = 1) %>% data.frame('cope' = .) %>% mutate(cope = replace(cope, MASKstud[,k] == 0, NA)) %>% as.matrix(.,ncol = 1)

  VARCOPEstud[,k] <- readNIfTI(paste(Base_Input,simulID,'/NSTUD_',k,'.gfeat/cope1.feat/stats/varcope1.nii.gz',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,] %>%
  matrix(.,ncol = 1) %>% data.frame('varcope' = .) %>% mutate(varcope = replace(varcope, MASKstud[,k] == 0, NA)) %>% as.matrix(.,ncol = 1)

  TMAPstud[,k] <- readNIfTI(paste(Base_Input,simulID,'/NSTUD_',k,'.gfeat/cope1.feat/stats/tstat1.nii.gz',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,] %>%
  matrix(.,ncol = 1) %>% data.frame('tstat' = .) %>% mutate(tstat = replace(tstat, MASKstud[,k] == 0, NA)) %>% as.matrix(.,ncol = 1)

  DOFstud[,k] <- readNIfTI(paste(Base_Input,simulID,'/NSTUD_',k,'.gfeat/cope1.feat/stats/tdof_t1.nii.gz',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,] %>%
  matrix(.,ncol = 1) %>% data.frame('tdof' = .) %>% mutate(tdof = replace(tdof, MASKstud[,k] == 0, NA)) %>% as.matrix(.,ncol = 1)

}
print("Preparing for third level")
# For MA approach: convert to Hedges' g and varHedge
# Transform to one vector data with second column the amount of subjects
TMAPstud <- data.frame('tvalue' = matrix(TMAPstud, ncol = 1), 'N' = rep(NSUBstud, each = prod(DIM)))

# Transform to an ES using hedgeG function, all studies in one vector
HedgeG <- hedgeG(t = TMAPstud$tvalue, N = TMAPstud$N)

# Calculate variance of ES
VarianceHedgeG <- varHedge(g = HedgeG, N = TMAPstud$N)

# Weights of each study
weigFix <- 1/VarianceHedgeG

# Now put in a vector with col = study
STHEDGE <- matrix(HedgeG, ncol = NSTUD)
STWEIGHTS <- matrix(weigFix, ncol = NSTUD)

# Clean up objects
rm(weigFix, VarianceHedgeG, HedgeG, TMAPstud)


########################################################################################################################################################################
########################################################################################################################################################################
print("THIRD LEVEL")
####************####
#### META-ANALYSIS: classical approach
####************####
# Calculate weighted average.
MA.WeightedAvg <- (apply((STHEDGE*STWEIGHTS),1,sum))/(apply(STWEIGHTS,1,sum))

# CI for weighted average based on weighted variance CI
CI.MA.weightedVariance <- (apply((STWEIGHTS*(STHEDGE - MA.WeightedAvg)^2),c(1),sum))/((NSTUD - 1) * apply(STWEIGHTS,1,sum))
CI.MA.upper.weightVar <- matrix(MA.WeightedAvg,ncol=1) + (qt(0.975,df=NSTUD-1) * sqrt(matrix(CI.MA.weightedVariance,ncol=1)))
CI.MA.lower.weightVar <- matrix(MA.WeightedAvg,ncol=1) - (qt(0.975,df=NSTUD-1) * sqrt(matrix(CI.MA.weightedVariance,ncol=1)))


########################################################################################################################################################################
########################################################################################################################################################################

####************####
#### IBMA: 3e level using FLAME
####************####

# Write auxiliarly files to Base_Output,simulID. We need:
  # STCOPE in nifti
  # STVARCOPE in nifti
  # 4D mask
  # design.mat file
  # design.grp file
  # design.con file

  #----- 1 ----#
  ### Design.mat
  fileCon <- paste(Base_Output,simulID,"/STdesign.mat",sep="")
  # Text to be written to the file
  cat('/NumWaves\t1
  /NumPoints\t',paste(NSTUD,sep=''),'
  /PPheights\t\t1.000000e+00

  /Matrix
  ',rep("1.000000e+00\n",NSTUD),file=fileCon)

  #----- 2 ----#
  ### Design.con
  fileCon <- file(paste(Base_Output,simulID,"/STdesign.con", sep=""))
  	writeLines('/ContrastName1	Group Average
  /NumWaves	1
  /NumContrasts	1
  /PPheights		1.000000e+00
  /RequiredEffect		5.034

  /Matrix
  1.000000e+00
  ',fileCon)
  close(fileCon)

  #----- 3 ----#
  ### Design.grp
  fileCon <- paste(Base_Output,simulID,"/STdesign.grp",sep="")
  # Text to be written to the file
  cat('/NumWaves\t1
  /NumPoints\t',paste(NSTUD,sep=''),'

  /Matrix
  ',rep("1\n",NSTUD),file=fileCon)

  #----- 4 ----#
  ### STCOPE.nii
  STCOPE4D <- nifti(img=array(COPEstud,dim=c(DIM,NSTUD)),dim=c(DIM,NSTUD),datatype = 16)
  writeNIfTI(STCOPE4D, filename = paste(Base_Output,simulID,'/STCOPE',sep=''),gzipped=FALSE)

  #----- 5 ----#
  ### VARCOPE.nii
  STVARCOPE4D <- nifti(img=array(VARCOPEstud,dim=c(DIM,NSTUD)),dim=c(DIM,NSTUD),datatype = 16)
  writeNIfTI(STVARCOPE4D, filename = paste(Base_Output,simulID,'/STVARCOPE',sep=''),gzipped=FALSE)

  #----- 6 ----#
  ### mask.nii
  ### This is intersection of all masks from previous subjects (flame crashes if voxels outside mask are included).
  MASKstud <- apply(MASKstud, 1, prod)
  mask <- nifti(img=array(MASKstud, dim=c(DIM,NSTUD)), dim=c(DIM,NSTUD), datatype=2)
  writeNIfTI(mask, filename = paste(Base_Output,simulID,'/mask',sep=''),gzipped=FALSE)

  #----- 7 ----#
  ### tdof_t1.nii
  TDOF4D <- nifti(img=array(DOFstud,dim=c(DIM,NSTUD)),dim=c(DIM,NSTUD),datatype = 16)
  writeNIfTI(TDOF4D, filename = paste(Base_Output,simulID,'/tdof_4D',sep=''),gzipped=FALSE)


# FSL TIME!
setwd(paste(Base_Output,simulID, sep = ''))
command <- paste(fslpath, 'flameo --cope=STCOPE --vc=STVARCOPE --dvc=tdof_4D --mask=mask --ld=GLM_stats --dm=STdesign.mat --cs=STdesign.grp --tc=STdesign.con --runmode=flame1', sep='')
Sys.setenv(FSLOUTPUTTYPE="NIFTI")
system(command)

### Now CI around COPE
GLM.COPE <- matrix(readNIfTI(paste(Base_Output,simulID,"/GLM_stats/cope1.nii",sep=""), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,],ncol=1)
GLM.SE <- sqrt(matrix(readNIfTI(paste(Base_Output,simulID,"/GLM_stats/varcope1.nii",sep=""), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,],ncol=1))
  # Degrees of freedom:
  tdof_t1 <- matrix(readNIfTI(paste(Base_Output,simulID,"/GLM_stats/tdof_t1.nii",sep=""), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,],ncol=1)

CI.GLM.upper.t <- GLM.COPE +  (qt(0.975,df=tdof_t1) * GLM.SE)
CI.GLM.lower.t <- GLM.COPE -  (qt(0.975,df=tdof_t1) * GLM.SE)


########################################################################################################################################################################
########################################################################################################################################################################


# Write objects to separate .txt files
write.table(CI.MA.upper.weightVar, file = paste(Base_Output,simulID, '/CI_MA_upper_weightVar.txt', sep = ''),
              row.names = FALSE, col.names = FALSE)
write.table(CI.MA.lower.weightVar, file = paste(Base_Output,simulID, '/CI_MA_lower_weightVar.txt', sep = ''),
              row.names = FALSE, col.names = FALSE)
write.table(MA.WeightedAvg, file = paste(Base_Output,simulID, '/MA_WeightedAvg.txt', sep = ''),
            row.names = FALSE, col.names = FALSE)
write.table(CI.GLM.upper.t, file = paste(Base_Output,simulID, '/CI_GLM_upper_t.txt', sep = ''),
            row.names = FALSE, col.names = FALSE)
write.table(CI.GLM.lower.t, file = paste(Base_Output,simulID, '/CI_GLM_lower_t.txt', sep = ''),
              row.names = FALSE, col.names = FALSE)
write.table(GLM.COPE, file = paste(Base_Output,simulID, '/GLM_COPE.txt', sep = ''),
              row.names = FALSE, col.names = FALSE)
write.table(CI.MA.weightedVariance, file = paste(Base_Output,simulID, '/CI_MA_weightedVariance.txt', sep = ''),
              row.names = FALSE, col.names = FALSE)
write.table(STHEDGE, file = paste(Base_Output,simulID, '/STHEDGE.txt', sep = ''),
              row.names = FALSE, col.names = FALSE)
write.table(STWEIGHTS, file = paste(Base_Output,simulID, '/STWEIGHTS.txt', sep = ''),
              row.names = FALSE, col.names = FALSE)

# Put objects in list
# ObjectsRestMAvsGLM <- list(
#   'CI.MA.upper.weightVar' = CI.MA.upper.weightVar,
#   'CI.MA.lower.weightVar' = CI.MA.lower.weightVar,
#   'MA.WeightedAvg' = MA.WeightedAvg,
#   'CI.GLM.upper.t' = CI.GLM.upper.t,
#   'CI.GLM.lower.t' = CI.GLM.lower.t,
#   'GLM.COPE' = GLM.COPE,
#   'CI.MA.weightedVariance' = CI.MA.weightedVariance,
#   'STHEDGE' = STHEDGE,
#   'STWEIGHTS' = STWEIGHTS
# )
#
# # Save list in output folder
# save(ObjectsRestMAvsGLM, file = paste(Base_Output,simulID,'/ObjectsRestMAvsGLM_',simulID,'.RData',sep=''))

# Clean up objects
rm(MA.WeightedAvg, CI.MA.weightedVariance, CI.MA.upper.weightVar, CI.MA.lower.weightVar)






