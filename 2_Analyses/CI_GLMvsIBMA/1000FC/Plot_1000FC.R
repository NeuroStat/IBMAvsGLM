####################
#### TITLE:     Plotting results of 1000FC IBMA versus GLM resting state fMRI.
#### Contents:
####
#### Source Files:
#### First Modified: 15/02/2017
#### Notes:
#################



##
###############
### Notes
###############
##




##
###############
### Preparation
###############
##

# Reset working directory
rm(list=ls())
gc(verbose = FALSE)

# Date of today
date <- Sys.Date()

# Set starting seed
set.seed(11121990)

# Set WD
wd <- "/Users/hanbossier/Dropbox/PhD/PhDWork/Meta Analysis/R Code/Studie_Simulation/SimulationGit/2_Analyses/CI_GLMvsIBMA/1000FC"
setwd(wd)

# Directories of the data for different scenario's
DATAwd <- list(
  'Take[8mmBox10]' = "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IBMAvsGLM/Results/Cambridge/ThirdLevel/8mm/boxcar10"
	)
NUMDATAwd <- length(DATAwd)
currentWD <- 1

# Number of conficence intervals
CIs <- c('MA-weightVar','GLM-t')
NumCI <- length(CIs)

# Number of executed runs
nruns.tmp <- matrix(c(
                1,2500
              ), ncol=2, byrow=TRUE)
nruns <- nruns.tmp[currentWD,2]


# Number of subjects and studies
nsub <- 20
nstud <- 5

# Dimension of brain
DIM <- c(91,109,91)

# True value
trueVal <- 0

# Load in libraries
library(oro.nifti)
library(dplyr)

library(foreach)
library(doParallel)

library(AnalyzeFMRI)
library(fmri)
library(lattice)
library(gridExtra)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(Hmisc)
library(devtools)
library(neuRosim)
library(scatterplot3d)

# Function for data wrangling: indicator for CI and true value
indicating <- function(UPPER, LOWER, trueVal){
	IND <- trueVal >= LOWER & trueVal <= UPPER
	COVERAGE <- apply(IND, 1, mean, na.rm=TRUE)
	return(COVERAGE)
}

# Load in functions from FixRan study
#source('~/Dropbox/PhD/PhDWork/Meta\ Analysis/R\ Code/Studie_FixRan/FixRanStudyGit.git/Development/functions.R')


##
###############
### Data Wrangling
###############
##

######################################################
# First we create a universal mask over all iterations
######################################################

# Set warnings off
options(warn = -1)

# Vector to check progress
CheckProgr <- floor(seq(1,nruns,length.out=10))

# Do you want to make an universal mask again?
WRITEMASK <- FALSE
if(isTRUE(WRITEMASK)){
  # Vector with all masks in it
  AllMask <- c()

  # Load in the masks
  for(i in 1:nruns){
    # Print progress
    if(i %in% CheckProgr) print(paste('LOADING MASKS. NOW AT ', (i/nruns)*100, '%', sep = ''))

    # Try reading in mask, then go to one column and convert to data frame.
    CheckMask <- try(readNIfTI(paste(DATAwd[[currentWD]], '/', i,'/mask.nii', sep = ''))[,,,1] %>%
              matrix(.,ncol = 1) %>% data.frame(), silent = TRUE)
      # If there is no mask, skip iteration
      if(class(CheckMask) == "try-error") next

      # Some masks are broken: if all values are zero: REPORT
      if(all(CheckMask == 0)){print(paste("CHECK MASK AT ITERATION ", i, sep = "")); next}

    # Bind the masks of all iterations together
    AllMask <- bind_cols(AllMask, CheckMask)
    rm(CheckMask)
  }

  # Take product to have universal mask
  UnivMask <- apply(AllMask, 1, prod)

  # Better write this to folder
  niftiimage <- nifti(img=array(UnivMask, dim = DIM),dim=DIM)
  writeNIfTI(niftiimage,filename=paste(DATAwd[[currentWD]],'/universalMask',sep=''),gzipped=FALSE)
}
if(isTRUE(!WRITEMASK)){
  # Read in mask
  UnivMask <- readNIfTI(paste(DATAwd[[currentWD]],'/universalMask.nii', sep = ''))[,,] %>%
            matrix(.,ncol = 1)
}

###############################
#### Now load in all the data
###############################

# Load the naming structure of the data
load(paste(paste(DATAwd[[currentWD]], '/1/ObjectsRestMAvsGLM_1.RData',sep=''))); objects <- names(ObjectsRestMAvsGLM); rm(ObjectsRestMAvsGLM)
OBJ.ID <- c(rep(objects[!objects %in% c("STHEDGE","STWEIGHTS")], each=prod(DIM)), rep(c("STHEDGE","STWEIGHTS"), each=c(prod(DIM)*nstud)))

objects.CI <- objects[grepl(c('upper'), objects) | grepl(c('lower'), objects)]

# Due to computational constraints, we cannot keep all voxels of all parameters in one data frame.
# Hence we split up all parameters. For each object we have, we create an empty vector.
# Later on dimension will be: [MASKED VOXELS] x [NRUNS]
for(o in 1:length(objects)){
  assign(objects[o], c())
}

# Load in the data
for(i in 1:nruns){
  if(i %in% CheckProgr) print(paste('LOADING: AT ', (i/nruns)*100, '%', sep = ''))
  CheckObj <- try(load(paste(DATAwd[[currentWD]], '/', i, '/ObjectsRestMAvsGLM_', i, '.RData', sep = '')), silent = TRUE)
    if(class(CheckObj) == "try-error"){print(paste('Missing data in iteration ', i, sep = '')); next}

  # Loop over all parameters.
  for(o in 1:length(objects)){
    Obj <- objects[o]
    # Gather the masked voxels. Note: for STWEIGHTS and STHEDGE, we have values for each study.
    # This is why we stack all columns into one row before we bind them together.
    # Binding is for each run => NCOLS = NRUNS.
    values <- ObjectsRestMAvsGLM[[Obj]] %>% data.frame() %>% filter(., UnivMask == 1) %>%
            stack() %>% select(values) %>% bind_cols(get(objects[o]), .)
      assign(objects[o], values)
      rm(values,Obj)
      gc()
  }
}

prod(dim(CI.MA.upper.weightVar))
object.size(CI.MA.upper.weightVar)
print(object.size(CI.MA.upper.weightVar), units = "GB")


# ObjectsRestMAvsGLM[[Obj]] %>% data.frame() %>% filter(., maskIT == 1) %>%
#        stack() %>% select(values) %>% dim()
#   get(objects[9]) %>% dim()
#   dim(values)
#   values <-
#   ObjectsRestMAvsGLM[[Obj]] %>% data.frame() %>% filter(., maskIT == 1) %>% stack() %>% select(values) %>% head()
#   matrix(.,ncol = 1) %>% data.frame() %>% head()
#   bind_cols(get(objects[o]), .)
#   ?stack
#
# objects
#   str(ObjectsRestMAvsGLM)
#   ObjectsRestMAvsGLM$CI.MA.upper.weightVar
#     T1 <- is.na(ObjectsRestMAvsGLM$CI.MA.upper.weightVar)
#   ObjectsRestMAvsGLM$CI.GLM.upper.weightVar
#     T2 <- is.na(ObjectsRestMAvsGLM$CI.GLM.upper.t)
#     table(T1+T2)
# Split object into CI coverage emove NA values


# # Unlist data and bind to AllData
# AllData <- unlist(ObjectsRestMAvsGLM, use.names = FALSE) %>%
#   matrix(., ncol = 1) %>% data.frame() %>% bind_cols(AllData, .)



# # Remove first column
# AllData <- AllData[,-1]
#   head(AllData)
#
# # Does dimension of AllData match the lenght of the OBJ.ID?
# dim(AllData)[1]==length(OBJ.ID)

# Save this object
save(AllData, file = paste(DATAwd[[currentWD]], '/AllData.RData', sep = ''))





#########################################################
# Make one loop in which we calculate CI coverage, bias and length on each object.
#########################################################

# Function to count the number of instances in which true value is between lower and upper CI.
indicator <- function(UPPER, LOWER, trueval){
  IND <- trueval >= LOWER & trueval <= UPPER
  IND[is.na(IND)] <- 0
  return(IND)
}

# Funtion to count the number of recorded values
counting <- function(UPPER, LOWER){
  count <- (!is.na(UPPER) & !is.na(LOWER))
  return(count)
}

# Pre-define the CI coverage, length and bias vectors in which we sum the values
# After running nruns, divide by amount of obtained runs.
summed.coverage.IBMA <-
  summed.coverage.GLM <-
  array(0,dim=c(sum(UnivMask == 1),1))

# Keeping count of amount of values
counterMA <- counterGLM <- 0

registerDoParallel(cores = 3)
getDoParWorkers()

# Load in the data
t1 <- Sys.time()
foreach(i=1:50) %dopar% {
#for(i in 1:10){
  #if(i %in% CheckProgr) print(paste('PROCESSING. NOW AT ', (i/nruns)*100, '%', sep = ''))
  CheckObj <- try(load(paste(DATAwd[[currentWD]], '/', i, '/ObjectsRestMAvsGLM_', i, '.RData', sep = '')), silent = TRUE)
    if(class(CheckObj) == "try-error"){print(paste('Missing data in iteration ', i, sep = '')); next}

    # CI coverage: loop over the two procedures
    for(p in 1:2){
      objUP <- objects.CI[grepl(c('upper'), objects.CI)][p]
      objLOW <- objects.CI[grepl(c('lower'), objects.CI)][p]
      UP <- ObjectsRestMAvsGLM[[objUP]] %>% data.frame() %>% filter(., UnivMask == 1)
      LOW <- ObjectsRestMAvsGLM[[objLOW]] %>% data.frame() %>% filter(., UnivMask == 1)
      if(grepl('MA', x = objUP)){
        summed.coverage.IBMA[,1] <- summed.coverage.IBMA[,1] +
          indicator(UPPER = UP, LOWER = LOW, trueval = 0)
        counterMA <- counterMA + counting(UPPER = UP, LOWER = LOW)
      }else{
        summed.coverage.GLM[,1] <- summed.coverage.GLM[,1] +
          indicator(UPPER = UP, LOWER = LOW, trueval = 0)
        counterGLM <- counterGLM + counting(UPPER = UP, LOWER = LOW)
      }
      rm(objUP, objLOW, UP, LOW)
    }

    # # CI length: loop over the two procedures
    # for(p in 1:2){
    #   objUP <- objects.CI[grepl(c('upper'), objects.CI)][p]
    #   objLOW <- objects.CI[grepl(c('lower'), objects.CI)][p]
    #   UP <- ObjectsRestMAvsGLM[[objUP]] %>% data.frame() %>% filter(., UnivMask == 1)
    #   LOW <- ObjectsRestMAvsGLM[[objLOW]] %>% data.frame() %>% filter(., UnivMask == 1)
    #   if(grepl('MA', x = objUP)){
    #     summed.coverage.IBMA[,1] <- summed.coverage.IBMA[,1] +
    #       indicator(UPPER = UP, LOWER = LOW, trueval = 0)
    #   }else{
    #     summed.coverage.GLM[,1] <- summed.coverage.GLM[,1] +
    #       indicator(UPPER = UP, LOWER = LOW, trueval = 0)
    #   }
    #   rm(objUP, objLOW, UP, LOW)
    # }


    rm(CheckObj)
}
Sys.time() - t1


stopCluster()

mean(summed.coverage.IBMA/counterMA)
table(counterMA)
mean(summed.coverage.GLM/counterGLM)
table(counterGLM)

foreach(i=1:50, .combine='c') %dopar% {
  print(i)
}

# TRY MAKING FUNCTION INSIDE THE DOPAR LOOP



#########################################################
###################### CI COVERAGE ######################
#########################################################

# Calculate coverage
mean.coverage.weightVar.MA <-
  mean.coverage.t.GLM <-
  array(NA,dim=c(prod(DIM),1))
mean.coverage.weightVar.MA[,1] <- indicating(UPPER = AllData[which(OBJ.ID=='CI.MA.upper.weightVar'),],LOWER = AllData[which(OBJ.ID=='CI.MA.lower.weightVar'),],trueVal = trueVal)
mean.coverage.t.GLM[,1] <- indicating(UPPER = AllData[which(OBJ.ID=='CI.GLM.upper.t'),],LOWER = AllData[which(OBJ.ID=='CI.GLM.lower.t'),],trueVal = trueVal)


# Put the 2 coverages in a list
mean.coverages <- list('MA' = mean.coverage.weightVar.MA,'GLM' = mean.coverage.t.GLM)

# Mean over all voxels
CI.coverages <- data.frame(
	'Mean' = matrix(sapply(mean.coverages, FUN=function(...){apply(...,2,mean, na.rm = TRUE)}),ncol=1),
	'SD' = matrix(sapply(mean.coverages, FUN=function(...){apply(...,2,sd, na.rm = TRUE)}),ncol=1),
	'CI' = factor(CIs, levels=CIs, labels=CIs)
	)



#########################################################
####################### CI LENGTH #######################
#########################################################

# Calculate CI length
mean.length.weightVar.MA <-
mean.length.t.GLM <-
array(NA,dim=c(prod(DIM),1))

mean.length.weightVar.MA[,1] <- apply(AllData[which(OBJ.ID=='CI.MA.upper.weightVar'),] - AllData[which(OBJ.ID=='CI.MA.lower.weightVar'),],1,mean, na.rm = TRUE)
mean.length.t.GLM[,1] <- apply(AllData[which(OBJ.ID=='CI.GLM.upper.t'),] - AllData[which(OBJ.ID=='CI.GLM.lower.t'),],1,mean, na.rm = TRUE)


# Put the 2 lengths in a list
mean.lengths <- list('MA' = mean.length.weightVar.MA,'GLM' = mean.length.t.GLM)

# Average over all voxels
CI.lengths <- data.frame(
	'Mean' = matrix(sapply(mean.lengths, FUN=function(...){apply(...,2,mean,na.rm = TRUE)}),ncol=1),
	'SD' = matrix(sapply(mean.lengths, FUN=function(...){apply(...,2,sd,na.rm = TRUE)}),ncol=1),
	'CI' = factor(CIs, levels=CIs, labels=CIs)
	)



#########################################################
################### STANDARDIZED BIAS ###################
#########################################################

MA.SDBETA <- apply(AllData[which(OBJ.ID=='MA.WeightedAvg'),],1,sd,na.rm = TRUE)
MA.MEANBETA <- apply(AllData[which(OBJ.ID=='MA.WeightedAvg'),],1,mean,na.rm = TRUE)

GLM.SDBETA <- apply(AllData[which(OBJ.ID=='GLM.COPE'),],1,sd, na.rm = TRUE)
GLM.MEANBETA <- apply(AllData[which(OBJ.ID=='GLM.COPE'),],1,mean,na.rm = TRUE)

mean.bias.MA <- matrix(((abs(MA.MEANBETA)-trueVal)/(MA.SDBETA))*100,ncol=1)
mean.bias.GLM <- matrix(((abs(GLM.MEANBETA)-trueVal)/(GLM.SDBETA))*100,ncol=1)

# Put the 2 bias values in a list
mean.bias <- list('MA' = mean.bias.MA,'GLM' = mean.bias.GLM)

# Average over all voxels
CI.bias <- data.frame(
	'Mean' = matrix(sapply(mean.bias, FUN=function(...){apply(...,2,mean,na.rm = TRUE)}),ncol=1),
	'SD' = matrix(sapply(mean.bias, FUN=function(...){apply(...,2,sd,na.rm = TRUE)}),ncol=1),
	'CI' = factor(CIs, levels=CIs, labels=c('MA', 'GLM'))
	)
