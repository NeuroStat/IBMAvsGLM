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
  'Take[8mmBox10]' = "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IBMAvsGLM/Results/Cambridge/ThirdLevel/8mm/boxcar10",
  'Take[8mmEvent2]' = "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IBMAvsGLM/Results/Cambridge/ThirdLevel/8mm/event2"
	)
NUMDATAwd <- length(DATAwd)
currentWD <- 2

# Number of conficence intervals
CIs <- c('MA-weightVar','GLM-t')
NumCI <- length(CIs)

# Number of executed runs
nruns.tmp <- matrix(c(
                1,2500,
                2,500
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
library(lattice)
library(grDevices)
library(ggplot2)
library(data.table)

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

##
###############
### Data Wrangling
###############
##

######################################################
# First we create a universal mask over all iterations
######################################################

# Set warnings off
# options(warn = -1)

# Vector to check progress
CheckProgr <- floor(seq(1,nruns,length.out=10))

# Vector of simulations where we have a missing mask
missingMask <- c()

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
      if(class(CheckMask) == "try-error"){ missingMask <- c(missingMask, i); next}

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



#########################################################
# Make one loop in which we calculate CI coverage, bias and length on each object.
#########################################################

# Load the naming structure of the data
load(paste(paste(DATAwd[['Take[8mmBox10]']], '/1/ObjectsRestMAvsGLM_1.RData',sep=''))); objects <- names(ObjectsRestMAvsGLM); rm(ObjectsRestMAvsGLM)
OBJ.ID <- c(rep(objects[!objects %in% c("STHEDGE","STWEIGHTS")], each=prod(DIM)), rep(c("STHEDGE","STWEIGHTS"), each=c(prod(DIM)*nstud)))

objects.CI <- objects[grepl(c('upper'), objects) | grepl(c('lower'), objects)]


# Pre-define the CI coverage and length vectors in which we sum the values
# After running nruns, divide by amount of obtained runs.
# For bias, we work with VAR(X) = E(X**2) - E(X)**2 and a vector in which we sum the bias.
# Hence, we need to sum X**2 and X in a separate vector.
summed.coverage.IBMA <- summed.coverage.GLM <-
summed.length.IBMA <- summed.length.GLM <-
summed.X.IBMA <- summed.X.GLM <-
summed.X2.IBMA <- summed.X2.GLM <-
    array(0,dim=c(sum(UnivMask == 1),1))

# Keeping count of amount of values
counterMA <- counterGLM <- 0

# Load in the data
t1 <- Sys.time()
for(i in 1:20){
  if(i %in% CheckProgr) print(paste('PROCESSING. NOW AT ', (i/nruns)*100, '%', sep = ''))

    # CI coverage: loop over the two procedures
    for(p in 1:2){
      objUP <- objects.CI[grepl(c('upper'), objects.CI)][p] %>% gsub(".", "_",.,fixed = TRUE)
      objLOW <- objects.CI[grepl(c('lower'), objects.CI)][p] %>% gsub(".", "_",.,fixed = TRUE)

      UP <- try(fread(file = paste(DATAwd[[currentWD]], '/', i, '/', objUP, '.txt', sep = ''), header = FALSE) %>% filter(., UnivMask == 1), silent = TRUE)
              if(class(UP) == "try-error"){print(paste('Missing data in iteration ', i, sep = '')); next}
      LOW <- fread(file = paste(DATAwd[[currentWD]], '/',i, '/', objLOW, '.txt', sep = ''), header = FALSE) %>% filter(., UnivMask == 1)
      if(grepl('MA', x = objUP)){
        # CI coverage: add when true value in CI
        summed.coverage.IBMA[,1] <- summed.coverage.IBMA[,1] +
          indicator(UPPER = UP, LOWER = LOW, trueval = 0)
        # CI length: sum the length
        summed.length.IBMA[,1] <- summed.length.IBMA[,1] + as.matrix(UP - LOW)
        # Add one to the count (if data is available)
        counterMA <- counterMA + counting(UPPER = UP, LOWER = LOW)
      }else{
        # GLM procedure: CI coverage
        summed.coverage.GLM[,1] <- summed.coverage.GLM[,1] +
          indicator(UPPER = UP, LOWER = LOW, trueval = 0)
        # CI length: sum the length
        summed.length.GLM[,1] <- summed.length.GLM[,1] + as.matrix(UP - LOW)
        # Count
        counterGLM <- counterGLM + counting(UPPER = UP, LOWER = LOW)
      }
      rm(objUP, objLOW, UP, LOW)
    }

    # Standardized bias: read in weighted average / cope
    WAVG <- fread(file = paste(DATAwd[[currentWD]], '/', i, '/MA_WeightedAvg.txt', sep = ''), header = FALSE) %>% filter(., UnivMask == 1)
    GLMCOPE <- fread(file = paste(DATAwd[[currentWD]], '/', i, '/GLM_COPE', '.txt', sep = ''), header = FALSE) %>% filter(., UnivMask == 1)
    # Sum X
    summed.X.IBMA[,1] <- summed.X.IBMA[,1] + as.matrix(WAVG)
    summed.X.GLM[,1] <- summed.X.GLM[,1] + as.matrix(GLMCOPE)
    # Sum X**2
    summed.X2.IBMA[,1] <- summed.X2.IBMA[,1] + as.matrix(WAVG ** 2)
    summed.X2.GLM[,1] <- summed.X2.GLM[,1] + as.matrix(GLMCOPE ** 2)

}
Sys.time() - t1

# Calculate the average (over nsim) CI coverage, length and bias
Coverage.IBMA <- summed.coverage.IBMA/counterMA
Coverage.GLM <- summed.coverage.GLM/counterGLM

Length.IBMA <- summed.length.IBMA/counterMA
Length.GLM <- summed.length.GLM/counterGLM

# Formula: Var(X) = E(X**2) - [E(X)]**2
  # E(X**2) = sum(X**2) / n
  # E(X) = sum(X) / n
  # \hat{var(X)} = var(X) * (N / N-1)
  # \hat{SD} = sqrt(\hat{var(X)})
samplingSD.IBMA <- sqrt(((summed.X2.IBMA/(counterMA)) - ((summed.X.IBMA/counterMA)**2)) * (counterMA / (counterMA - 1)))
samplingSD.GLM <- sqrt(((summed.X2.GLM/(counterGLM)) - ((summed.X.GLM/counterGLM)**2)) * (counterGLM / (counterGLM - 1)))

# Standardized bias: true beta = 0
Bias.IBMA <- ((summed.X.IBMA / counterMA) - 0) / samplingSD.IBMA
Bias.GLM <- ((summed.X.GLM / counterGLM) - 0) / samplingSD.GLM


# Quick summary
mean(summed.coverage.IBMA/counterMA)
table(counterMA)
summary(summed.coverage.IBMA/counterMA)

mean(summed.coverage.GLM/counterGLM, na.rm = TRUE)
table(counterGLM)
# Histogram
hist(summed.coverage.IBMA/counterMA)
hist(summed.coverage.GLM/counterGLM)

# Violin plot
ViolinPlot <- data.frame(Coverages = c((summed.coverage.IBMA/counterMA),
                                    (summed.coverage.GLM/counterGLM)),
                        Method = rep(c('IBMA', 'GLM'), each = length(summed.coverage.IBMA)))
ViolinAll <- ggplot(ViolinPlot, aes(x = Method, y=Coverages)) + geom_violin() + coord_flip() +
  ggtitle("Plotting all voxels")

# Heatmap of the coverages
emptBrainIBMA <- emptBrainGLM <- array(NA, dim = prod(DIM))
emptBrainIBMA[UnivMask == 1] <- c(summed.coverage.IBMA/counterMA)
emptBrainGLM[UnivMask == 1] <- c(summed.coverage.GLM/counterGLM)

LevelPlotMACoV <- levelplot(array(emptBrainIBMA, dim = DIM)[,,c(36:46)], col.regions = topo.colors,
                	xlim=c(0,DIM[1]),ylim=c(0,DIM[2]), xlab = 'x', ylab = 'y',
                  main = 'CI coverage meta-analysis')
LevelPlotGLMCoV <- levelplot(array(emptBrainGLM, dim = DIM)[,,c(36:46)], col.regions = topo.colors,
                  xlim=c(0,DIM[1]),ylim=c(0,DIM[2]), xlab = 'x', ylab = 'y',
                  main = 'CI coverage GLM')

# Bias
emptBrainIBMA <- emptBrainGLM <- array(NA, dim = prod(DIM))
emptBrainIBMA[UnivMask == 1] <- Bias.IBMA
emptBrainGLM[UnivMask == 1] <- Bias.GLM
LevelPlotMABias <- levelplot(array(emptBrainIBMA, dim = DIM)[,,c(36:46)], col.regions = topo.colors,
                	xlim=c(0,DIM[1]),ylim=c(0,DIM[2]), xlab = 'x', ylab = 'y',
                  main = 'Standardized bias Meta-Analysis')
LevelPlotGLMBias <- levelplot(array(emptBrainGLM, dim = DIM)[,,c(36:46)], col.regions = topo.colors,
                  xlim=c(0,DIM[1]),ylim=c(0,DIM[2]), xlab = 'x', ylab = 'y',
                  main = 'Standardized bias GLM')
DifferenceBias <- levelplot(array(emptBrainIBMA - emptBrainGLM, dim = DIM)[,,c(36:46)], col.regions = topo.colors,
                	xlim=c(0,DIM[1]),ylim=c(0,DIM[2]), xlab = 'x', ylab = 'y',
                  main = 'Bias MA - GLM')

# CI length
emptBrainIBMA <- emptBrainGLM <- array(NA, dim = prod(DIM))
emptBrainIBMA[UnivMask == 1] <- Length.IBMA
emptBrainGLM[UnivMask == 1] <- Length.GLM
LevelPlotMACL <- levelplot(array(emptBrainIBMA, dim = DIM)[,,c(36:46)], col.regions = topo.colors,
                	xlim=c(0,DIM[1]),ylim=c(0,DIM[2]), xlab = 'x', ylab = 'y',
                  main = 'CI length Meta-Analysis')
LevelPlotGLMCL <- levelplot(array(emptBrainGLM, dim = DIM)[,,c(36:46)], col.regions = topo.colors,
                  xlim=c(0,DIM[1]),ylim=c(0,DIM[2]), xlab = 'x', ylab = 'y',
                  main = 'CI length GLM')



##
###############
### One voxel in 8mm, boxcar 10
###############
##

# Accordingly to the paradigm you fit on the resting state, we noticed that some
# voxels have low coverage accross all simulations.
# For instance in 8mm smoothing, boxcar = 10 SEC, look at row 48554 when you
# only take the masked voxels in one vector.
summedVoxel <- 0
voxel48554 <- array(NA, dim = c(500,2))
for(i in 1:500){
  UP <- try(fread(file = paste(DATAwd[[1]], '/', i, '/CI_MA_upper_weightVar.txt', sep = ''), header = FALSE) %>% filter(., UnivMask == 1), silent = TRUE)
        if(class(UP) == "try-error"){print(paste('Missing data in iteration ', i, sep = '')); next}
  LOW <- fread(file = paste(DATAwd[[1]], '/', i, '/CI_MA_lower_weightVar.txt', sep = ''), header = FALSE) %>% filter(., UnivMask == 1)
    voxel48554[i,1] <- summedVoxel + indicator(UPPER = UP, LOWER = LOW, trueval = 0)[48554]
    voxel48554[i,2] <- i
}

# Plot the cumulative sum of times it is included in the CI
data.frame(voxel48554) %>% filter(., !is.na(X1)) %>% mutate(., CumSum = cumsum(X1)) %>%
  rename('indicator' = X1, 'iteration' = X2, 'CumSum' = CumSum) %>%
  ggplot(., aes(x = iteration, y = CumSum)) + geom_line() +
  geom_abline(slope = 1, linetype = 2) +
  scale_x_continuous(limits = c(0,500)) +
  scale_y_continuous(limits = c(0,500))

# Where is this voxel located in 3D space?
rowID <- seq(1,length(UnivMask)) %>% data.frame(rowID = .)
voxID <- UnivMask %>% data.frame(UnivMask = .) %>% bind_cols(.,rowID) %>%
  filter(., UnivMask == 1) %>% slice(48554)
  # Now that we have the position of the voxel in the non-masked brain, convert to 3D
  emptBrain <- array(0, dim = prod(DIM))
  emptBrain[voxID$rowID] <- 1
  which(array(emptBrain, dim = DIM) == 1, arr.ind = TRUE)

# Levelplot
levelplot(array(emptBrain, dim = DIM))
# On top of a mask
MaskVoxel <- UnivMask
MaskVoxel[voxID$rowID] <- 2
levelplot(array(MaskVoxel, dim = DIM), col.regions = rainbow(n = 3),cuts = 2)
# zoom in
levelplot(array(MaskVoxel, dim = DIM)[,,30], col.regions = rainbow(n = 3),cuts = 2)

# Check values of this voxel in some simulations
i<-1
CheckVoxel <- try(load(paste(DATAwd[[currentWD]], '/', i, '/ObjectsRestMAvsGLM_', i, '.RData', sep = '')), silent = TRUE)
str(ObjectsRestMAvsGLM)

# Compare upper and lower limit of MA
ObjectsRestMAvsGLM[[1]][voxID$rowID]
ObjectsRestMAvsGLM[[1]][(voxID$rowID+1:10)]

ObjectsRestMAvsGLM[[2]][voxID$rowID]
ObjectsRestMAvsGLM[[2]][(voxID$rowID+1:10)]

# Look at weighted average
ObjectsRestMAvsGLM[[3]][voxID$rowID]
ObjectsRestMAvsGLM[[3]][(voxID$rowID+1:10)]

# Look at cope
ObjectsRestMAvsGLM[[6]][voxID$rowID]
ObjectsRestMAvsGLM[[6]][(voxID$rowID+1:10)]

# Look at upper and lower GLM CI
ObjectsRestMAvsGLM[[4]][voxID$rowID]
ObjectsRestMAvsGLM[[4]][(voxID$rowID+1:10)]

ObjectsRestMAvsGLM[[5]][voxID$rowID]
ObjectsRestMAvsGLM[[5]][(voxID$rowID+1:10)]



# Read in all simulations and gather info on this voxel, compared to 10 voxels furhter.
voxCOPE <- c()

for(i in sample(x = 1:500, size = 20)){
  if(i == 16) next
  if(i == 18) next
  if(i == 28) next
  tmp <- readNIfTI(paste(DATAwd[[currentWD]], '/', i, '/STCOPE.nii', sep = ''))[,,,1]
  voxCOPE <- c(voxCOPE, array(tmp, dim = prod(DIM))[voxID$rowID])
}

voxCOPE
t.test(voxCOPE, alternative = 'two.sided')




##
###############
### Voxels with above chance low/high coverage
###############
##

# Location of max COPE value at third level
emptBrainMax <- array(0, dim = prod(DIM))
emptBrainMin <- array(0, dim = prod(DIM))
for(i in 1:500){
  GCOPE <- fread(file = paste(DATAwd[[currentWD]], '/', i, '/GLM_COPE.txt', sep = ''), header = FALSE)
  posMAX <- which(GCOPE == max(GCOPE), arr.ind = TRUE)[1]
  posMIN <- which(GCOPE == min(GCOPE), arr.ind = TRUE)[1]

  emptBrainMax[posMAX] <- emptBrainMax[posMAX] + 1
  emptBrainMin[posMIN] <- emptBrainMin[posMIN] + 1
}
table(emptBrainMax)
table(emptBrainMin)

emptBrainMax[UnivMask == 0] <- emptBrainMin[UnivMask == 0] <- NA

zMAX <- which(array(emptBrainMax, dim = DIM) == max(array(emptBrainMax, dim = DIM), na.rm = TRUE),
arr.ind = TRUE)[3]
zMIN <- which(array(emptBrainMin, dim = DIM) == max(array(emptBrainMin, dim = DIM), na.rm = TRUE),
arr.ind = TRUE)[3]

levelplot(array(emptBrainMax, dim = DIM)[,,zMAX], main = 'Position out of 500 simulations of the maximum COPE value',
          col.regions = rainbow(max(emptBrainMax, na.rm = TRUE)), cuts = max(emptBrainMax, na.rm = TRUE))
levelplot(array(emptBrainMax, dim = DIM)[,,42], main = 'Position out of 500 simulations of the maximum COPE value',
          col.regions = rainbow(max(emptBrainMax, na.rm = TRUE)), cuts = max(emptBrainMax, na.rm = TRUE))
levelplot(array(emptBrainMin, dim = DIM)[,,zMIN], main = 'Position out of 500 simulations of the minimum COPE value',
          col.regions = rainbow(max(emptBrainMin, na.rm = TRUE)), cuts = max(emptBrainMin, na.rm = TRUE))


# Heatmaps of the COPE maps
for(i in 1:100){
  tmp <- readNIfTI(paste(DATAwd[[currentWD]], '/', i, '/STCOPE.nii', sep = ''))[,,,]
}

GCOPE_ALL <- data.frame('NA' = array(NA, dim = prod(DIM))) %>% tbl_df()
for(i in 1:500){
  GCOPE <- fread(file = paste(DATAwd[[currentWD]], '/', i, '/GLM_COPE.txt', sep = ''), header = FALSE, col.names = as.character(i))
  GCOPE_ALL <- bind_cols(GCOPE_ALL, GCOPE)
}
GCOPE_ALL <- select(GCOPE_ALL, -1)

meanGCOPE <- apply(GCOPE_ALL, 1, mean)
meanGCOPE[UnivMask==0] <- NA
levelplot(array(meanGCOPE, dim = DIM), main = 'Average third level COPE value over 500 iterations',
          col.regions = topo.colors)








#########################################################
# Parallel version
#########################################################

library(foreach)
library(doParallel)

registerDoParallel(cores = 3)
getDoParWorkers()


foreach(i=1:50, .combine='c') %dopar% {
  print(i)
}

# TRY MAKING FUNCTION INSIDE THE DOPAR LOOP



stopCluster()



###############################
#### Attempt to load in all data
###############################

# In the following piece of code, I tried to load in all data at once and bind them in one object.
# It turned out that this is not efficient/possible.

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

# Save this object
save(AllData, file = paste(DATAwd[[currentWD]], '/AllData.RData', sep = ''))




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
    # Write all data to text files
    write.table(ObjectsRestMAvsGLM[[Obj]], file = paste(DATAwd[[currentWD]], '/', i, '/',gsub('.', '_', x = Obj, fixed = TRUE),'.txt', sep = ''),
                row.names = FALSE, col.names = FALSE)
    rm(Obj)
  }
}

head(ObjectsRestMAvsGLM[[Obj]])
write.table(ObjectsRestMAvsGLM$CI.MA.upper.weightVar, file = 'CI_MA_upper_weightVar.txt', row.names = FALSE, col.names = FALSE)



###############################
#### Code that was used to caluclate CI coverage, length and bias
###############################



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




##
###############
### Older Code
###############
##

# CheckObj <- try(load(paste(DATAwd[[currentWD]], '/', i, '/ObjectsRestMAvsGLM_', i, '.RData', sep = '')), silent = TRUE)
#   if(class(CheckObj) == "try-error"){print(paste('Missing data in iteration ', i, sep = '')); next}

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



#
# MA.SDBETA <- apply(AllData[which(OBJ.ID=='MA.WeightedAvg'),],1,sd,na.rm = TRUE)
# MA.MEANBETA <- apply(AllData[which(OBJ.ID=='MA.WeightedAvg'),],1,mean,na.rm = TRUE)
#
# GLM.SDBETA <- apply(AllData[which(OBJ.ID=='GLM.COPE'),],1,sd, na.rm = TRUE)
# GLM.MEANBETA <- apply(AllData[which(OBJ.ID=='GLM.COPE'),],1,mean,na.rm = TRUE)
#
# mean.bias.MA <- matrix(((abs(MA.MEANBETA)-trueVal)/(MA.SDBETA))*100,ncol=1)
# mean.bias.GLM <- matrix(((abs(GLM.MEANBETA)-trueVal)/(GLM.SDBETA))*100,ncol=1)
#
#
# X <- c(1,4,2,5,6,3,2,7,8,5,4)
# var(X) * (10) / 11
# (sum(X**2) / length(X)) - ((sum(X)/length(X))**2)
# (mean(X**2) - (mean(X)**2)) * (11 / 10)
#
# var(X)
# ((sum(X**2) / 11) - ((sum(X)/11)**2)) * (11 / 10)
# sum((X - (sum(X)/11))**2)/10
#
# (sum(X**2)/10**2) - (X/10)
