####################
#### TITLE:     Process the files from the SimpDesGrid.R simulation.
#### Contents:
####
#### Source Files: HPC - Version
#### First Modified: 24/02/2016
#### Notes:
#################



##
###############
### Notes
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
wd <- "/Users/hanbossier/Dropbox/PhD/PhDWork/Meta Analysis/R Code/Studie_Simulation/SimulationGit"
setwd(wd)

# Directories of the data for different takes
DATAwd <- list(
  'Take[8]' = "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/Take8"
	)

# Prefixes
prefix <- list(
  'Take[8]' = 'SDG_'
)

NUMDATAwd <- length(DATAwd)
currentWD <- 1

# Number of scenarios
NumScen.tmp <- matrix(c(
                8,30
              ), ncol=2, byrow=TRUE)
NumScen <- NumScen.tmp[currentWD,2]


# Number of conficence intervals
CIs <- c('norm','t','weightVar')
NumCI <- length(CIs)


# Number of executed simulations
nsim.tmp <- matrix(c(
                8,3000
              ), ncol=2, byrow=TRUE)
nsim <- nsim.tmp[currentWD,2]


# Dimension of brain
DIM.tmp <- array(NA, dim=c(NUMDATAwd,3))
	DIM.tmp[c(1),] <- c(16,16,16)
DIM <- DIM.tmp[currentWD,]

# Number of subjects and studies
TablesOverview <- list(
	'[1]' = '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/Take6/OverView.txt'
	)
overview.tmp <- matrix(c(			# This takes the element from TablesOverview
                8,1
              ), ncol=2, byrow=TRUE)
OverView.Sel <- read.table(file=TablesOverview[[overview.tmp[currentWD,2]]], header=TRUE)


# Load in libraries
library(AnalyzeFMRI)
library(fmri)
library(lattice)
library(gridExtra)
library(oro.nifti)
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
source('~/Dropbox/PhD/PhDWork/Meta\ Analysis/R\ Code/Studie_FixRan/FixRanStudyGit.git/Development/functions.R')


# Function for data wrangling: indicator for CI and true value
indicating <- function(UPPER, LOWER, trueVal){
	IND <- trueVal >= LOWER & trueVal <= UPPER
	COVERAGE <- apply(IND, 1, mean, na.rm=TRUE)
	return(COVERAGE)
}

# True value
trueVal <- 0

##
###############
### Data Wrangling
###############
##
# These are the objects that will be loaded in and saved.
objects <- c('CI.upper.norm','CI.lower.norm',
            'CI.upper.t','CI.lower.t',
            'CI.upper.weightVar','CI.lower.weightVar')

# For computational efficiency, we add all values that we load in into one long vector.
# Afterwards, we split them up again so that each column corresponds to a simulation.
# To do so, we need to know what the number of rows is that we will load in.
TEMPOBJ <- load(paste(DATAwd[[currentWD]],'/1/S1_',prefix[[currentWD]],'AllObjects_K1',sep=''))
NumberRows <- length(matrix(unlist(get(TEMPOBJ)[objects]),ncol=1))
  rm(TEMPOBJ)

# Load in R objects and combine them.
NumErr <- c()
OBJ.SIM <- matrix(NA,nrow=NumberRows,ncol=nsim)
OBJ.SCEN <- list()
print("Let's go!")
for(s in 1:NumScen){
  print(s)
  for(i in 1:nsim){
		# Load in objects
    OBJ.TMP <- try(load(paste(DATAwd[[currentWD]],'/',i,'/S',s,'_',prefix[[currentWD]],'AllObjects_K',i,sep='')),silent=TRUE)
    if(class(OBJ.TMP)=='try-error') {NumErr <- rbind(NumErr,c(s,i)); next}
    OBJ.SIM[,i] <- NeededData <- matrix(unlist(get(OBJ.TMP)[objects]),ncol=1)
    }
    OBJ.SCEN[[s]] <- matrix(OBJ.SIM,nrow=NumberRows)
    rm(OBJ.SIM)
    gc(verbose = FALSE)
    OBJ.SIM <- matrix(NA,nrow=NumberRows,ncol=nsim)
}
str(OBJ.SCEN)

# Save this object
save(OBJ.SCEN, file = paste(DATAwd[[currentWD]],'/Take8-CI.All', sep=''))


# Also save the missing simulations
if(is.null(NumErr)){
  cat('NULL', file = paste(DATAwd[[currentWD]], '/Take',currentWD, '-NumErr.txt', sep=''))
}else{
  colnames(NumErr) <- c('Scenario', 'Simulation')
  write.table(NumErr, file = paste(DATAwd[[currentWD]], '/Take8-NumErr.txt', sep=''),sep='\t',row.names=FALSE,quote=FALSE)
}


##
###############
### Compute CI coverages
###############
##
# Load in the data
load(paste(DATAwd[[currentWD]],'/Take8-CI.All', sep=''))

# CI coverages over all voxels and simulations
mean.coverage.norm <-
mean.coverage.t <-
mean.coverage.weightVar <-
		array(NA,dim=c(prod(DIM),NumScen))

OBJ.ID <- rep(objects,each=prod(DIM))

# For loop over all scenarios
for(s in 1:NumScen){
	print(s)
  OBJ.SCEN.ID <- OBJ.SCEN[[s]]
	mean.coverage.norm[,s] <- indicating(OBJ.SCEN.ID[which(OBJ.ID=='CI.upper.norm'),],OBJ.SCEN.ID[which(OBJ.ID=='CI.lower.norm'),],trueVal = trueVal)
	mean.coverage.t[,s] <- indicating(OBJ.SCEN.ID[which(OBJ.ID=='CI.upper.t'),],OBJ.SCEN.ID[which(OBJ.ID=='CI.lower.t'),],trueVal = trueVal)
	mean.coverage.weightVar[,s] <- indicating(OBJ.SCEN.ID[which(OBJ.ID=='CI.upper.weightVar'),],OBJ.SCEN.ID[which(OBJ.ID=='CI.lower.weightVar'),],trueVal = trueVal)
}


# Save CI coverage objects
save(mean.coverage.norm, file = paste(DATAwd[[currentWD]],'/Take8-mean.coverage.norm', sep=''))
save(mean.coverage.t, file = paste(DATAwd[[currentWD]],'/Take8-mean.coverage.t', sep=''))
save(mean.coverage.weightVar, file = paste(DATAwd[[currentWD]],'/Take8-mean.coverage.weightVar', sep=''))




