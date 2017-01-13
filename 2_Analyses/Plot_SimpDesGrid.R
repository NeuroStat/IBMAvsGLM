####################
#### TITLE:     Plotting results of the lowLevelToMetaSimpDesGrid.R data.
#### Contents:
####
#### Source Files:
#### First Modified: 23/02/2016
#### Notes:
#################



##
###############
### Notes
###############
##

# Blocked design for individual subjects.
# One condition.
# These N subjects are pooled using simple OLS pooling.
# The resulting images are converted to Hedges' g and pooled using fixed/random effects meta-analysis.



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


##
###############
### Data Wrangling
###############
##

# Data is being preprocessed in PreProcessSimpDesGrid.R!
  # Let's load the data in.
load(file=paste(DATAwd[[currentWD]],'/Take8-mean.coverage.norm', sep=''))
load(file = paste(DATAwd[[currentWD]],'/Take8-mean.coverage.t', sep=''))
load(file = paste(DATAwd[[currentWD]],'/Take8-mean.coverage.weightVar', sep=''))

# Put the 3 CI coverages in a list
mean.coverages <- list('norm' = mean.coverage.norm,'t' = mean.coverage.t, 'weightVar' = mean.coverage.weightVar)

CI.coverages <- data.frame(
	'Mean' = matrix(sapply(mean.coverages, FUN=function(...){apply(...,2,mean)}),ncol=1),
	'SD' = matrix(sapply(mean.coverages, FUN=function(...){apply(...,2,sd)}),ncol=1),
	'Scenario' = rep(seq(1,NumScen),NumCI),
	'CI' = rep(CIs, each = NumScen),
	'SampleSize' = rep(OverView.Sel[,'Subjects'],NumCI),
	'Studies' = rep(OverView.Sel[,'Studies'],NumCI)
	)

CI.coverages$CI <- factor(CI.coverages$CI, labels=CIs)
CI.coverages$SampleSize <- factor(CI.coverages$SampleSize)
CI.coverages$Studies <- factor(CI.coverages$Studies)


##
###############
### Plotting: coverage and weighted averages
###############
##

# Plotting the coverages
colours <- c('#1b9e77','#d95f02','#7570b3')
ggplot(CI.coverages, aes(x=SampleSize, y=Mean, colour=Studies)) +
  geom_point(aes(colour=Studies),size=1.3) +
  geom_line(aes(group=Studies), size=1) +
	geom_hline(yintercept=0.95,colour='red') +
	facet_wrap( ~ CI) +
  scale_x_discrete(name="Sample Size") +
  scale_colour_manual(values = colours, name='Amount of studies', labels = c(2,5,10)) +
  ggtitle(label='Coverages of 3 types of CI over all voxels and simulations') +
  theme(plot.title = element_text(lineheight=.4,size=13, face="plain"),
    axis.text.x = element_text(angle = 315, hjust = 0, vjust = 0.95),
		legend.key = element_rect(fill='#d9d9d9', colour = '#d9d9d9'),
		legend.background = element_rect(colour = '#d9d9d9', fill = '#d9d9d9'),
		strip.background = element_rect(fill='#d9d9d9'))






















