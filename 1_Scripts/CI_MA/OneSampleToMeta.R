####################
#### TITLE:     Simulate null data for one sample test to meta-analysis.
#### Contents:
####
#### Source Files: Local - Version
#### First Modified: 04/02/2016
#### Notes:
#################



##
###############
### Notes
###############
##

# Here, we try to check coverages of a one sample test to a meta-analysis.


# Reset working directory
rm(list=ls())
gc(verbose = FALSE)


# Set starting seed.
starting.seed <- 1112
set.seed(starting.seed)

# Set WD
wd <- "/Users/hanbossier/Dropbox/PhD/PhDWork/Meta Analysis/R Code/Studie_Simulation/SimulationGit"
setwd(wd)


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



# Load in functions from FixRan study: THIS HAS TO COME AFTER ALL LIBRARIES ARE LOADED AS WE SOMETIMES FIX FUNCTIONS THAT ARE BUGGED IN THE PACKAGES
source_url('https://raw.githubusercontent.com/HBossier/FixRanStudyGit/master/Development/functions.R',sha1='c4c3b98288ab8a9bdf0d081f2ace902d5cd13e18')




##
###############
### Simulation steps
###############
##

nsim <- 50000
nsub <- 80
nstud <- 50
trueVal <- 0

# Number of conficence intervals
CIs <- c('norm','t','weightVar')
NumCI <- length(CIs)

# The CI's and the vector containing the indicators, which we will take the average of
ind.coverage.norm <- CI.upper.norm <- CI.lower.norm <- array(NA,dim=nsim)
ind.coverage.t <- CI.upper.t <- CI.lower.t <- array(NA,dim=nsim)
ind.coverage.weightVar <- CI.upper.weightVar <- CI.lower.weightVar <- array(NA,dim=nsim)

# The weighted mean values
wMean <- array(NA,dim=nsim)

# Indicator for progress in the for loop
INDICATOR <- seq(1,nsim,by=nsim/10)

# FOR LOOP
for(i in 1:nsim){
	if(i %in% INDICATOR)print(i)
	Stud <- array(NA,dim=nstud)
	varStud <- array(NA,dim=nstud)
	for(s in 1:nstud){
		# Subjects in group
		X1 <- trueVal + rnorm(n=nsub,mean=0,sd=2)
		S2_X1 <- var(X1)

		# Hedge g
		t <- (mean(X1))/sqrt(S2_X1/nsub)
		g <- hedgeG(t,nsub)
		Vg <- varHedge(g,nsub)

		# Save in vector
		Stud[s] <- g
		varStud[s] <- Vg

	}

	# Now meta-analysis
	# Weights
	W <- 1/varStud
	# Weighted means
	wMean[i] <- sum(W*Stud)/sum(W)

	# Variance of summary effect
	varSum <- 1/sum(W)

	# Confidence intervals
	CI.upper.norm[i] <- wMean[i] + (1.96 * sqrt(varSum))
	CI.lower.norm[i] <- wMean[i] - (1.96 * sqrt(varSum))
	CI.upper.t[i] <- wMean[i] + (qt(0.975,df=nsub-1)  * sqrt(varSum))
	CI.lower.t[i] <- wMean[i] - (qt(0.975,df=nsub-1)  * sqrt(varSum))
		# Weighted variance
		CI.weightedVariance <- sum(W*(Stud - wMean[i])**2)/((nstud - 1) * sum(W))
	CI.upper.weightVar[i] <- wMean[i] + (qt(0.975,df=nsub-1) * sqrt(CI.weightedVariance))
	CI.lower.weightVar[i] <- wMean[i] - (qt(0.975,df=nsub-1) * sqrt(CI.weightedVariance))

	# Coverage
	ind.coverage.norm[i] <- ifelse(trueVal >= CI.lower.norm[i] & trueVal <= CI.upper.norm[i], 1, 0)
	ind.coverage.t[i] <- ifelse(trueVal >= CI.lower.t[i] & trueVal <= CI.upper.t[i], 1, 0)
	ind.coverage.weightVar[i] <- ifelse(trueVal >= CI.lower.weightVar[i] & trueVal <= CI.upper.weightVar[i], 1, 0)

}

# Take the averages of the coverage indicators
mean.coverage.norm <- mean(ind.coverage.norm);mean.coverage.norm
mean.coverage.t <- mean(ind.coverage.t);mean.coverage.t
mean.coverage.weightVar <- mean(ind.coverage.weightVar);mean.coverage.weightVar

	# Indicators
	ind.coverages <- list('norm' = ind.coverage.norm, 't' = ind.coverage.t, 'weightVar' = ind.coverage.weightVar)




	##
	###############
	### Plotting the number of simulations needed
	###############
	##



	# Number of calculations
	BY <- 10
	Calcs <- seq(600,nsim,by=BY)
	  PreCalcs <- seq(501,(nsim-(BY-1)),by=BY)
	NumCalc <- length(Calcs)

	# MEANS's
	MEAN <- data.frame(
	    'Mean' = array(NA,dim=c(NumCI*NumCalc)),
	    'CI' = array(NA,dim=c(NumCI*NumCalc))
	    )

	# For loop
	for(j in 1:NumCI){
	 # Pre-LOOP
	  DataMEAN <- ind.coverages[[CIs[j]]]
	   dat.tmp <- DataMEAN[c(1:500)]
	   dat.mean <- dat.tmp
	    rm(dat.tmp)
	  for(i in 1:NumCalc){
	    Index <- ((j-1) * NumCalc) + i
	    dat.tmp <- DataMEAN[c(PreCalcs[i]:Calcs[i])]
	    dat.mean <- c(dat.mean,dat.tmp)

	    # Calculate MEAN
	    MEANsim.tmp <- mean(dat.mean,na.rm=TRUE)
	    MEAN[Index,'Mean'] <- MEANsim.tmp
	    MEAN[Index, 'CI'] <- j

	    # Remove objects
	    rm(MEANsim.tmp,dat.tmp)
	  }
	  rm(dat.mean)
	}

	# Data frame
	MEAN$Time <- rep(c(1:NumCalc),times=NumCI)
	  MEAN$CI <- factor(MEAN$CI, labels = CIs)
	head(MEAN)

	# Plot
	LabelTime <- seq(600,nsim,length.out=6)
	ggplot(MEAN, aes(x = Time, y = Mean, group = CI)) +
	  geom_line(aes(colour = CI)) +
		scale_x_continuous(name='Number of simulations', labels = LabelTime) + 
  scale_y_continuous(name='Average emperical coverage')





