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

# Date of today
date <- Sys.Date()

# Set starting seed
set.seed(11121990)

# Set WD
wd <- "/Users/hanbossier/Dropbox/PhD/PhDWork/Meta Analysis/R Code/Studie_Simulation/SimulationGit"
# Next are data for different scenarios
DATAwd <- list(
	'Smooth[1]' = "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/Smoothed",
	'NonSmooth[2]' = "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/NonSmoothed",
  'WhiteNonSmooth[3]' = "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/WhiteNonSmooth"
	)

prefix <- list(
  'Smooth[1]' = "",
  'NonSmooth[2]' = "NSm",
  'WhiteNonSmooth[3]' = "WNSm"
  )


NUMDATAwd <- length(DATAwd)
currentWD <- 3

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
  # Function to be fixed in neuRosim
  stimfunction<-function (totaltime, onsets, durations, accuracy)
{
  if (max(onsets) > totaltime) {
    stop("Mismatch between onsets and totaltime")
  }
  s <- rep(0, totaltime/accuracy)
  os <- onsets/accuracy
  dur <- durations/accuracy
  if (length(durations) == 1) {
    dur <- dur * rep(1, length(onsets))
  }
  else if (length(durations) != length(onsets)) {
    stop("Mismatch between number of onsets and number of durations.")
  }
  for (i in (1:length(onsets))) {
    if ((os[i] + dur[i]) <= totaltime/accuracy) {
      s[c(os[i]:(os[i] + dur[i]-1))] <- 1
    }
    else {
      s[c(os[i]:(totaltime/accuracy))] <- 1
    }
  }
  return(s)
}

# Load in functions from FixRan study
source('~/Dropbox/PhD/PhDWork/Meta\ Analysis/R\ Code/Studie_FixRan/FixRanStudyGit.git/Development/functions.R')


# Dimension of the 2D images
DIM <- c(64,64)


# True value
trueVal <- 0
# Within study variability
wstudSD <- 2

# Number of subjects, if used in simulation
nsub <- 20

##
###############
### Simulation steps
###############
##

########################################
## Start off with a normal distribution:
#   * 10.000 simulations
#   * Treat 64x64 image as 64*64 subjects from a normal distribution
#   * Construct CI around mean within simulation.
#   * Check coverage over all simulations

# Vector for mean coverage within simulation
nsim <- 10000
mean.coverage.norm <- array(NA,nsim)
DIM <- c(64,64)

# Start for loop
for(i in 1:nsim){
  # Beta image with true value
  TrueBetaImage <- array(trueVal,dim=c(DIM))
  # Add white noise
  BetaImage <- TrueBetaImage + rnorm(n=prod(DIM),trueVal,wstudSD)

  # Now calculate confidence interval
  CI.upper.norm <- mean(BetaImage) + (1.96 * (sd(BetaImage)/sqrt(prod(DIM))))
  CI.lower.norm <- mean(BetaImage) - (1.96 * (sd(BetaImage)/sqrt(prod(DIM))))

  # Put mean in vector
  mean.coverage.norm[i] <- ifelse(trueVal > CI.lower.norm & trueVal < CI.upper.norm, 1, 0)
}
# Get the average.
mean(mean.coverage.norm)





########################################
## Normal distribution of group maps based on OLS pooling of subjects:
#   * 1.000 simulations
#   * Create 64x64 beta-images for N subjects.
#   * Pool these with ordinary OLS pooling method: T-maps
#   * Construct CI in each voxel, around T-value (using t-distr).
#   * Check coverage in each voxel over all simulations


# Get vector for mean coverage within simulation
nsim <- 1000
DIM <- c(64,64)
mean.coverage.norm <- array(NA,dim=c(prod(DIM),nsim))

# Vector with percentage where we are
perc <- seq(0,1,by=0.1)

# Start for loop
for(i in 1:nsim){
  # Beta image with true value for N subjects
  TrueBetaImage <- array(trueVal,dim=c(DIM,nsub))
  # Add white noise
  BetaImage <- TrueBetaImage + rnorm(n=prod(c(DIM,nsub)),trueVal,wstudSD)
  # Calculate variance
  varBeta <- apply(array(BetaImage,dim=c(prod(DIM),nsub)),1,var)
  # Calculate mean beta value
  meanBeta <- apply(array(BetaImage,dim=c(prod(DIM),nsub)),1,mean)
  # Calculate T-map
  TMap <- meanBeta/sqrt(varBeta/nsub)

  # Now calculate confidence interval
  CI.upper.norm <- meanBeta + (qt(0.975,nsub-1) * sqrt(varBeta/nsub))
  CI.lower.norm <- meanBeta - (qt(0.975,nsub-1) * sqrt(varBeta/nsub))

  # Check if true value is in the CI
  mean.coverage.norm[,i] <- ifelse(trueVal > CI.lower.norm & trueVal < CI.upper.norm, 1, 0)
}
# Check mean in voxel and overall mean.
mean.coverage.norm
dim(mean.coverage.norm)
mean.coverage.norm.vox <- apply(mean.coverage.norm,1,mean)
mean(mean.coverage.norm.vox)

  # As this takes quit long, let us save this.
  save(mean.coverage.norm, file=paste(wd,'/RObjects/',date,'-mean_coverage_TMap',sep=''))
  # Load in object
  load(paste(wd,'/RObjects/',date,'-mean_coverage_TMap',sep=''))

############
## PLOTTING

# Test levelplot
levelplot(array(mean.coverage.norm.vox,dim=DIM))
# Add histogram
mean.coverage.norm.frame <- data.frame('coverage' = matrix(mean.coverage.norm.vox,ncol=1))

# Arrange in 1 frame
levelplot <- levelplot(array(mean.coverage.norm.vox,dim=DIM),xlab='X',ylab='Y')
hist <- ggplot(mean.coverage.norm.frame, aes(x=coverage)) + geom_histogram() +
scale_x_continuous(name="") +
geom_vline(xintercept=0.95,colour='red') +
ggtitle(paste('Coverage for T-value in 64x64 voxels. Mean = ', round(mean(mean.coverage.norm.vox),3),sep='')) +
theme(plot.title = element_text(lineheight=.6, face="bold"))

grid.arrange(levelplot,hist,ncol=2)





########################################
## CI for group maps after transformation to ES. Normal distribution approximation.
#   * 1.000 simulations
#   * Create 64x64 beta-images for N subjects, only 1 study.
#   * Pool these with ordinary OLS pooling method: T-maps
#   * Transform to ES with formula used in FixRan study.
#   * Construct CI in each voxel, around ES, with normal approximation from Borenstein et al. (2009)
#   * Check coverage in each voxel over all simulations


# Number of simulations
nsim <- 1000
# Get vector for mean coverage within simulation
mean.coverage.norm <- array(NA,dim=c(prod(DIM),nsim))

# Start for loop
for(i in 1:nsim){
  # Keeping track of progress
  if(c(i/nsim) %in% perc) print(paste('At ', i/nsim*100, '%', sep=''))
  # Beta image with true value for N subjects
  TrueBetaImage <- array(trueVal,dim=c(DIM,nsub))
  # Add white noise
  BetaImage <- TrueBetaImage + rnorm(n=prod(c(DIM,nsub)),trueVal,wstudSD)
  # Calculate variance
  varBeta <- apply(array(BetaImage,dim=c(prod(DIM),nsub)),1,var)
  # Calculate mean beta value
  meanBeta <- apply(array(BetaImage,dim=c(prod(DIM),nsub)),1,mean)
  # Calculate T-map
  TMap <- meanBeta/sqrt(varBeta/nsub)

  # Transform to an ES using hedgeG function
  HedgeG <- apply(matrix(TMap,ncol=1),1,FUN=hedgeG,N=nsub)
    # Calculate variance of ES
  VarianceHedgeG <- apply(matrix(HedgeG,ncol=1),1,FUN=varHedge,N=nsub)


  # Now calculate confidence intervals for ES based on assumption of normal distributed ES
    CI.upper.norm <- matrix(HedgeG,ncol=1) + (1.96 * sqrt(matrix(VarianceHedgeG,ncol=1)))
    CI.lower.norm <- matrix(HedgeG,ncol=1) - (1.96 * sqrt(matrix(VarianceHedgeG,ncol=1)))
    # Other possibilities...?

  # Now check coverage and save into vector
  mean.coverage.norm[,i] <- ifelse(trueVal > CI.lower.norm & trueVal < CI.upper.norm, 1, 0)
}
# What is the mean?
mean.coverage.norm
dim(mean.coverage.norm)
mean.coverage.norm.vox <- apply(mean.coverage.norm,1,mean);mean.coverage.norm.vox
mean(mean.coverage.norm.vox)

  # As this takes quit long, let us save this.
  save(mean.coverage.norm, file=paste(wd,'/RObjects/',date,'-mean_coverage_norm_ES',sep=''))
  # Load in object
  load(paste(wd,'/RObjects/',date,'-mean_coverage_norm_ES',sep=''))

############
## PLOTTING

# Test levelplot
levelplot(array(mean.coverage.norm.vox,dim=DIM))
# Add histogram
mean.coverage.norm.frame <- data.frame('coverage' = matrix(mean.coverage.norm.vox,ncol=1))

# Arrange in 1 frame
levelplot <- levelplot(array(mean.coverage.norm.vox,dim=DIM),xlab='X',ylab='Y')
hist <- ggplot(mean.coverage.norm.frame, aes(x=coverage)) + geom_histogram() +
scale_x_continuous(name="") +
geom_vline(xintercept=0.95,colour='red') +
ggtitle(paste('Coverage for ES in 64x64 voxels. Mean = ', round(mean(mean.coverage.norm.vox),3),sep='')) +
  theme(plot.title = element_text(lineheight=.6, face="bold"))

grid.arrange(levelplot,hist,ncol=2)




########################################
## CI for group maps after transformation to ES. Normal distribution approximation.
#   * 500 simulations
#   * Create 64x64 beta-images for N subjects and 15 studies.
#   * Pool each study with ordinary OLS pooling method: T-maps.
#   * Transform each study to ES with formula used in FixRan study.
#   * Aggregate studies using fixed effects meta-analysis.
#   * Construct CI in each voxel, around weighted average of ES, with normal approximation from Borenstein et al. (2009)
#   * Check coverage in each voxel over all simulations


# Number of simulations
nsim <- 500
# Get vector for mean coverage within simulation
mean.coverage.norm <- array(NA,dim=c(prod(DIM),nsim))
# number of studies (k)
nstud <- 15

# Vector with percentage where we are
perc <- seq(0,1,by=0.1)

# Start for loop
for(i in 1:nsim){
  if(c(i/nsim) %in% perc) print(paste('At ', i/nsim*100, '%', sep=''))
  # Beta image with true value for N subjects and k studies.
  TrueBetaImage <- array(trueVal,dim=c(DIM,nsub,nstud))
  # Add white noise
  BetaImage <- TrueBetaImage + rnorm(n=prod(c(DIM,nsub,nstud)),trueVal,wstudSD)
  # Now switch to study levels
    # Vector for g-values
    HedgeGStud <- array(NA,dim=c(prod(DIM),nstud))
    # Vector for weights
    weightsStud <- array(NA,dim=c(prod(DIM),nstud))
  for(s in 1:nstud){
    BetaImage.stud <- BetaImage[,,,s]
    # Calculate variance of contrast for each study.
    varBeta <- apply(array(BetaImage.stud,dim=c(prod(DIM),nsub)),1,var)
    # Calculate mean beta value for each study.
    meanBeta <- apply(array(BetaImage.stud,dim=c(prod(DIM),nsub)),1,mean)
    # Calculate T-map for each study.
    TMap <- meanBeta/sqrt(varBeta/nsub)

    # Transform to an ES using hedgeG function, for each study
    HedgeG <- apply(matrix(TMap,ncol=1),1,FUN=hedgeG,N=nsub)
    # Calculate variance of ES
    VarianceHedgeG <- apply(matrix(HedgeG,ncol=1),1,FUN=varHedge,N=nsub)
    # Weights of this study
    weigFix <- 1/VarianceHedgeG

      # Put hedge's G and weights in vector
      HedgeGStud[,s] <- HedgeG
      weightsStud[,s] <- weigFix
    }

    # Now calculate weighted average.
    WeightedAvg <- (apply((HedgeGStud*weightsStud),1,sum))/(apply(weightsStud,1,sum))

    # Calculate variance of weighted average
    varWeightAvg <- 1/apply(weightsStud,1,sum)

    # Now calculate confidence intervals for weighted average based on assumption of normal distributed ES
    CI.upper.norm <- matrix(WeightedAvg,ncol=1) + (1.96 * sqrt(matrix(varWeightAvg,ncol=1)))
    CI.lower.norm <- matrix(WeightedAvg,ncol=1) - (1.96 * sqrt(matrix(varWeightAvg,ncol=1)))

  # Now check coverage and save into vector
  mean.coverage.norm[,i] <- ifelse(trueVal > CI.lower.norm & trueVal < CI.upper.norm, 1, 0)

  }
# Results
mean.coverage.norm
dim(mean.coverage.norm)
mean.coverage.norm.vox <- apply(mean.coverage.norm,1,mean);mean.coverage.norm.vox
mean(mean.coverage.norm.vox)

  # As this takes quit long, let us save this.
  save(mean.coverage.norm, file=paste(wd,'/RObjects/',date,'-mean_coverage_norm_wmean',sep=''))
  # Load in object
  load(paste(wd,'/RObjects/',date,'-mean_coverage_norm_wmean',sep=''))


############
## PLOTTING

# Test levelplot
levelplot(array(mean.coverage.norm.vox,dim=DIM))
# Add histogram
mean.coverage.norm.frame <- data.frame('coverage' = matrix(mean.coverage.norm.vox,ncol=1))

# Arrange in 1 frame
levelplot <- levelplot(array(mean.coverage.norm.vox,dim=DIM),xlab='X',ylab='Y')
hist <- ggplot(mean.coverage.norm.frame, aes(x=coverage)) + geom_histogram() +
scale_x_continuous(name="") +
geom_vline(xintercept=0.95,colour='red') +
ggtitle(paste('Coverage for weighted average in 64x64 voxels. Mean = ', round(mean(mean.coverage.norm.vox),3),sep='')) +
  theme(plot.title = element_text(lineheight=.6, face="bold"))

grid.arrange(levelplot,hist,ncol=2)






########################################
## CI for classical meta-analysis. Standardized mean difference.
#   * 500 simulations
#   * Create 20 studies (each N=30), each looking at difference between group 1 and 2.
#   * Transform each study to ES with formula used in Borenstein.
#   * Aggregate studies using fixed effects meta-analysis.
#   * Construct CI normal distribution from Borenstein et al. (2009)
#   * Check coverage over all simulations

# Reset seed
set.seed(1112)

nsim <- 5000
nsub <- 30
nstud <- 20
trueVal <- 0

CI.upper <- array(NA,dim=nsim)
CI.lower <- array(NA,dim=nsim)
coverage <- array(NA,dim=nsim)

for(i in 1:nsim){
	dStud <- array(NA,dim=nstud)
	varStud <- array(NA,dim=nstud)
	for(s in 1:nstud){
		N1 <- nsub/2
		N2 <- nsub/2
		# Subjects in first group
		X1 <- trueVal + rnorm(n=N1,mean=0,sd=2)
		S2_X1 <- var(X1)
		# Second group
		X2 <- trueVal + rnorm(n=N2,mean=0,sd=2)
		S2_X2 <- var(X2)
		# Sample variance
		Swithin <- sqrt((((N1-1)*S2_X1)+((N2-1)*S2_X2))/(N1 + N2 - 2))
		# Cohen d
		d <- (mean(X1) - mean(X2))/Swithin
		# Variance
		Vd <- (N1+N2)/(N1*N2) + (d^2/(2*(N1+N2)))
		# Save in vector
		dStud[s] <- d
		varStud[s] <- Vd
	}

	# Now meta-analysis
	# Weights
	W <- 1/varStud
	# Weighted mean
	wMean <- sum(W*dStud)/sum(W)
	# Variance of summary effect
	varSum <- 1/sum(W)
	# Confidence interval
	CI.upper[i] <- wMean + (1.96 * sqrt(varSum))
	CI.lower[i] <- wMean - (1.96 * sqrt(varSum))

	# Coverage
	coverage[i] <- ifelse(trueVal >= CI.lower[i] & trueVal <= CI.upper[i], 1, 0)
}


mean(coverage)




########################################
## CI for classical meta-analysis. One sample effect-size
#   * 500 simulations
#   * Create 20 studies (each N=30), each looking at difference between group and 0.
#   * Transform each study to ES with formula used in Borenstein.
#   * Aggregate studies using fixed effects meta-analysis.
#   * Construct CI normal distribution from Borenstein et al. (2009)
#   * Check coverage over all simulations


nsim <- 50000
nsub <- 30
nstud <- 20
trueVal <- 0

dCI.upper <- gCI.upper <- gdCI.upper <- array(NA,dim=nsim)
dCI.lower <- gCI.lower <- gdCI.lower <- array(NA,dim=nsim)
dcoverage <- gcoverage <- gdcoverage <- array(NA,dim=nsim)
dwMean <- gwMean <- gdwMean <- array(NA,dim=nsim)

for(i in 1:nsim){
	dStud <- array(NA,dim=nstud)
	dvarStud <- array(NA,dim=nstud)
	gStud <- array(NA,dim=nstud)
	gvarStud <- array(NA,dim=nstud)
	gdvarStud <- array(NA,dim=nstud)
	for(s in 1:nstud){
		# Subjects in group
		X1 <- trueVal + rnorm(n=nsub,mean=0,sd=2)
		S2_X1 <- var(X1)

		# Cohen d
		d <- (mean(X1))/sqrt(S2_X1)
		# Variance
		Vd <- (1/nsub) + (1 - (gamma((nsub - 2) / 2) / gamma((nsub - 1) / 2))^2 * (nsub - 3) / 2) * d^2

		# Hedge g
		t <- (mean(X1))/sqrt(S2_X1/nsub)
		g <- hedgeG(t,nsub)
		Vg <- varHedge(g,nsub)

		# Hedge g with variance of Radua
		Vgd <- (1/nsub) + (1 - (gamma((nsub - 2) / 2) / gamma((nsub - 1) / 2))^2 * (nsub - 3) / 2) * g^2

		# Save in vector
		dStud[s] <- d
		dvarStud[s] <- Vd
		gStud[s] <- g
		gvarStud[s] <- Vg
		gdvarStud[s] <- Vgd
	}

	# Now meta-analysis
	# Weights
	dW <- 1/dvarStud
	gW <- 1/gvarStud
	gdW <- 1/gdvarStud
	# Weighted means
	dwMean[i] <- sum(dW*dStud)/sum(dW)
	gwMean[i] <- sum(gW*gStud)/sum(gW)
	gdwMean[i] <- sum(gdW*gStud)/sum(gdW)
	# Variance of summary effect
	dvarSum <- 1/sum(dW)
	gvarSum <- 1/sum(gW)
	gdvarSum <- 1/sum(gdW)

	# Confidence interval
	dCI.upper[i] <- dwMean[i] + (1.96 * sqrt(dvarSum))
	dCI.lower[i] <- dwMean[i] - (1.96 * sqrt(dvarSum))
	gCI.upper[i] <- gwMean[i] + (1.96 * sqrt(gvarSum))
	gCI.lower[i] <- gwMean[i] - (1.96 * sqrt(gvarSum))
	gdCI.upper[i] <- gdwMean[i] + (1.96 * sqrt(gdvarSum))
	gdCI.lower[i] <- gdwMean[i] - (1.96 * sqrt(gdvarSum))

	# Coverage
	dcoverage[i] <- ifelse(trueVal >= dCI.lower[i] & trueVal <= dCI.upper[i], 1, 0)
	gcoverage[i] <- ifelse(trueVal >= gCI.lower[i] & trueVal <= gCI.upper[i], 1, 0)
	gdcoverage[i] <- ifelse(trueVal >= gdCI.lower[i] & trueVal <= gdCI.upper[i], 1, 0)

}

mean(dcoverage);mean(gcoverage);mean(gdcoverage)

par(mfrow=c(2,2))
hist(dwMean, main='Cohen d', xlab='Weighted average')
hist(gwMean, main='Hedge g, ours', xlab='Weighted average')
hist(gdwMean, main='Hedge g, Radua', xlab='Weighted average')



########################################
## CI for classical meta-analysis. Based on linear regression.
#		* Linear regression with X = blocked design.
#   * Create 20 studies (each N=30), each looking at difference between group and 0.
#   * Transform each study to ES with formula used in Borenstein.
#   * Aggregate studies using fixed effects meta-analysis.
#   * Construct CI normal distribution from Borenstein et al. (2009)
#   * Check coverage over all simulations


nsim <- 10000
nsub <- 30
nstud <- 20
trueVal <- 0

gCI.upper <- gdCI.upper <- array(NA,dim=nsim)
gCI.lower <- gdCI.lower <- array(NA,dim=nsim)
gcoverage <- gdcoverage <- array(NA,dim=nsim)
gwMean <- gdwMean <- array(NA,dim=nsim)

for(i in 1:nsim){
	dStud <- array(NA,dim=nstud)
	dvarStud <- array(NA,dim=nstud)
	gStud <- array(NA,dim=nstud)
	gvarStud <- array(NA,dim=nstud)
	gdvarStud <- array(NA,dim=nstud)
	Beta <- c()
	for(s in 1:nstud){
		# X values
		X <- sample(c(0:50),nsub)
		# Data
		Y <- trueVal + rnorm(n=length(X),mean=0,sd=3)
		# Fit model
		fit <- lm(Y~X)

		# T-value, transform to hedge G and two types of variance.
		t <- summary(fit)$coefficients['X','t value']
		g <- hedgeG(t,nsub)
		Vg <- varHedge(g,nsub)
		Vgd <- (1/nsub) + (1 - (gamma((nsub - 2) / 2) / gamma((nsub - 1) / 2))^2 * (nsub - 3) / 2) * g^2

		# Save in vector
		gStud[s] <- g
		gvarStud[s] <- Vg
		gdvarStud[s] <- Vgd
	}
	# Now meta-analysis
	# Weights
	gW <- 1/gvarStud
	gdW <- 1/gdvarStud
	# Weighted means
	gwMean[i] <- sum(gW*gStud)/sum(gW)
	gdwMean[i] <- sum(gdW*gStud)/sum(gdW)
	# Variance of summary effect
	gvarSum <- 1/sum(gW)
	gdvarSum <- 1/sum(gdW)

	# Confidence interval
	gCI.upper[i] <- gwMean[i] + (1.96 * sqrt(gvarSum))
	gCI.lower[i] <- gwMean[i] - (1.96 * sqrt(gvarSum))
	gdCI.upper[i] <- gdwMean[i] + (1.96 * sqrt(gdvarSum))
	gdCI.lower[i] <- gdwMean[i] - (1.96 * sqrt(gdvarSum))

	# Coverage
	gcoverage[i] <- ifelse(trueVal >= gCI.lower[i] & trueVal <= gCI.upper[i], 1, 0)
	gdcoverage[i] <- ifelse(trueVal >= gdCI.lower[i] & trueVal <= gdCI.upper[i], 1, 0)

}

mean(gcoverage);mean(gdcoverage)

hist(gwMean)
hist(gdwMean)






########################################
## CI for group maps after low level simulation of fMRI data using neuRosim.
#   * Create 16x16x16 null-images for N subjects and 15 studies.
#   * Each image is created using the same design.
#   * No between study changes in the SNR.
#   * Pool each study with ordinary OLS pooling method: T-maps.
#   * Transform each study to ES with formula used in FixRan study.
#   * Aggregate studies using fixed effects meta-analysis.
#   * Construct several types of CI in each voxel
#   * Check coverage in each voxel over all simulations


###########################
# SEE FILE lowLevelToMeta.R
###########################


################### // ###################
## This code now gets simulated on HPC ##
################### // ###################
nsim <- 8000
DIM <- c(16,16,16)
mean.coverage.norm <- uppr.coverage.norm <- lowr.coverage.norm <-
	mean.coverage.t <- uppr.coverage.t <- lowr.coverage.t <-
	mean.coverage.weightAvg <- uppr.coverage.weightAvg <- lowr.coverage.weightAvg <-
			array(NA,dim=c(prod(DIM),nsim))

trueVal <- 0

####************####
#### Results
####************####
# Load in R objects and calculate coverage.
for(i in 1:nsim){
  # Load in CI.upper.norm and CI.lower.norm
  load(paste(DATAwd[[currentWD]],'/',i,'/',prefix[[currentWD]],'CI.upper.norm_',i,sep=''))
  load(paste(DATAwd[[currentWD]],'/',i,'/',prefix[[currentWD]],'CI.lower.norm_',i,sep=''))
  # Calculate indicator for coverage and save into vector
  mean.coverage.norm[,i] <- ifelse(trueVal >= CI.lower.norm & trueVal <= CI.upper.norm, 1, 0)
  rm(CI.upper.norm,CI.lower.norm)
}
mean.coverage.norm
dim(mean.coverage.norm)
mean.coverage.norm.vox <- apply(mean.coverage.norm,1,mean);head(mean.coverage.norm.vox)
mean(mean.coverage.norm.vox)

  # As this takes some time, let us save this.
  save(mean.coverage.norm, file=paste(wd,'/RObjects/',date,'-mean_coverage_norm_lowToMeta',sep=''))
  # Load in object
  load(paste(wd,'/RObjects/2016-01-18-mean_coverage_norm_lowToMeta',sep=''))


####************####
#### Plotting
####************####
levelplot(array(mean.coverage.norm.vox,dim=DIM))

# Add histogram
mean.coverage.norm.frame <- data.frame('coverage' = matrix(mean.coverage.norm.vox,ncol=1))

# Arrange in 1 frame
levelplot <- levelplot(array(mean.coverage.norm.vox,dim=DIM),xlab='X',ylab='Y')
hist <- ggplot(mean.coverage.norm.frame, aes(x=coverage)) + geom_histogram(binwidth=0.002225) +
scale_x_continuous(name="") +
geom_vline(xintercept=0.95,colour='red') +
ggtitle(label='') +
  theme(plot.title = element_text(lineheight=.6, face="bold"))

grid.arrange(levelplot,hist,ncol=2,
  main=textGrob(paste('Coverage for simulated data from first level to weighted mean average in 16x16x16 voxels. Mean = ', round(mean(mean.coverage.norm.vox),3),sep=''),
  gp=gpar(cex=1.5),vjust=2))


## Plotting the weighted averages
# Load in R objects.
mean.wAvg.norm <- array(NA,dim=c(prod(DIM),nsim))
for(i in 1:nsim){
  # Load in CI.upper.norm and CI.lower.norm
  load(paste(DATAwd[[currentWD]],'/',i,'/',prefix[[currentWD]],'WeightedAvg_',i,sep=''))
  mean.wAvg.norm[,i] <- WeightedAvg
  rm(WeightedAvg)
}

head(mean.wAvg.norm)
dim(mean.wAvg.norm)

# Take some random voxel and look at weighted mean distribution
ToSample <- sample(c(1:dim(mean.wAvg.norm)[1]),20)
par(mfrow=c(4,5))
for(i in 1:length(ToSample)){
  index <- ToSample[i]
  hist(mean.wAvg.norm[index,], main = paste('Voxel ', index,sep=''),xlab='Weighted average',ylab='')
}

# Look at distribution of weighted averages across all voxels (and simulations).
mean.wAvg.norm.vox <- apply(mean.wAvg.norm,1,mean)
  length(mean.wAvg.norm.vox)
mean.wAvg.norm.frame <- data.frame('wAvg' = matrix(mean.wAvg.norm.vox,ncol=1))
ggplot(mean.wAvg.norm.frame, aes(x=wAvg)) + geom_histogram() + geom_density(size=1.4,colour='#016450') +
scale_x_continuous(name="")
# Or density plot
ggplot(mean.wAvg.norm.frame, aes(x=wAvg)) + geom_density(size=1.5,fill='#1c9099',colour='#016450') +
scale_x_continuous(name="")


# Look at weighted mean distribution of voxels with larger coverage
ToLook <- which(mean.coverage.norm.vox>0.986)
length(ToLook)
par(mfrow=c(5,6))
for(i in 1:length(ToLook)){
  index <- ToLook[i]
  hist(mean.wAvg.norm[index,], main = paste('Voxel ', index,sep=''),xlab='Weighted average',ylab='')
}


# Try to make a funnel plot of the individual voxels with their CI bar around
ggplot(mean.wAvg.norm.frame) + geom_points

# Horizontal plot
hor.wAvg <- data.frame('x' = seq(1:length(mean.wAvg.norm.vox)),
                          'y' = mean.wAvg.norm.vox)

ggplot(hor.wAvg, aes(x, y)) + geom_point()
  # Where are the voxels with these values > 0.005
  IDval <- mean.wAvg.norm.frame > .005
  checkVox <- array(0,dim=length(mean.wAvg.norm.vox))
  checkVox[IDval] <- 1
checkVox <- array(checkVox,dim=DIM)
  levelplot(checkVox)
checkVox[4,4,4] <- 5
checkVox[10,10,10] <- 5
checkVox[10,10,11] <- 2.5
  levelplot(checkVox)


