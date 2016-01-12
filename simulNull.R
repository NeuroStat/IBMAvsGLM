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

# Set starting seed
set.seed(11121990)

# Set WD
wd <- "/Users/hanbossier/Dropbox/PhD/PhDWork/Meta Analysis/R Code/Studie_Simulation/SimulationGit"
setwd(wd)


# Load in libraries
library(lattice)
library(gridExtra)
library(oro.nifti)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(Hmisc)
library(fslr)
library(bootES)

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







########################################
## CI for group maps after transformation to ES. Normal distribution approximation.
#   * 1.000 simulations
#   * Create 64x64 beta-images for N subjects and 15 studies.
#   * Pool each study with ordinary OLS pooling method: T-maps.
#   * Transform each study to ES with formula used in FixRan study.
#   * Aggregate studies using fixed effects meta-analysis.
#   * Construct CI in each voxel, around weighted average of ES, with normal approximation from Borenstein et al. (2009)
#   * Check coverage in each voxel over all simulations


# Number of simulations
nsim <- 1000
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













