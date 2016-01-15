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
## CI for group maps after low level simulation of fMRI data using neuRosim.
#   * 500 simulations
#   * Create 16x16x16 null-images for N subjects and 15 studies.
#   * Each image is created using the same design.
#   * No between study changes in the SNR.
#   * Pool each study with ordinary OLS pooling method: T-maps.
#   * Transform each study to ES with formula used in FixRan study.
#   * Aggregate studies using fixed effects meta-analysis.
#   * Construct CI in each voxel, around weighted average of ES, with normal approximation from Borenstein et al. (2009)
#   * Check coverage in each voxel over all simulations

## Reset seed
set.seed(11121990)

####************####
#### Global options
####************####
nstud <- 15
nsub <- 15
TR <- 2
nscan <- 200
total <- TR*nscan
on1 <- seq(1,total,40)
on2 <- seq(20,total,40)
onsets <- list(on1,on2)
duration <- list(20,20)
effect.null <- list(0,0)                              ## No effect
effect <- list(1,1) 			                            ## Effect of 1 for designmatrix
DIM <- c(16,16,16)


####************####
#### Subject/Study specific simulation details
####************####
# Subject parameters
TrueLoc1 <- c(4,4,4)
TrueLoc2 <- c(10,10,10)
TrueWhiteNoise <- 0.7
TrueRadius <- 1
locations1 <- locations2 <- weights <- radius <- c()

COPE <- VARCOPE <- TMAP <- array(NA,dim=c(DIM,nsub))

# Loop over nsub to get the weights and the locations
for(s in 1:nsub){
  # Random locations
  loc.tmp1 <- TrueLoc1 + round(rnorm(3,0,2),0)
    while(any(loc.tmp1<1)){
      loc.tmp1 <- TrueLoc1 + round(rnorm(3,0,2),0)
    }
  loc.tmp2 <- TrueLoc2 + round(rnorm(3,0,2),0)
    while(any(loc.tmp2>16)){
      loc.tmp2 <- TrueLoc2 + round(rnorm(3,0,2),0)
    }
  locations1 <- rbind(locations1,loc.tmp1)
  locations2 <- rbind(locations2,loc.tmp2)

  # Radius
  radius.tmp <- TrueRadius + sample(c(0,1),size=2)
  radius <- rbind(radius,'rad' = radius.tmp)

  # Random noise components
  whiteNoise <- round(TrueWhiteNoise + rnorm(1,0,0.5),2)
  while(whiteNoise > 1 || whiteNoise < 0.5){
    whiteNoise <- round(TrueWhiteNoise + rnorm(1,0,0.5),2)
  }
  temporalNoise <- round((1-whiteNoise)/2 + rnorm(1,0,0.25),2)
  while(temporalNoise > c(1-whiteNoise) || temporalNoise < 0){
    temporalNoise <- round(TrueWhiteNoise + rnorm(1,0,0.5),2)
  }
  spatialNoise <- 1-whiteNoise-temporalNoise

  weights <- rbind(weights, c(whiteNoise, temporalNoise,0 ,0, 0, spatialNoise))
  rm(loc.tmp1,loc.tmp2,whiteNoise,temporalNoise,spatialNoise)
}

# Study parameters
SWEIGHTS <- SHEDGE <- SCOPE <- SVARCOPE <- STMAP <- array(NA,dim=c(DIM,nstud))



####************####
#### Design matrices
####************####
# Design Matrices via neuRosim:
#     * We need three design vectors:
#     * The first two have an intercept (needed for analysis).
#        * These will be the two columns (1 -1 contrast) of the design matrix in the analysis.
#     * The third one is used to generate data with a NULL effect.
designC1 <- simprepTemporal(onsets = list(on1), durations = list(duration[[1]]),
                         hrf = "double-gamma", TR = TR, totaltime = total,
                         effectsize = list(effect[[1]]))

designC2 <- simprepTemporal(onsets = list(on2), durations = list(duration[[2]]),
                        hrf = "double-gamma", TR = TR, totaltime = total,
                        effectsize = list(effect[[1]]))

design.null <- simprepTemporal(regions = 2, onsets = onsets, durations = duration,
                         hrf = "double-gamma", TR = TR, totaltime = total,
                         effectsize = effect.null)

# X-matrix in order to fit the model later on (combination of C1 and C2).
x <- fmri.design(matrix(c(simTSfmri(designC1, nscan=nscan, TR=TR, noise="none"),
        simTSfmri(designC2, nscan=nscan, TR=TR, noise="none")),ncol=2),0)



####************####
#### GENERATE DATA: INCLUDE STUDIES --> SUBJECTS
####************####

# For loop over studies
for(t in 1:nstud){
  print(paste('------------------------- STUDY ', t,' -------------------------',sep=''))
  # For loop over nsub
  for(s in 1:nsub){
    print(paste('At subject, ', s, sep=''))
    coordinates <- list(c(locations1[s,]),c(locations2[s,]))
    # Define two regions (which does nothing as there is no effect, )
    regions <- simprepSpatial(regions = 2, coord = coordinates, radius = list(radius[s,1],radius[s,2]), form ="cube", fading = 0)
      rm(coordinates)

    # Weighting structure: white, temporal and spatial noise.
    #   * Order = white, temporal, low-frequency, physyiological, task related and spatial.
    w <- weights[s,]
    # Base value
    base <- 5

    # Actual simulated data
    sim.data <-  simVOLfmri(design=design.null, image=regions, base=base, dim=DIM, SNR=0.5,
                 type ="gaussian", noise= "mixture", spat="gaussRF", FWHM=2, weights=w, verbose = TRUE)
        rm(w)

    # 3D Gaussian Kernel over the 4D data
    fwhm <- 2
    sigma <- fwhm/(sqrt(8)*log(2))
    smoothint <- -round(2*sigma):round(2*sigma)

    for(i in 1:nscan){
     sim.data[,,,i] <- GaussSmoothArray(sim.data[,,,i], voxdim=c(1,1,1), ksize = length(smoothint), sigma = diag(sigma, 3))
    }

    # Reshape the data for fMRI analysis and make it the correct class
    datafmri <- list(ttt=writeBin(c(sim.data), raw(),4), mask=array(1,dim=DIM), dim = c(DIM, nscan))
    class(datafmri) <- "fmridata"
      rm(sim.data)
    ####************####
    #### ANALYZE DATA
    ####************####


    # Fitting GLM model: estimated AR(1)-coefficients are used to whiten data, may produce warnings because data is pre-smoothed.
    model <- fmri.lm(datafmri,x, actype = "accalc", keep="all",contrast=c(1,-1))

    # Estimated contrast of parameter beta's from model
    COPE.sub <- model$cbeta
      COPE[,,,s] <- COPE.sub
    VARCOPE.sub <- model$var
      VARCOPE[,,,s] <- VARCOPE.sub
    # Constructing t-map
    TMAP.sub <- array(c(COPE.sub)/sqrt(c(VARCOPE.sub)), dim=c(DIM))
      TMAP[,,,s] <- TMAP.sub

      # Need a plot to check?
      PLOT <- FALSE
      if(isTRUE(PLOT)){
        levelplot(TMAP[,,,1])
        PVAL <- array(1-pt(TMAP,df=nscan-1),dim=DIM)
        levelplot(PVAL[,,,1])
        IDsign <- PVAL[,,,1] < 0.05
        PVAL[IDsign] <- 1
        PVAL[!IDsign] <- 0
        levelplot(PVAL)
      }
      rm(COPE.sub,VARCOPE.sub,TMAP.sub)
  }

  ####************####
  #### GROUP ANALYSIS
  ####************####

  # Group COPE (average)
  GCOPE <- apply(COPE,c(1,2,3),mean)

  # Now we will do the OLS estimation of the variance
  GVARCOPE <- apply(COPE,c(1,2,3),var)

  # TMAP
  GTMAP <- GCOPE/sqrt(GVARCOPE/(nsub-1))


  ####************####
  #### TRANSFORM TO ES
  ####************####
  # Transform to an ES using hedgeG function, for each study
  HedgeG <- apply(matrix(GTMAP,ncol=1),1,FUN=hedgeG,N=nsub)
  # Calculate variance of ES
  VarianceHedgeG <- apply(matrix(HedgeG,ncol=1),1,FUN=varHedge,N=nsub)
  # Weights of this study
  weigFix <- 1/VarianceHedgeG

    # Put GCOPE, GVARCOPE, GTMAP, hedge's G and weights in vector
    SCOPE[,,,t] <- GCOPE
    SVARCOPE[,,,t] <- GVARCOPE
    STMAP[,,,t] <- GTMAP
    SHEDGE[,,,t] <- HedgeG
    SWEIGHTS[,,,t] <- weigFix

    rm(GCOPE,GVARCOPE,HedgeG,weigFix)
    gc(verbose = FALSE)
}


# Re-format arrays
SHEDGE.mat <- matrix(SHEDGE,ncol=nstud)
SWEIGHTS.mat <- matrix(SWEIGHTS,ncol=nstud)

# Now calculate weighted average.
WeightedAvg <- (apply((SHEDGE.mat*SWEIGHTS.mat),1,sum))/(apply(SWEIGHTS.mat,1,sum))

# Calculate variance of weighted average
varWeightAvg <- 1/apply(SWEIGHTS.mat,1,sum)

# Now calculate confidence intervals for weighted average based on assumption of normal distributed ES
CI.upper.norm <- matrix(WeightedAvg,ncol=1) + (1.96 * sqrt(matrix(varWeightAvg,ncol=1)))
CI.lower.norm <- matrix(WeightedAvg,ncol=1) - (1.96 * sqrt(matrix(varWeightAvg,ncol=1)))





# Now check coverage and save into vector
mean.coverage.norm[,i] <- ifelse(trueVal > CI.lower.norm & trueVal < CI.upper.norm, 1, 0)

















