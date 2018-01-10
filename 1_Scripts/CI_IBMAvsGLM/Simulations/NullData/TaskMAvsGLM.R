####################
#### TITLE:     MA vs Mixed Model GLM using FSL's FLAME on activation data.
#### Contents:
####
#### Source Files:
#### First Modified: 03/05/2017
#### Notes:
#################



##
###############
### Notes
###############
##

# Now simulate activation data with realistic noise, sample sizes, etc, ...
# In this scenario, we will compare the outcome of transforming the second level GLM to an ES and execute a meta-analsysis with the scenario in which you use a third level GLM.
# MEASURES:
#   * CI coverage
#   * Standardized bias


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
  # K'th simulation
  K <- try(as.numeric(as.character(input)[1]),silent=TRUE)
  # Which scenario
  SCEN <- try(as.numeric(as.character(input)[2]),silent=TRUE)
  # Which machine
  MACHINE <- try(as.character(input)[3],silent=TRUE)
    # If no machine is specified, then it has to be this machine!
    if(is.na(MACHINE)){
      MACHINE <- 'MAC'
      K <- 1
      SCEN <- 1
    }
  # DataWrite directory: where all files are written to
  DataWrite <- try(as.character(input)[4],silent=TRUE)

# Set starting seed: it is the product of the amount of voxels, the number of studies and the number of subjects!
starting.seed <- 36865*K
set.seed(starting.seed)


# Set WD
if(MACHINE=='HPC'){
  wd <- '/user/scratch/gent/gvo000/gvo00022/vsc40728/IBMAvsMA'
}
if(MACHINE=='MAC'){
  wd <- '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/MAvsIBMA'
}
setwd(wd)

# Give path to FSL
if(MACHINE=='HPC'){
  fslpath <- ''
}
if(MACHINE=='MAC'){
  fslpath <- '/usr/local/fsl/bin/'
}

# DataWrite if machine = Mac
if(MACHINE=='MAC'){
  DataWrite <- '~/Desktop/IBMAtmp'
}


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
if(MACHINE == 'MAC'){
  source('~/Dropbox/PhD/PhDWork/Meta\ Analysis/R\ Code/Studie_FixRan/FixRanStudyGit.git/Development/functions.R')
}
if(MACHINE == 'HPC'){
  source('/user/scratch/gent/gvo000/gvo00022/vsc40728/Simulation/functions.R')
}



##
###############
### Simulation steps
###############
##

####************####
#### Global options
####************####
TR <- 2
nscan <- 200
total <- TR*nscan
on1 <- seq(1,total,40)              # Block design: 20 sec ON/OFF
onsets <- list(on1)
duration <- list(20)
effect <- list(1) 			            # Effect of 1 for designmatrix
# To find 3D grids that can be plotted as 2D, find the integers in:
  # ((5:250)**2)**(1/3)
DIM <- c(36,36,36)


#constant parameters
ES1 <- 0.01
ES2 <- 0.02
ES3 <- 0.03
ES4 <- 0.04
base <- 100
actES1 <- ES1*base
actES2 <- ES2*base
actES3 <- ES3*base
actES4 <- ES4*base


A <- array(0, dim = DIM)
  A[c(2:3), c(2:3), c(2:3)] <- 1
B <- matrix(A, ncol = 8)
B[c(3:5),c(3:5)] <- 1
matrix(B, ncol = 8)
array(B, dim = DIM)

####************####
#### Scenario specific simulation details
####************####

# Let us randomly draw sample sizes for the individual studies from the Poldrack et al. paper
# For this, source the data from David et al., as well as the added data
# scraping Neurosynth

"https://raw.githubusercontent.com/poldracklab/ScanningTheHorizon/master/SampleSize/neurosynth_study_data.txt"

"https://raw.githubusercontent.com/poldracklab/ScanningTheHorizon/master/SampleSize/david_sampsizedata.txt"





tbeta <- function(volume,design){

  #volume <- sim.data ; array with smoothed signal
  #design <- pred    ; design

  time <- length(design) #length of time series: 400
  dimens <- dim(volume) #dimension of signal array (in the example: 40 40 40 400)

  varact <- c(var(design)) # variance of x
  meanact <- c(mean(design)) # average of x
  # for every voxel, we estimate the covariance between the design (x) and the signal (y). Then divided by variance of design (x)
  b1 <- apply(volume,c(1,2,3),cov,y=design)/varact

}
dim(sim.data)
b1Jas <- tbeta(volume = sim.data, design = x)
data.frame(b1, c(b1Jas))

design <- x
volume <- sim.data

tsebeta <- function(b1, volume, design){

  time <- length(design) #length of time series: 400
  dimens <- dim(volume) #dimension of signal array (in the example: 40 40 40 400)

  varact <- c(var(design)) # variance of x
  meanact <- c(mean(design)) # average of x

  # mean signal for every voxel - b1 * average of x, the design
  b0 <- apply(volume,c(1,2,3),mean)-b1Jas*meanact
  # in every voxel: b1 * value of the design to calculate predicted value under the model
  predact <- array(array(b1Jas,dim=c(dimens[1],dimens[2],dimens[3],1)) %*% array(design,dim=c(1,time)),dim=dimens)
  # same, but with b0 added
  pred <- array(rep(b0,time),dim=c(dimens[1],dimens[2],dimens[3],time)) + predact
  # calculation of SSE (for all voxels together)
  help <- (pred-volume)^2
  # se2: calculated for every voxel through SSE/(n - 2)
  se2 <- apply(help,c(1,2,3),sum)/(time-2)
  # denominator for variance of b1
  help2 <- sum((design-mean(design))^2)
  # standard error of b1
  sb1 <- sqrt(se2/help2)

}

seBetaJas <- tsebeta(b1Jas, volume = sim.data, design = x)

# Fitting GLM model.
model.lm <- lm(Y.data ~ x)
b1 <- coef(model.lm)['x',]
  # Estimate residual (we need to extend the design matrix with an intercept)
  xIN <- cbind(1,x)
  BETA <- coef(model.lm)
  res <- (t(Y.data - xIN %*% BETA) %*% (Y.data - xIN %*% BETA))/(nscan - 2)
  res <- diag(res)
  # Contrast: not interested in intercept
  CONTRAST <- matrix(c(0,1),nrow=1)

VARCOPE <- CONTRAST %*% (solve(t(xIN) %*% xIN)) %*% t(CONTRAST) %*% res
data.frame( seBeta = c(sqrt(VARCOPE)),
          seBetaJas = c(seBetaJas))

tvalR <- c()
for(i in 1:64){
  tvalR <- c(tvalR,
    summary(lm(Y.data[,i] ~ x))[['coefficients']][2,3]
    )
}

sqrt(diag(solve(t(xIN) %*% xIN) * c(res)))
summary(model.lm)

data.frame(tValJas = c(b1Jas / seBetaJas),
          tValHan = c(b1 / sqrt(VARCOPE)),
          tValR = tvalR)


# Actual simulated data
w <- c(0.5,0.5,0,0,0,0)
sim.data <- simVOLfmri(design=design.null, image=regions, base=base, dim=DIM, SNR=0.5,
             type ="gaussian", noise= "mixture", spat="gaussRF", FWHM=2, weights=w, verbose = TRUE)
  # Transform it to correct dimension (Y = t x V)
  Y.data <- t(matrix(sim.data,ncol=nscan))

plot(Y.data[,1], type = 'l')
ar(Y.data[,1], AIC = FALSE, order.max = 1)
ar.ols(Y.data[,1])

plot(sunspot.year, type = 'l')
plot(lh, type = 'l')
ar(lh)
ar(sunspot.year)

ar(Y.data[,1])$ar


solve((t(xIN) %*% xIN)) %*% t(xIN) %*% Y.data[,1]

cov(x, Y.data[,1]) / var(x)
model.lm <- lm(Y.data ~ x)
b1 <- coef(model.lm)['x',]
b1[1]


times <- 1:5
rho <- 0.5
sigma <- 2
###############
H <- abs(outer(times, times, "-"))
V <- sigma * rho^H
p <- nrow(V)
V[cbind(1:p, 1:p)] <- V[cbind(1:p, 1:p)] * sigma
V

##########################################################################################
AR <- ar(Y.data[,1])$ar
H <- abs(outer(1:nscan, 1:nscan, "-"))
V <- 1 * AR^H

solve((t(xIN) %*% solve(V) %*% xIN)) %*% t(xIN) %*% solve(V) %*% Y.data[,1]
cov(x, Y.data[,1]) / var(x)
##########################################################################################





#constant parameters
es1 <- 0.01
es2 <- 0.02
base <- 100
actes1 <- es1*base
actes2 <- es2*base

fwhm <- 8
sigma <- fwhm/sqrt(8*log(2))
width <- 5


#Image characteristics
dim1 <- 64
dim2 <- 64
dim3 <- 40
imdim <- c(dim1,dim2,dim3) #image dimensions
voxdim <- c(3.5,3.5,3.51) #voxelsize
ext <- 4
TR <- 2
nregio <- 2
coord1 <- c(22,17,11)
coord2 <- c(40,18,11)
coord <- list(coord1,coord2)

#Mask
noisemask <- array(1,dim=imdim)
nobrain <- which(noisemask == 0)
nvoxels <- sum(noisemask==1)



#CREATING GROUND TRUTH

regions1 <- simprepSpatial(regions = 1,
          coord = coord1, radius = ext, form = "sphere", fading = 0)
truthdesign <- simprepTemporal(totaltime = 1, regions = 1, onsets=1, effectsize = 1, durations=1, TR=1, acc=0.1)
truth1 <- simVOLfmri(design = truthdesign, image = regions1,
          dim = imdim, SNR = 1, noise = "none", template = noisemask)[,,,1]
truth1 <- ifelse(truth1 > 0, 1, 0)

regions2 <- simprepSpatial(regions=1, coord=coord2, radius=ext, form="sphere", fading=0)
truthdesign <- simprepTemporal(1,1,onsets=1,effectsize = 1, durations=1,TR=1,acc=0.1)
truth2 <- simVOLfmri(design=truthdesign,image=regions2,dim=imdim,SNR=1,noise="none", template=noisemask)[,,,1]
truth2 <- ifelse(truth2 > 0, 1, 0)

truth <- truth1 + truth2
truth <- ifelse(truth > 0, 1, 0)
notruth <- which(truth==0)
actvox <- sum(truth==1)
levelplot(truth)


## CREATING DESIGN AND SIMULATED fMRI TIME SERIES

#putting in the temporal signal: block design 20s ON/OFF
total.time <- nscan * TR
dur <- 20
onsets <- seq(1, total.time, 40)
# Generating a design matrix
X <- simprepTemporal(total.time,1,onsets=onsets,effectsize = 100,durations=dur,TR=TR,acc=0.1, hrf="double-gamma")
# Generate time series for ONE active voxel
pred <- simTSfmri(design=X, base=100, SNR=1, noise="none", verbose=FALSE)

design1 <- es1 * (pred-base) + base
design2 <- es2 * (pred-base) + base


#creating the signal in the anatomical mask
smtruth1 <- GaussSmoothArray(truth1,voxdim,ksize=width,sigma=diag(sigma,3))
smtruth2 <- GaussSmoothArray(truth2,voxdim,ksize=width,sigma=diag(sigma,3))
smtruth <- smtruth1 + smtruth2
rawsignal1 <- smtruth1 %o% design1
dim(rawsignal1)
rawsignal2 <- smtruth2 %o% design2
rawsignal <- rawsignal1 + rawsignal2

signal <- array(NA,dim=c(imdim,nscan))

for (p in 1:nscan) { # Make sure the smoothed signal remains within the anatomical region

    slice <- rawsignal[,,,p]
    slice[notruth] <- 0
    signal[,,,p] <- slice
    rm(slice)

}


# creating gaussian noise in the brain mask
sigmnoise <- 2
noisim <- array(rnorm(prod(imdim)*nscan,0,sigmnoise),dim=c(imdim,nscan))
snoisim <- GaussSmoothArray(noisim, voxdim=voxdim, ksize = width, sigma = diag(sigma,3))

gaussnoise <- array(NA,dim=c(imdim,nscan))
for (p in 1:nscan) { # Make sure the smoothed signal remains within the anatomical region

slice <- snoisim[,,,p]
slice[nobrain] <- 0
gaussnoise[,,,p] <- slice
rm(slice)

}


#creating the final image
data <- gaussnoise + signal

b1 <- tbeta(data,pred)
sb1 <- tsebeta(b1,data,pred)
tmap <- b1/sb1

# Grand mean scaling
GrandMean <- mean(data, na.rm = TRUE)
for(i in 1:nscan){
  data[,,,i] <- data[,,,i] - GrandMean + (100**2)
}
mean(data)

b1_scale <- tbeta(data,pred)
sb1_scale <- tsebeta(b1_scale,data,pred)
tmap_scale <- b1_scale/sb1_scale

summary(b1 - b1_scale)
summary(sb1 - sb1_scale)
summary(tmap - tmap_scale)
summary(tmap)
summary(tmap_scale)

tmap[is.na(tmap)] <- 0
levelplot(tmap)


x <- round(runif(n = 30, min = 10,max = 50))
y <- x + rnorm(n = 30, sd = 20)

cov(x,y)
var(x)
cov(x,y)/var(x)
cor(x,y)

x2 <- x+100000
cov(x2,y)
var(x2)
cov(x2,y)/var(x2)
cor(x2,y)

var(x)
var(x+10)
data.frame(x, (x*10))
data.frame(x, (x+270))
xTimes10 <- x*10
xPlus300 <- x+300
plot(y ~ x, type = 'p', col = 'black', pch = 15,xlim = c(10,500))
points(y ~ xTimes10, type = 'p', col = 'red', pch = 15)
points(y ~ xPlus300, type = 'p', col = 'green', pch = 15)


# Test featquerry
meanFunc <- readNIfTI("/Volumes/2_TB_WD_Elements_10B8_Han/PhD/HCP/Mixed_Effects.gfeat/1_8.gfeat/cope1.feat/mean_func.nii.gz")[,,]
meanFunc <- readNIfTI("/Volumes/2_TB_WD_Elements_10B8_Han/PhD/HCP/80_language/1/MNINonLinear/Results/tfMRI_LANGUAGE/tfMRI_LANGUAGE_hp200_s4.gfeat/mean_func.nii.gz")[,,]
meanFunc <- readNIfTI("/Volumes/2_TB_WD_Elements_10B8_Han/PhD/HCP/80_language/1/MNINonLinear/Results/tfMRI_LANGUAGE/tfMRI_LANGUAGE_hp200_s4.gfeat/cope1.feat/mean_func.nii.gz")[,,]
meanFunc <- readNIfTI("/Volumes/2_TB_WD_Elements_10B8_Han/PhD/HCP/80_language/1/MNINonLinear/Results/tfMRI_LANGUAGE/tfMRI_LANGUAGE_hp200_s4.gfeat/cope1.feat/mean_func.nii.gz")[,,]
mean(meanFunc)


# Note: featquerry command =
    # featquery 1 "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/HCP/Mixed_Effects.gfeat/25_32.gfeat/cope1.feat" 1 stats/cope1 featquery -p -w mask.nii.gz
# Run at pwd:
    # /Volumes/2_TB_WD_Elements_10B8_Han/PhD/HCP/Mixed_Effects.gfeat/25_32.gfeat/cope1.feat
# Result is max change of 3.273%

meanFuncTest <- readNIfTI("/Volumes/2_TB_WD_Elements_10B8_Han/PhD/HCP/Mixed_Effects.gfeat/25_32.gfeat/cope1.feat/mean_func.nii.gz")[,,]
cope1 <- readNIfTI("/Volumes/2_TB_WD_Elements_10B8_Han/PhD/HCP/Mixed_Effects.gfeat/25_32.gfeat/cope1.feat/stats/cope1.nii.gz")[,,]
max(cope1) / 3.273  * 100
max(cope1) / 5535.542  * 100
mean(meanFuncTest[which(mask==1)])
max(cope1) / 9121.696 * 100

# there are 8 subjects, hence '8 values of the maximum cope 1'

tmpts <- readNIfTI("/Volumes/2_TB_WD_Elements_10B8_Han/PhD/HCP/Mixed_Effects.gfeat/25_32.gfeat/cope1.feat/featquery/tmpts.nii.gz")[,,,]
mask <- readNIfTI("/Volumes/2_TB_WD_Elements_10B8_Han/PhD/HCP/Mixed_Effects.gfeat/25_32.gfeat/cope1.feat/featquery/mask.nii.gz")[,,]
apply(tmpts, c(1,2,3,4), mean)
TMPT1 <- tmpts[,,,1]
mean(TMPT1[which(mask == 1)])

# This is file with all % BOLD changes
tmp <- readNIfTI("/Volumes/2_TB_WD_Elements_10B8_Han/PhD/HCP/Mixed_Effects.gfeat/25_32.gfeat/cope1.feat/featquery/tmp.nii.gz")[,,]
max(tmp)

# It is retreived by: cope1 * 128.05834375 / mean_func
  # This 128.05... changes in other analyses. Not sure, but should come from design or something? Or amount of subjects as it approaches to 100. Or correlation between EVs?
  # See:
  # /Volumes/2_TB_WD_Elements_10B8_Han/PhD/EklundReplication/1/NSTUD_1.gfeat/cope1.feat/featquery


