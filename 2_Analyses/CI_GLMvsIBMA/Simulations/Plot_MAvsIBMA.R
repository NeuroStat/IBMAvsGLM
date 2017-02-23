####################
#### TITLE:     Plotting results of the MAvsIBMA.R data.
#### Contents:
####
#### Source Files:
#### First Modified: 25/02/2016
#### Notes:
#################



##
###############
### Notes
###############
##

# Blocked design for individual subjects.
# One condition.
# These N subjects are pooled using FLAME pooling
# The resulting images are converted to either Hedges' g and pooled using fixed effects meta-analysis.
# OR using 3e level GLM with again FLAME1.



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
  'Take[MAvsIBMA]' = "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/MAvsIBMA"
	)
NUMDATAwd <- length(DATAwd)
currentWD <- 1

# Number of conficence intervals
CIs <- c('MA-weightVar','GLM-t')
NumCI <- length(CIs)

# Number of executed simulations
nsim.tmp <- matrix(c(
                1,3000
              ), ncol=2, byrow=TRUE)
nsim <- nsim.tmp[currentWD,2]


# Number of subjects and studies
nsub <- 100
nstud <- 5

# Dimension of brain
DIM.tmp <- array(NA, dim=c(NUMDATAwd,3))
	DIM.tmp[c(1),] <- c(4,4,4)
DIM <- DIM.tmp[currentWD,]

# True value
trueVal <- 0

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


#### First load in all the data
AllData <- c()
# Load in the data
for(i in 1:nsim){
  load(paste(DATAwd[[currentWD]],'/',i,'/SCEN_1/ObjectsMAvsIBMA_',i,sep=''))
  AllData <- c(AllData, matrix(unlist(ObjectsMAvsIBMA), ncol=1))
}
# Make nsim number of columns
AllData <- matrix(AllData,ncol=nsim)

# Load the naming structure of the data
load(paste(DATAwd[[currentWD]],'/1/SCEN_1/ObjectsMAvsIBMA_1',sep='')); objects <- names(ObjectsMAvsIBMA); rm(ObjectsMAvsIBMA)
OBJ.ID <- c(rep(objects[!objects %in% c("STHEDGE","STWEIGHTS")], each=prod(DIM)), rep(c("STHEDGE","STWEIGHTS"), each=c(prod(DIM)*nstud)))

# Does dimension of AllData match the lenght of the OBJ.ID?
dim(AllData)[1]==length(OBJ.ID)




#########################################################
###################### CI COVERAGE ######################
#########################################################

# Calculate coverage
mean.coverage.weightVar.MA <-
  mean.coverage.t.IBMA <-
  array(NA,dim=c(prod(DIM),1))
mean.coverage.weightVar.MA[,1] <- indicating(UPPER = AllData[which(OBJ.ID=='CI.MA.upper.weightVar'),],LOWER = AllData[which(OBJ.ID=='CI.MA.lower.weightVar'),],trueVal = trueVal)
mean.coverage.t.IBMA[,1] <- indicating(UPPER = AllData[which(OBJ.ID=='CI.IBMA.upper.t'),],LOWER = AllData[which(OBJ.ID=='CI.IBMA.lower.t'),],trueVal = trueVal)


# Put the 2 coverages in a list
mean.coverages <- list('MA' = mean.coverage.weightVar.MA,'IBMA' = mean.coverage.t.IBMA)

# Mean over all voxels
CI.coverages <- data.frame(
	'Mean' = matrix(sapply(mean.coverages, FUN=function(...){apply(...,2,mean)}),ncol=1),
	'SD' = matrix(sapply(mean.coverages, FUN=function(...){apply(...,2,sd)}),ncol=1),
	'CI' = factor(CIs, levels=CIs, labels=CIs)
	)



#########################################################
####################### CI LENGTH #######################
#########################################################

# Calculate CI length
mean.length.weightVar.MA <-
mean.length.t.IBMA <-
array(NA,dim=c(prod(DIM),1))

mean.length.weightVar.MA[,1] <- apply(AllData[which(OBJ.ID=='CI.MA.upper.weightVar'),] - AllData[which(OBJ.ID=='CI.MA.lower.weightVar'),],1,mean)
mean.length.t.IBMA[,1] <- apply(AllData[which(OBJ.ID=='CI.IBMA.upper.t'),] - AllData[which(OBJ.ID=='CI.IBMA.lower.t'),],1,mean)


# Put the 2 lengths in a list
mean.lengths <- list('MA' = mean.length.weightVar.MA,'IBMA' = mean.length.t.IBMA)

# Average over all voxels
CI.lengths <- data.frame(
	'Mean' = matrix(sapply(mean.lengths, FUN=function(...){apply(...,2,mean)}),ncol=1),
	'SD' = matrix(sapply(mean.lengths, FUN=function(...){apply(...,2,sd)}),ncol=1),
	'CI' = factor(CIs, levels=CIs, labels=CIs)
	)



#########################################################
################### STANDARDIZED BIAS ###################
#########################################################

MA.SDBETA <- apply(AllData[which(OBJ.ID=='MA.WeightedAvg'),],1,sd)
MA.MEANBETA <- apply(AllData[which(OBJ.ID=='MA.WeightedAvg'),],1,mean)

IBMA.SDBETA <- apply(AllData[which(OBJ.ID=='IBMA.COPE'),],1,sd)
IBMA.MEANBETA <- apply(AllData[which(OBJ.ID=='IBMA.COPE'),],1,mean)

mean.bias.MA <- matrix(((abs(MA.MEANBETA)-trueVal)/(MA.SDBETA))*100,ncol=1)
mean.bias.IBMA <- matrix(((abs(IBMA.MEANBETA)-trueVal)/(IBMA.SDBETA))*100,ncol=1)

# Put the 2 bias values in a list
mean.bias <- list('MA' = mean.bias.MA,'IBMA' = mean.bias.IBMA)

# Average over all voxels
CI.bias <- data.frame(
	'Mean' = matrix(sapply(mean.bias, FUN=function(...){apply(...,2,mean)}),ncol=1),
	'SD' = matrix(sapply(mean.bias, FUN=function(...){apply(...,2,sd)}),ncol=1),
	'CI' = factor(CIs, levels=CIs, labels=c('MA', 'GLM'))
	)


##
###############
### Plotting: coverage and weighted averages
###############
##

# Function for levelplot
ValuesOnLevelPlot2D <- function(x, y, z, ...) {
    panel.levelplot(x,y,z,...)
    panel.text(x, y, round(z,3),col='red')
}

# Plotting the coverages
CCI1 <- levelplot(array(mean.coverages[['MA']], dim=rep(sqrt(prod(DIM)),2)),
      col.regions = gray(0:100/100), at=seq(0.75,1,by=0.05), main='Meta-Analysis',xlab='',ylab='',
        colorkey = list(space = "bottom"),
               panel=ValuesOnLevelPlot2D)
CCI2 <- levelplot(array(mean.coverages[['IBMA']], dim=rep(sqrt(prod(DIM)),2)),
      col.regions = gray(0:100/100), at=seq(0.75,1,by=0.05), main='3 level GLM',xlab='',ylab='',
              colorkey = list(space = "bottom"),
              panel=ValuesOnLevelPlot2D)
grid.arrange(CCI1,CCI2,nrow=1,top = textGrob('CI - Coverage of each voxel over 3000 simulations.', gp=gpar(fontsize=20,font=1)))


# Plotting the lengths
LCI1 <- levelplot(array(mean.lengths[['MA']], dim=rep(sqrt(prod(DIM)),2)),
      col.regions = gray(100:0/100), at=seq(0,1,by=0.05), main='Meta-Analysis',xlab='',ylab='',
        colorkey = list(space = "bottom"),
               panel=ValuesOnLevelPlot2D)
LCI2 <- levelplot(array(mean.lengths[['IBMA']], dim=rep(sqrt(prod(DIM)),2)),
      col.regions = gray(100:0/100), at=seq(0,1,by=0.05), main='3 level GLM',xlab='',ylab='',
              colorkey = list(space = "bottom"),
              panel=ValuesOnLevelPlot2D)
grid.arrange(LCI1,LCI2,nrow=1,top = textGrob('CI - Length of each voxel over 3000 simulations.', gp=gpar(fontsize=20,font=1)))


# Plotting first the distribution of either weighted average, or COPE-value
WA_density <- ggplot(data.frame('value' = MA.MEANBETA), aes(x=value)) + geom_density(fill='#1b9e77')
COPE_density <- ggplot(data.frame('value' = IBMA.MEANBETA), aes(x=value)) + geom_density(fill='#7570b3')
grid.arrange(WA_density, COPE_density, nrow=1, top = textGrob('Average density (over all simulations) of weighted average/COPE', gp=gpar(fontsize=20,font=1)))

# Plotting the standardized bias
BCI1 <- levelplot(array(mean.bias[['MA']], dim=rep(sqrt(prod(DIM)),2)),
      col.regions = gray(100:0/100), at=seq(0,ceiling(max(unlist(mean.bias))),length.out=100), main='Meta-Analysis',xlab='',ylab='',
        colorkey = list(space = "bottom"),
               panel=ValuesOnLevelPlot2D)
BCI2 <- levelplot(array(mean.bias[['IBMA']], dim=rep(sqrt(prod(DIM)),2)),
      col.regions = gray(100:0/100), at=seq(0,ceiling(max(unlist(mean.bias))),length.out=100), main='3 level GLM',xlab='',ylab='',
              colorkey = list(space = "bottom"),
              panel=ValuesOnLevelPlot2D)
grid.arrange(BCI1,BCI2,nrow=1, top = textGrob('Standardized bias (%) of each voxel over 3000 simulations.', gp=gpar(fontsize=20,font=1)))




##
###############
### Check CI using the files in MA_stats folder
###############
##

zflame <- c()
# Load in the data
for(i in 1:nsim){
  zflameUpper <- readNIfTI(paste(DATAwd[[currentWD]],'/',i,'/SCEN_1/MA_stats/zflame1uppertstat1',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
  zflameLower <- readNIfTI(paste(DATAwd[[currentWD]],'/',i,'/SCEN_1/MA_stats/zflame1lowertstat1',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
  zflame <- c(zflame, matrix(c(zflameUpper,zflameLower), ncol=1))
}
  # Make nsim number of columns
  zflame.mat <- matrix(zflame,ncol=nsim)

mean.coverage.zflame.IBMA <-
array(NA,dim=c(prod(DIM),1))
OBJ.ZFLAME <- rep(c('UpperZ','LowerZ'), each=prod(DIM))

mean.coverage.zflame.IBMA[,1] <- indicating(UPPER = zflame.mat[which(OBJ.ZFLAME=='UpperZ'),],LOWER = zflame.mat[which(OBJ.ZFLAME=='LowerZ'),],trueVal = trueVal)
mean.coverage.zflame.IBMA
# These values are super close together, not exaclty a CI...
zflame.mat[which(OBJ.ZFLAME=='UpperZ'),3000]
zflame.mat[which(OBJ.ZFLAME=='LowerZ'),3000]




##
###############
### Check the VARCOPES of the second and third level
###############
##

VARCOPELVL3 <- c()
# Load in the data
for(i in 1:1){
  dat.tmp <- readNIfTI(paste(DATAwd[[currentWD]],'/',i,'/SCEN_1/MA_stats/varcope1',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
  VARCOPELVL3 <- c(VARCOPELVL3, matrix(dat.tmp, ncol=1))
}

# Compare this with the variance used in the meta-analysis (in the weighted variance, only saved for simulation 1)
load('/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/MAvsIBMA/1/SCEN_1/MA_stats/CI.MA.weightedVariance')
cbind(VARCOPELVL3,CI.MA.weightedVariance)
summary(VARCOPELVL3)
summary(CI.MA.weightedVariance)

  # Larger values in flame...




##
###############
### Calculate a CI around the second level cope-map
###############
##

lvl2.CI.upper.t <- lvl2.CI.lower.t <- array(NA, dim=c(prod(DIM)*nstud,30))

# This means we go into a simulation and look at the K studies.
# We can then construct a CI around each of these studies

# Start with 30 simulations
for(i in 1:30){
  for(k in 1:nstud){
    INDEX <- (k-1) * prod(DIM)
      START <- INDEX + 1
      END <- INDEX + prod(DIM)
    cope <- readNIfTI(paste(DATAwd[[currentWD]],'/',i,'/SCEN_1/study',k,'_stats/cope1.nii',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
    varcope <- readNIfTI(paste(DATAwd[[currentWD]],'/',i,'/SCEN_1/study',k,'_stats/varcope1.nii',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]

    lvl2.CI.upper.t[START:END,i] <- matrix(cope,ncol=1) + (qt(0.975, df=c((nsub/nstud)-1)) * sqrt(matrix(varcope,ncol=1)))
    lvl2.CI.lower.t[START:END,i] <- matrix(cope,ncol=1) - (qt(0.975, df=c((nsub/nstud)-1)) * sqrt(matrix(varcope,ncol=1)))}
}




##
###############
### Determine the df in FLAME
###############
##

# So FLAME1 uses a noncentral multivariate t-distribution approximation for inference at the second level.
# However, it is unclear how the DOF are calculated, the UPPER or the LOWER limit?
# So we will load in a z-map from an analysis. These z-maps are calculated through: t to p to z. So now we will go from z to p to t. And then try to calculate t to p, starting from t.

# Rounding digits
digRound <- 8
# Location:
testingStudy <- '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/MAvsIBMA/1/SCEN_1/study1_stats'
zmap <- round(readNIfTI(paste(testingStudy,'/zstat1.nii',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,],digits=digRound)
tmap <- round(readNIfTI(paste(testingStudy,'/tstat1.nii',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,],digits=digRound)


# The p-map
pmap <- round(pnorm(zmap),digits=digRound)

# Lower limit of DOF = N-1
lowerDOF <- nsub-1
# Upper limit of DOF = +Inf
upperDOF <- Inf

# Now go from t to P* (p star)
LowerPmapS <- round(pt(tmap,df=lowerDOF),digits=digRound)
UpperPmapS <- round(pt(tmap,df=upperDOF),digits=digRound)

# Now have both maps next to the pmap we got from z to p
comparisson <- data.frame(
      'True p-map' = array(pmap, dim=prod(DIM)),
      'Lower df' = array(LowerPmapS, dim=prod(DIM)),
      'Upper df' = array(UpperPmapS, dim=prod(DIM)),
      'Diff.true.low' =array(pmap, dim=prod(DIM)) -  array(LowerPmapS, dim=prod(DIM)),
      'Diff.true.low' =array(pmap, dim=prod(DIM)) -  array(UpperPmapS, dim=prod(DIM)))


# Or maybe directer: start with the t-map we get from FSL, then calculate the p to z-map transformation and compare this with zmap from FSL.
LowerTrans <- round(array(qnorm(pt(tmap, df=lowerDOF)), dim=prod(DIM)),digits=digRound)
UpperTrans <- round(array(qnorm(pt(tmap, df=upperDOF)), dim=prod(DIM)),digits=digRound)
compTran <- data.frame(
  'True z-map' = array(zmap, dim=prod(DIM)),
  'Lower df' = LowerTrans,
  'Upper df' = UpperTrans,
  'Diff.z.low' = array(zmap, dim=prod(DIM)) - LowerTrans,
  'Diff.z.upp' = array(zmap, dim=prod(DIM)) - UpperTrans)
compTran


# ==> Probably the LOWER limit


##
###############
### Performing a MA on the COPE and VARCOPE images
###############
##

# We can calculate a weighted average of the COPE values using the VARCOPES from the second level.
# Furthermore we can calculate the CI using the weighted variance, also based on the COPE and VARCOPE.

weightedAverage <- upper.weightVar <- lower.weightVar <- array(NA, dim=c(prod(DIM),nsim))
# We will go to all the cope and varcopes of each study
for(i in 1:nsim){
  if(i==c(nsim/2)) print('At 50%')
  copeStud <- array(NA, dim=c(prod(DIM),nstud))
  varcopeStud <- array(NA, dim=c(prod(DIM),nstud))
  for(j in 1:nstud){
    copeStud[,j] <- readNIfTI(paste(DATAwd[[currentWD]],'/',i,'/SCEN_1/study',j,'_stats/cope1.nii', sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
    varcopeStud[,j] <- readNIfTI(paste(DATAwd[[currentWD]],'/',i,'/SCEN_1/study',j,'_stats/varcope1.nii', sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
  }
  # First caluclate the weighted average
  weightedAverage[,i] <- rowSums(copeStud*(1/varcopeStud))/rowSums((1/varcopeStud))
  # Now the weighted variance and CI
  weightVar <- (apply(((1/varcopeStud)*(copeStud - weightedAverage[,i])^2),c(1),sum))/((nstud - 1) * apply((1/varcopeStud),1,sum))
  upper.weightVar[,i] <- matrix(weightedAverage[,i],ncol=1) + (qt(0.975,df=nstud-1) * sqrt(matrix(weightVar,ncol=1)))
  lower.weightVar[,i] <- matrix(weightedAverage[,i],ncol=1) - (qt(0.975,df=nstud-1) * sqrt(matrix(weightVar,ncol=1)))
}

# Average the weighted average over the simulations
MACV.MEANBETA <- apply(weightedAverage,1,mean,na.rm=TRUE)

# Stack both distributions next to each other
colours <- c('#238b45','#2171b5')
MEANBETAS <- data.frame('Value' = c(MA.MEANBETA, MACV.MEANBETA), 'source' = rep(c('MA', 'COPE-VARCOPE'), each=prod(DIM)))
  ggplot(MEANBETAS, aes(x=Value, group=source)) + geom_density(aes(colour=source),size=2) +
  scale_x_continuous(name="Weighted average") +
  scale_colour_manual(values = colours, name='Source', labels = c('COPE-VARCOPE -> MA', "g-variance -> MA"))



#### ------ Standardized bias ------ ####
MACV.SDBETA <- apply(weightedAverage,1,sd)
mean.bias.MACV <- matrix(((abs(MACV.MEANBETA)-trueVal)/(MACV.SDBETA))*100,ncol=1)

# Let us compare this bias with the previously calculated ones
MACV.mean.bias <- list('MACV' = mean.bias.MACV,'MA' = mean.bias.MA, 'GLM' = mean.bias.IBMA)

# Average over all voxels
CVCI.bias <- data.frame(
	'Mean' = matrix(sapply(MACV.mean.bias, FUN=function(...){apply(...,2,mean)}),ncol=1),
	'SD' = matrix(sapply(MACV.mean.bias, FUN=function(...){apply(...,2,sd)}),ncol=1),
	'CI' = factor(c('COPE-VARCOPE->MA',CIs), levels=c('COPE-VARCOPE->MA',CIs), labels=c('MACV','MA', 'GLM'))
	)

# Make an extra plot
BMACV <- levelplot(array(MACV.mean.bias[['MACV']], dim=rep(sqrt(prod(DIM)),2)),
      col.regions = gray(100:0/100), at=seq(0,ceiling(max(unlist(MACV.mean.bias))),length.out=100), main='COPE-VARCOPE -> MA',xlab='',ylab='',
              colorkey = list(space = "bottom"),
              panel=ValuesOnLevelPlot2D)
grid.arrange(BCI1,BCI2,BMACV,nrow=1, top = textGrob('Standardized bias (%) of each voxel over 3000 simulations.', gp=gpar(fontsize=20,font=1)))



#### ------ CI length ------ ####
mean.length.MACV <- array(NA,dim=c(prod(DIM),1))
mean.length.MACV[,1] <- apply(upper.weightVar - lower.weightVar,1,mean)
  # Compare this with previously calculated CI lengths
  MACV.mean.length <- list('MACV' = mean.length.MACV, 'MA' = mean.lengths$MA, 'GLM' = mean.lengths$IBMA)
  # Average over all voxels
  CVCI.length <- data.frame(
    'Mean' = matrix(sapply(MACV.mean.length, FUN=function(...){apply(...,2,mean)}),ncol=1),
    'SD' = matrix(sapply(MACV.mean.length, FUN=function(...){apply(...,2,sd)}),ncol=1),
    'CI' = factor(c('COPE-VARCOPE->MA',CIs), levels=c('COPE-VARCOPE->MA',CIs), labels=c('MACV','MA', 'GLM'))
    )
  # Make an extra plot
  LEMACV <- levelplot(array(MACV.mean.length[['MACV']], dim=rep(sqrt(prod(DIM)),2)),
        col.regions = gray(100:0/100), at=seq(0,1,by=0.05), main='COPE-VARCOPE -> MA',xlab='',ylab='',
                colorkey = list(space = "bottom"),
                panel=ValuesOnLevelPlot2D)
  grid.arrange(LCI1,LCI2,LEMACV,nrow=1, top = textGrob('Average CI length of each voxel over 3000 simulations.', gp=gpar(fontsize=20,font=1)))



#### ------ CI coverage ------ ####
mean.coverage.MACV <- array(NA,dim=c(prod(DIM),1))
mean.coverage.MACV[,1] <- indicating(UPPER = upper.weightVar,LOWER = lower.weightVar,trueVal = trueVal)
  # Compare this with previously calculated CI coverages.
  MACV.mean.CI <- list('MACV' = mean.coverage.MACV,'MA' = mean.coverages$MA, 'GLM' = mean.coverages$IBMA)
  # Average over all voxels
  CVCI.bias <- data.frame(
  	'Mean' = matrix(sapply(MACV.mean.CI, FUN=function(...){apply(...,2,mean)}),ncol=1),
  	'SD' = matrix(sapply(MACV.mean.CI, FUN=function(...){apply(...,2,sd)}),ncol=1),
  	'CI' = factor(c('COPE-VARCOPE->MA',CIs), levels=c('COPE-VARCOPE->MA',CIs), labels=c('MACV','MA', 'GLM'))
  	)
  # Make an extra plot
  CIMACV <- levelplot(array(MACV.mean.CI[['MACV']], dim=rep(sqrt(prod(DIM)),2)),
        col.regions = gray(0:100/100), at=seq(0.75,1,by=0.05), main='COPE-VARCOPE -> MA',xlab='',ylab='',
                colorkey = list(space = "bottom"),
                panel=ValuesOnLevelPlot2D)
  grid.arrange(CCI1,CCI2,CIMACV,nrow=1, top = textGrob('Average coverage of each voxel over 3000 simulations.', gp=gpar(fontsize=20,font=1)))



#### ------ QQ Plots ------ ####
# I would like to compare the MACV values with the 3 level GLM and the fixed effect MA
  # Let us do this through QQ-plots
  # Start with just one voxel and compare all values over the simulations
  MACV.ONEVOX <- sort(weightedAverage[1,])
  tmp<-AllData[which(OBJ.ID=='MA.WeightedAvg'),];MA.ONEVOX <- sort(tmp[1,]);rm(tmp)
  tmp<-AllData[which(OBJ.ID=='IBMA.COPE'),];GLM.ONEVOX <- sort(tmp[1,]);rm(tmp)

  qqplot(MA.ONEVOX,GLM.ONEVOX,main='Comparing MA versus GLM approach')
  qqplot(MACV.ONEVOX,MA.ONEVOX,main='Comparing MACV versus MA approach')
  qqplot(MACV.ONEVOX,GLM.ONEVOX,main='Comparing MACV versus GLM approach')

  qqnorm(GLM.ONEVOX, main='Normal Q-Q plot with GLM approach')
  qqnorm(MA.ONEVOX, main='Normal Q-Q plot with MA approach')
  qqnorm(MACV.ONEVOX, main='Normal Q-Q plot with MACV approach')


  #######
  # Now all values
  qqplot(x=AllData[which(OBJ.ID=='MA.WeightedAvg'),], y=AllData[which(OBJ.ID=='IBMA.COPE'),], xlab='MA approach', ylab='GLM approach', main='Q-Q plot over all voxels AND simulations.')

  #######
  # All 16 voxels, averaged over all simulations
  qqplot(x=MA.MEANBETA, y=IBMA.MEANBETA, main='Q-Q plot of the voxels between MA and GLM approach')
  qqplot(x=MA.MEANBETA, y=MACV.MEANBETA, main='Q-Q plot of the voxels between MA and MACV approach')
  qqplot(x=MACV.MEANBETA, y=IBMA.MEANBETA, main='Q-Q plot of the voxels between MACV and GLM approach')




##
###############
### Comparing maps on level 2
###############
##

# I want to compare COPE maps with hedges' g on second level. Then the varcope and variance of hedges's
# Start with COPE vs Hedges'g of the first study
copeStud <- hedgeStud <- array(NA, dim=c(prod(DIM),nsim))
varcopeStud <- stweights <- array(NA, dim=c(prod(DIM),nsim))
for(i in 1:nsim){
  if(i==c(nsim/2)) print('At 50%')
  copeStud[,i] <- readNIfTI(paste(DATAwd[[currentWD]],'/',i,'/SCEN_1/study1_stats/cope1.nii', sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
    hedgeAllStud <- AllData[which(OBJ.ID=='STHEDGE'),i]
  hedgeStud[,i] <- hedgeAllStud[c(1:prod(DIM))]

  varcopeStud[,i] <- readNIfTI(paste(DATAwd[[currentWD]],'/',i,'/SCEN_1/study1_stats/varcope1.nii', sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
    stweightsAllStud <- AllData[which(OBJ.ID=='STWEIGHTS'),i]
  stweights[,i] <- stweightsAllStud[c(1:prod(DIM))]
}

str(stweights)

# Compare COPE with hedges g
  # Start with comparing over all values
  par(mfrow=c(1,2), oma=c(0,0,2,0))
  qqplot(x=copeStud, y= hedgeStud, xlab='Level 2 in GLM - all values', ylab='Meta-analysis approach')
  # Averaged over all simulations.
  qqplot(x=apply(copeStud,1,mean,na.rm=TRUE), y=apply(hedgeStud,1,mean,na.rm=TRUE),xlab='Level 2 in GLM - averaged', ylab='Meta-analysis approach')
  title("Q-Q plot comparing 2e level COPE and hedges g", outer=TRUE)
  par(mfrow=c(1,1))

  # Histogram
  frameLvl2 <- data.frame('value' = matrix(c(hedgeStud,copeStud),ncol=1), 'source' = rep(c('g', 'COPE'), each = prod(dim(copeStud))))
  ggplot(frameLvl2, aes(x=value)) + geom_density(aes(colour=source),size=2)

# Compare VARCOPE with variance of Weighted average
  # Over all values
  par(mfrow=c(1,2), oma=c(0,0,2,0))
  qqplot(x=varcopeStud, y= c(1/stweights), xlab='Level 2 in GLM', ylab='Meta-analysis approach')
  # Averaged over all simulations
  qqplot(x=apply(varcopeStud,1,mean,na.rm=TRUE), y=c(1/apply(stweights, 1, mean, na.rm=TRUE)),
    xlab='Level 2 in GLM', ylab = 'Meta-analysis approach')
  title("Q-Q plot comparing 2e level VARCOPE and variance of g in the first study", outer=TRUE)
  par(mfrow=c(1,1))

  # Histogram
  frameLvl2 <- data.frame('value' = matrix(c(as.numeric(1/stweights),as.numeric(varcopeStud)),ncol=1), 'source' = rep(c('variance g', 'VARCOPE'), each = prod(dim(copeStud))))
  ggplot(frameLvl2, aes(x=value)) + geom_histogram(aes(fill=source)) +
  theme(legend.position='top')



##
###############
### Comparing maps on level 3
###############
##

# We will go to all the cope and varcopes of each study
copeMA <- array(NA, dim=c(prod(DIM),nsim))
varcopeMA <- array(NA, dim=c(prod(DIM),nsim))
for(i in 1:nsim){
  if(i==c(nsim/2)) print('At 50%')
  copeMA[,i] <- readNIfTI(paste(DATAwd[[currentWD]],'/',i,'/SCEN_1/MA_stats/cope1.nii', sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
  varcopeMA[,i] <- readNIfTI(paste(DATAwd[[currentWD]],'/',i,'/SCEN_1/MA_stats/varcope1.nii', sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
}

# Compare COPE with Weighted average
  # Start with comparing over all values
  par(mfrow=c(1,2), oma=c(0,0,2,0))
  qqplot(x=copeMA, y= AllData[which(OBJ.ID=='MA.WeightedAvg'),], xlab='Level 3 in GLM - all values', ylab='Meta-analysis approach')
  # Averaged over all simulations (which we already did in previous part).
  qqplot(x=apply(copeMA,1,mean,na.rm=TRUE), y=MA.MEANBETA,xlab='Level 3 in GLM - averaged over simulations', ylab = 'Meta-analysis approach')
  title("Q-Q plot comparing 3e level COPE and weighted average in the MA", outer=TRUE)

  # Histogram
  frameLvl3 <- data.frame('value' = matrix(c(AllData[which(OBJ.ID=='MA.WeightedAvg'),],copeMA),ncol=1), 'source' = rep(c('WAvg', 'COPE3'), each = prod(dim(copeMA))))
  ggplot(frameLvl3, aes(x=value)) + geom_density(aes(colour=source),size=2) +
  theme(legend.position='top')

# Compare VARCOPE with variance of Weighted average
  # Over all values
  par(mfrow=c(1,2), oma=c(0,0,2,0))
  qqplot(x=varcopeMA, y= AllData[which(OBJ.ID=='CI.MA.weightedVariance'),], xlab='Level 3 in GLM - all values', ylab='Meta-analysis approach')
  # Averaged over all simulations
  qqplot(x=apply(varcopeMA,1,mean,na.rm=TRUE), y=apply(AllData[which(OBJ.ID=='CI.MA.weightedVariance'),], 1, mean, na.rm=TRUE),
    xlab='Level 3 in GLM - averaged over simulations', ylab = 'Meta-analysis approach')
  title("Q-Q plot comparing 3e level VARCOPE and weighted variance in MA", outer=TRUE)

  # Histogram
  frameLvl3 <- data.frame('value' = matrix(c(AllData[which(OBJ.ID=='CI.MA.weightedVariance'),],varcopeMA),ncol=1), 'source' = rep(c('WVar', 'VARCOPE3'), each = prod(dim(copeMA))))
  ggplot(frameLvl3, aes(x=value)) + geom_density(aes(colour=source),size=2) +
  theme(legend.position='top')






