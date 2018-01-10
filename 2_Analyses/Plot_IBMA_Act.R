####################
#### TITLE:     Plotting results of the MAvsIBMA_Act.R data.
#### Contents:
####
#### Source Files:
#### First Modified: 24/09/2017
#### Notes:
#################



##
###############
### Notes
###############
##

# Activation data!

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
  'Take[MAvsIBMA_Act]' = "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/MAvsIBMA_act/Results"
)
NUMDATAwd <- length(DATAwd)
currentWD <- 1

# Number of conficence intervals
CIs <- c('MA-weightVar','GLM-t')
NumCI <- length(CIs)

# Number of executed simulations
nsim.tmp <- matrix(c(
  1,1000
), ncol=2, byrow=TRUE)
nsim <- nsim.tmp[currentWD,2]


# Number of subjects and studies
nsub <- 100
nstud <- 5

# Dimension of brain
DIM.tmp <- array(NA, dim=c(NUMDATAwd,3))
DIM.tmp[c(1),] <- c(9,9,9)
DIM <- DIM.tmp[currentWD,]


# Load in libraries
library(AnalyzeFMRI)
library(lattice)
library(grid)
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

################
#### TRUE VALUES: NEED TO RECREATE DESIGN MATRIX
################
## BETA 1 ##
# Signal characteristics
TR <- 2
nscan <- 200
total <- TR*nscan
on1 <- seq(1,total,40)
onsets <- list(on1)
duration <- list(20)
whiteSigma <- 7

# Generating a design matrix
X <- simprepTemporal(total,1,onsets=onsets,effectsize = 100, durations=duration,
                     TR = TR, acc=0.1, hrf="double-gamma")

# Generate design
pred <- simTSfmri(design=X, base=100, SNR=1, noise="none", verbose=FALSE)

# We know that: range(signal) = Beta_1 * range(X), so we solve for Beta_1 to get its true value.
es1 <- 0.01
base <- 100
signal <- es1 * (pred-base) + base
trueBeta1 <- (max(signal) - min(signal)) / (max(pred) - min(pred))

## VAR(BETA1) ##
# Variance of beta 1 = variance of beta 1 at subject level / N.
# Variance of beta1_subj = whiteSigma^2/SSX
sigma2_Beta1_subj <- whiteSigma^2 / (sum((pred - mean(pred))^2))
sigma2_Beta1 <- sigma2_Beta1_subj / nsub
## TVALUE ##
trueTval <- trueBeta1 / sqrt(sigma2_Beta1)

## HEDGES' G ##
trueG <- hedgeG(t = trueTval, N = nsub)

## LOCATION ##
# True (unsmoothed) values (object called GroundTruth, created in MAvsIBMA_Act.R)
load('/Users/hanbossier/Dropbox/PhD/PhDWork/Meta Analysis/R Code/Studie_Simulation/SimulationGit/2_Analyses/GroundTruth_Act.Rda')

# True values: either for Beta 1 or hedges' g
GroundTruthCope <- GroundTruth
GroundTruthCope[GroundTruth == 1] <- trueBeta1
GroundTruthCope[GroundTruth == 0] <- NA
GroundTruthES <- GroundTruth
GroundTruthES[GroundTruth == 1] <- trueG
GroundTruthES[GroundTruth == 0] <- NA
TrueLocations <- c(4,4,5)


cop <- readNIfTI(paste(DataWrite,"/study",t,"_stats/cope1.nii",sep=""), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
varc <- readNIfTI(paste(DataWrite,"/study",t,"_stats/varcope1.nii",sep=""), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]

dim(cop)
array(apply(COPE, 1, mean), dim = DIM)[4,4,5]
cop[4,4,5]
array(apply(COPE, 1, t.test), dim = DIM)[4,4,5]
array(apply(COPE, 1, var), dim = DIM)[4,4,5]

varc[4,4,5]
2.099544e-06 / nsub

array(apply(COPE, 1, mean), dim = DIM)[4,4,5] /
  sqrt(2.099544e-06 / nsub)
STMAP[4,4,5]


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


# Temp save
#save(AllData, file = '/Volumes/1_5_TB_Han_HDD/Temp/AllData.Rda')
#save(objects, file = '/Volumes/1_5_TB_Han_HDD/Temp/objects.Rda')

#########################################################
###################### CI COVERAGE ######################
#########################################################

# Calculate coverage
mean.coverage.weightVar.MA <- data.frame(UPPER = c(AllData[which(OBJ.ID=='CI.MA.upper.weightVar'),]),
          LOWER = c(AllData[which(OBJ.ID=='CI.MA.lower.weightVar'),]),
          TrueValue = rep(array(GroundTruthES, dim = prod(DIM)), nsim),
          sim = rep(seq(1,nsim), each = prod(DIM)),
          voxel = rep(seq(1, prod(DIM)), nsim)) %>%
  mutate(COVERAGE = ifelse(TrueValue <= UPPER & TrueValue >= LOWER, 1, 0)) %>%
  group_by(voxel) %>% summarise(mean.coverage.weightVar.MA = mean(COVERAGE, na.rm = TRUE)) %>% tbl_df()

mean.coverage.t.IBMA <- data.frame(UPPER = c(AllData[which(OBJ.ID=='CI.IBMA.upper.t'),]),
                                         LOWER = c(AllData[which(OBJ.ID=='CI.IBMA.lower.t'),]),
                                         TrueValue = rep(array(GroundTruthCope, dim = prod(DIM)), nsim),
                                         sim = rep(seq(1,nsim), each = prod(DIM)),
                                         voxel = rep(seq(1, prod(DIM)), nsim)) %>%
  mutate(COVERAGE = ifelse(TrueValue <= UPPER & TrueValue >= LOWER, 1, 0)) %>%
  group_by(voxel) %>% summarise(mean.coverage.weightVar.IBMA = mean(COVERAGE, na.rm = TRUE)) %>% tbl_df()


data.frame(UPPER = c(AllData[which(OBJ.ID=='CI.MA.upper.weightVar'),]),
           LOWER = c(AllData[which(OBJ.ID=='CI.MA.lower.weightVar'),]),
           TrueValue = rep(array(GroundTruthES, dim = prod(DIM)), nsim)) %>% tbl_df() %>% View()


summarise(mean.coverage.t.IBMA, mean(mean.coverage.weightVar.IBMA, na.rm = TRUE))

dim(AllData[which(OBJ.ID == 'IBMA.COPE'),50])
levelplot(array(AllData[which(OBJ.ID == 'IBMA.COPE'),50], dim = DIM))

select(mean.coverage.t.IBMA, mean.coverage.weightVar.IBMA) %>% unlist(.) %>% as.numeric(.) %>%
  array(., dim = DIM)

levelplot(array(unlist(mean.coverages[['MA']]), dim = DIM)[(TrueLocations[1] - 2):(TrueLocations[1] + 2),
                                                           (TrueLocations[2] - 2):(TrueLocations[2] + 2),
                                                           (TrueLocations[3] - 2):(TrueLocations[3] + 2)])

levelplot(array(unlist(mean.coverages[['IBMA']]), dim = DIM)[(TrueLocations[1] - 2):(TrueLocations[1] + 2),
                                                           (TrueLocations[2] - 2):(TrueLocations[2] + 2),
                                                           (TrueLocations[3] - 2):(TrueLocations[3] + 2)])

levelplot(array(AllData[which(OBJ.ID == 'IBMA.COPE'),1], dim = DIM)[(TrueLocations[1] - 2):(TrueLocations[1] + 2),
                                                           (TrueLocations[2] - 2):(TrueLocations[2] + 2),
                                                           (TrueLocations[3] - 2):(TrueLocations[3] + 2)])

levelplot(array(AllData[which(OBJ.ID == 'MA.WeightedAvg'),1], dim = DIM)[(TrueLocations[1] - 2):(TrueLocations[1] + 2),
                                                                  (TrueLocations[2] - 2):(TrueLocations[2] + 2),
                                                                  (TrueLocations[3] - 2):(TrueLocations[3] + 2)])


# Check hedge g
allG <- array(AllData[which(OBJ.ID == 'STHEDGE'),], dim = c(DIM, nstud, nsim))
avgG <- NULL
for(j in 1:nsim){
  for(i in 1:nstud){
    tmp <- allG[,,,i,j]
    tmp[GroundTruth == 0] <- NA
    avgG <- c(avgG, mean(tmp, na.rm = TRUE))
  }
}
mean(avgG)
trueG     # Hmmm



hist(AllData[which(OBJ.ID == 'IBMA.COPE'),])
hist(AllData[which(OBJ.ID == 'MA.WeightedAvg'),])

# Put the 2 coverages in a list
mean.coverages <- list('MA' = mean.coverage.weightVar.MA[,2],'IBMA' = mean.coverage.t.IBMA[,2])

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

mean.bias.MA <- matrix(((abs(MA.MEANBETA)-matrix(GroundTruth, ncol = 1))/(MA.SDBETA))*100,ncol=1)
mean.bias.IBMA <- matrix(((abs(IBMA.MEANBETA)-matrix(GroundTruth, ncol = 1))/(IBMA.SDBETA))*100,ncol=1)

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

# BETTER PLOT 3D
levelplot()

# Function for levelplot
ValuesOnLevelPlot2D <- function(x, y, z, ...) {
  panel.levelplot(x,y,z,...)
  panel.text(x, y, round(z,3),col='red')
}

# Plotting the coverages
levelplot(array(unlist(mean.coverages[['MA']]), dim = DIM))
levelplot(array(unlist(mean.coverages[['IBMA']]), dim=rep(sqrt(prod(DIM)),2)))


CCI1 <- levelplot(array(unlist(mean.coverages[['MA']]), dim= DIM),
                  col.regions = gray(0:100/100), at=seq(0,1,by=0.05), main='Meta-Analysis',xlab='',ylab='',
                  colorkey = list(space = "bottom"))
CCI2 <- levelplot(array(unlist(mean.coverages[['IBMA']]), dim= DIM),
                  col.regions = gray(0:100/100), at=seq(0,1,by=0.05), main='3 level GLM',xlab='',ylab='',
                  colorkey = list(space = "bottom"))
grid.arrange(CCI1,CCI2,nrow=1,top = textGrob('CI - Coverage of each voxel over 1000 simulations.', gp=gpar(fontsize=20,font=1)))


# Try with ggplot
plotCI <- data.frame(CI = rbind(mean.coverages[['MA']],mean.coverages[['IBMA']]),
                     Type = rep(c("MA", "mixed effects GLM"), each = 64),
                     x = factor(rep(rep(1:8,8),2)),
                     y = factor(rep(rep(1:8, each = 8),2)))
ggplot(plotCI, aes(x=x,y=y)) + geom_tile(aes(fill = CI)) +
  facet_wrap(~Type) +
  scale_fill_gradient2(name = "", low = scales::muted("blue"), mid = "white",
                       high = scales::muted("red"), midpoint = 0.95, space = "Lab",
                       na.value = "grey50", guide = "colourbar") +
  scale_x_discrete(name = "") + scale_y_discrete(name = "") + ggtitle("") +
  theme_minimal() +
  theme(strip.text = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.9,"cm"))



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

# Try with ggplot
PlotBias <- data.frame(bias = rbind(mean.bias[['MA']],mean.bias[['IBMA']]),
                       Type = rep(c("MA", "mixed effects GLM"), each = 64),
                       x = factor(rep(rep(1:8,8),2)),
                       y = factor(rep(rep(1:8, each = 8),2)))
ggplot(PlotBias, aes(x=x,y=y)) + geom_tile(aes(fill = bias)) +
  facet_wrap(~Type) +
  scale_fill_gradient2(name = "", low = scales::muted("blue"), mid = "white",
                       high = scales::muted("red"), midpoint = 0, space = "Lab",
                       na.value = "grey50", guide = "colourbar") +
  scale_x_discrete(name = "") + scale_y_discrete(name = "") + ggtitle("Average standardized bias") +
  theme_minimal()

# Maybe try the difference between the two, alongside table of average and SD
DiffBias <- data.frame(Diffbias = c(mean.bias[['MA']] - mean.bias[['IBMA']]),
                       x = factor(rep(1:8,8)),
                       y = factor(rep(1:8, each = 8)))
DiffBiasPlot <- ggplot(DiffBias, aes(x=x,y=y)) + geom_tile(aes(fill = Diffbias)) +
  scale_fill_gradient2(name = "", low = scales::muted("blue"), mid = "white",
                       high = scales::muted("red"), midpoint = 0, space = "Lab",
                       na.value = "grey50", guide = "colourbar") +
  scale_x_discrete(name = "") + scale_y_discrete(name = "") + ggtitle("") +       # Standardized bias (MA) - standardized bias (GLM)
  theme_minimal() +
  theme(strip.text = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.9,"cm"))


tableBias <- data.frame(Type = CI.bias$CI, average = CI.bias$Mean, SD = CI.bias$SD)
tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
tbl <- tableGrob(tableBias, rows=NULL, theme=tt)
# Plot chart and table into one object
grid.arrange(tbl , DiffBiasPlot,
             nrow=1,
             as.table=TRUE,
             widths=c(1,2))





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






