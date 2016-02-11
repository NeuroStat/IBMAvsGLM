####################
#### TITLE:     Plotting results of the lowLevelToMeta.R null data.
#### Contents:
####
#### Source Files:
#### First Modified: 01/02/2016
#### Notes:
#################



##
###############
### Notes
###############
##

# Blocked design for individual subjects.
# Location and noise varies over subjects (only white, temporal and spatial noise, but magnitude differs).
# Two conditions, contrast is 1 -1.
# These N subjects are pooled using simple OLS pooling.
# The resulting images are converted to Hedges' g and pooled using fixed/random effects meta-analysis.



# Different takes:
#		* Take 1: multivariate, 7 scenario's
#		* Take 2: univariate (1 voxel), 1 scenario (only white noise)
#		* Take 3: increasing sample size and amount of studies (white noise only)


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


# Directory for different takes
DATAwd <- list(
	'Take[1]' = "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/Take1",
	'Take[2]' = "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/Take2",
  'Take[3]' = "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/Take3"
	)

# Prefixes
prefix <- list(
	'Take[1]' = "WNSm",
	'Take[2]' = "UNI_",
  'Take[3]' = 'SK_'
)

# Suffix
suffix <- list(
	'Take[1]' = "weightAvg_",
	'Take[2]' = "weightVar_",
  'Take[3]' = "weightVar_"
)


NUMDATAwd <- length(DATAwd)
currentWD <- 3

setwd(wd)

# Number of scenarios
if(currentWD == 1) NumScen <- 7
if(currentWD==2) NumScen <- 1
if(currentWD == 3) NumScen <- 45

# Number of conficence intervals
CIs <- c('norm','t','weightAvg')
NumCI <- length(CIs)

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


# Load in functions from FixRan study
source('~/Dropbox/PhD/PhDWork/Meta\ Analysis/R\ Code/Studie_FixRan/FixRanStudyGit.git/Development/functions.R')


# For the third scenario, we need our data frame with the sample size and studies
# First make the data frame with the combinations of subjects and studies
is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
ss <- rep(seq(10,100,by=10),9)
k <- rep(seq(2,10),each=10)
ID <- is.wholenumber(ss/k)
OverView <- data.frame('Subjects' = ss, 'Studies' = k, 'Selection' = ID)
	OverView.Sel <- OverView[ID,-3]


##
###############
### Data Wrangling
###############
##


##############################################
# CI coverages over all voxels and simulations
nsim <- 350
if(currentWD %in% c(1,3)) DIM <- c(16,16,16)
if(currentWD==2) DIM <- c(1,1,1)
mean.coverage.norm <- uppr.coverage.norm <- lowr.coverage.norm <-
	mean.coverage.t <- uppr.coverage.t <- lowr.coverage.t <-
	mean.coverage.weightVar <- uppr.coverage.weightVar <- lowr.coverage.weightVar <-
			array(NA,dim=c(prod(DIM),nsim,NumScen))

trueVal <- 0

NumErr <- c()

# Load in R objects and calculate coverage.
for(s in 1:NumScen){
  print(s)
  for(i in 1:nsim){
    if(currentWD == 1 & s==7 & i==54)next
    # Load in upper and lower confidence bounds
    try(load(paste(DATAwd[[currentWD]],'/',i,'/','SCEN_',s,'/',prefix[[currentWD]],'CI.upper.norm_',i,sep='')),silent=TRUE)
    try(load(paste(DATAwd[[currentWD]],'/',i,'/','SCEN_',s,'/',prefix[[currentWD]],'CI.lower.norm_',i,sep='')),silent=TRUE)
    try(load(paste(DATAwd[[currentWD]],'/',i,'/','SCEN_',s,'/',prefix[[currentWD]],'CI.upper.t_',i,sep='')),silent=TRUE)
    try(load(paste(DATAwd[[currentWD]],'/',i,'/','SCEN_',s,'/',prefix[[currentWD]],'CI.lower.t_',i,sep='')),silent=TRUE)
    try(load(paste(DATAwd[[currentWD]],'/',i,'/','SCEN_',s,'/',prefix[[currentWD]],'CI.upper.',suffix[[currentWD]],i,sep='')),silent=TRUE)
    try(load(paste(DATAwd[[currentWD]],'/',i,'/','SCEN_',s,'/',prefix[[currentWD]],'CI.lower.',suffix[[currentWD]],i,sep='')),silent=TRUE)
			# Check if objects exist
			if(!exists(x=c('CI.upper.norm','CI.lower.norm','CI.lower.t','CI.lower.t','CI.lower.weightVar','CI.lower.weightVar'))) NumErr <- c(NumErr,1)
			if(!exists(x=c('CI.upper.norm','CI.lower.norm','CI.lower.t','CI.lower.t','CI.lower.weightVar','CI.lower.weightVar'))) next
		# Calculate indicators for coverages and save into vector
    mean.coverage.norm[,i,s] <- ifelse(trueVal >= CI.lower.norm & trueVal <= CI.upper.norm, 1, 0)
    mean.coverage.t[,i,s] <- ifelse(trueVal >= CI.lower.t & trueVal <= CI.upper.t, 1, 0)
		if(currentWD==1){
			mean.coverage.weightVar[,i,s] <- ifelse(trueVal >= CI.lower.weightAvg & trueVal <= CI.upper.weightAvg, 1, 0)
		}
		if(currentWD %in% c(2,3)){
			mean.coverage.weightVar[,i,s] <- ifelse(trueVal >= CI.lower.weightVar & trueVal <= CI.upper.weightVar, 1, 0)
		}
    rm(CI.upper.norm,CI.lower.norm,CI.upper.t,CI.lower.t)
			if(currentWD==1) rm(CI.upper.weightAvg,CI.lower.weightAvg)
			if(currentWD %in% c(2,3)) rm(CI.upper.weightVar,CI.lower.weightVar)
  }
}
# Put the 3 CI coverages in a list
mean.coverages <- list('norm' = mean.coverage.norm,'t' = mean.coverage.t, 'weightVar' = mean.coverage.weightVar)

head(mean.coverages)
str(mean.coverages)

  # As this takes some time, let us save this.
  save(mean.coverages, file=paste(wd,'/RObjects/',date,'-mean_coverages',sep=''))
  # Load in object
	if(currentWD == 1){			# RERUN: OK
  	load(paste(wd,'/RObjects/2016-02-09-mean_coverages',sep=''))
	}
	if(currentWD == 2){
  	load(paste(wd,'/RObjects/2016-02-08-mean_coverages',sep=''))
	}



# Make a data frame with the average coverage over all simulations and voxels with SD as well and variable for scenario
CI.coverages <- data.frame(
    'Mean' = rep(NA,c(NumScen * NumCI)),
    'SD' = rep(NA,c(NumScen * NumCI)),
    'Scenario' = rep(NA,c(NumScen * NumCI)),
    'CI' = rep(NA,c(NumScen * NumCI))
	)
	if(currentWD == 3){
		CI.coverages$SampleSize <- rep(NA,c(NumScen * NumCI))
		CI.coverages$Studies <- rep(NA,c(NumScen * NumCI))
	}

for(s in 1:NumScen){
  for(i in 1:NumCI){
    index <- (NumCI * (s-1)) + i
    data.CI <- mean.coverages[[i]]
    data.Scen <- data.CI[,,s]
      if(currentWD %in% c(1,3)) mean.Scen <- round(apply(data.Scen,1,mean,na.rm=TRUE),4)
			if(currentWD==2) mean.Scen <- round(mean(data.Scen,na.rm=TRUE),4)
        mean.vox <- round(mean(mean.Scen,na.rm=TRUE),4)
      sd <- round(sd(mean.Scen,na.rm=TRUE),4)
    # Put the values in the data frame
    CI.coverages[index,'Mean'] <- mean.vox
    CI.coverages[index,'SD'] <- sd
    CI.coverages[index,'Scenario'] <- s
    CI.coverages[index,'CI'] <- i
		CI.coverages[index,'SampleSize'] <- OverView.Sel[s,'Subjects']
		CI.coverages[index,'Studies'] <- OverView.Sel[s,'Studies']

  }
}


CI.coverages
  CI.coverages$CI <- factor(CI.coverages$CI, labels=CIs)
		if(currentWD == 3){
			CI.coverages$SampleSize <- factor(CI.coverages$SampleSize)
			CI.coverages$Studies <- factor(CI.coverages$Studies)
		}



##############################################
# Weighted averages

mean.wAvg <- array(NA,dim=c(prod(DIM),nsim,NumScen))

# Load in R objects
INDICATOR <- seq(1,nsim,by=c(nsim/10))
for(s in 1:NumScen){
  print(paste('Scenario ', s, sep=''))
  for(i in 1:nsim){
    if(s==7 & i==54)next
    if(i %in% INDICATOR)print(i)
    try(load(paste(DATAwd[[currentWD]],'/',i,'/','SCEN_',s,'/',prefix[[currentWD]],'WeightedAvg_',i,sep='')),silent=TRUE)
			if(!exists('WeightedAvg')) next
    mean.wAvg[,i,s] <- WeightedAvg
    rm(WeightedAvg)
  }
}

head(mean.wAvg)
dim(mean.wAvg)

# Average over all simulations
mean.wAvg.sim <- apply(mean.wAvg,c(1,3),mean,na.rm=TRUE)
  mean.wAvg.sim.F <- data.frame('WeightedAvg' = matrix(mean.wAvg.sim,ncol=1),
                                'Scenario' = rep(c(1:NumScen), each = prod(DIM)))
  mean.wAvg.sim.F$Scenario <- factor(mean.wAvg.sim.F$Scenario)
		if(currentWD == 3){
			mean.wAvg.sim.F$SampleSize <- factor(rep(OverView.Sel[,'Subjects'],each=prod(DIM)))
			mean.wAvg.sim.F$Studies <- factor(rep(OverView.Sel[,'Studies'],each=prod(DIM)))
		}


# One random voxel in a data frame
  VOX <- c(sample(x=c(1:prod(DIM)),size=1))
  VOX <- prod(c(5,5,5))
	if(currentWD==2) VOX <- c(1,1,1)
mean.wAvg.OneVox <- data.frame('WeightAvg' = matrix(mean.wAvg[VOX,,],ncol=1),
                               'Scenario' = rep(c(1:NumScen), each = nsim))
  mean.wAvg.OneVox$Scenario <- factor(mean.wAvg.OneVox$Scenario)


##
###############
### Plotting: coverage and weighted averages
###############
##

###---------------###
##### TAKE 1
###---------------###
limits <- aes(ymax = Mean + SD, ymin= Mean - SD)
labels <- c('White & Not smoothed',
    'White & Smoothed',
    '0.84,0.05,0.02,0.02,0.02,0.05',
    '0.64,0.15,0.02,0.02,0.02,0.15',
    '0.45,0.25,0.02,0.02,0.02,0.25',
    '0.24,0.35,0.02,0.02,0.02,0.35',
    '0.04,0.45,0.02,0.02,0.02,0.45')

# First quadrant: coverages
Q1 <- ggplot(CI.coverages, aes(x=factor(Scenario), y=Mean, colour=CI)) +
  geom_point(aes(colour=CI),size=3) +
  geom_line(aes(group=CI)) +
  scale_x_discrete(name="", labels=labels) +
  geom_errorbar(limits,width=0.15) +
  scale_colour_manual(values = c('#2b8cbe', '#016c59','#8c510a'), name='CI', labels = c('Normal', 't', 'weighted variance')) +
  ggtitle(label='Coverages of 3 types of CI over all voxels and simulations') +
  theme(plot.title = element_text(lineheight=.4,size=10, face="plain"),
    axis.text.x = element_text(angle = 315, hjust = 0, vjust = 0.95))



# Second quadrant: distribution of weighted average from random voxel over all simulations
colours <- sort(c('#00B6EB', '#00C094', '#53B400', '#A58AFF', '#C49A00', '#F8766D', '#FB61D7'))
Q2 <- ggplot(mean.wAvg.OneVox, aes(x=WeightAvg,group = Scenario)) + geom_density(aes(fill = Scenario),adjust=1,position='stack', colour='white') +
scale_x_continuous(name="") +
geom_vline(xintercept=0,colour='red') +
scale_fill_manual(values= colours, name="Scenario", labels = labels) +
ggtitle(label='Weighted averages over all simulations in ONE voxel') + theme_bw() +
  theme(plot.title = element_text(lineheight=.4,size=10, face="plain"),
    legend.position="left")



# Third quadrant: distribtution of weighted average over all voxels and simulations
Q3 <- ggplot(mean.wAvg.sim.F, aes(x=WeightedAvg,group = Scenario)) + geom_density(aes(fill = Scenario),adjust=1,position='stack', colour='white') +
scale_x_continuous(name="") +
geom_vline(xintercept=0,colour='red') +
scale_fill_manual(values= colours, name="Scenario", labels = labels) +
scale_colour_manual(values= 'white', name="Scenario", labels = labels) +
ggtitle(label='Weighted averages over all simulations in ALL voxels') + theme_bw() +
  theme(plot.title = element_text(lineheight=.4,size=10, face="plain"),
    legend.position="none")

# Fourth quadrant: levelplot of third scenario with the weighted averages
Q4 <- levelplot(array(mean.wAvg.sim[,3],dim=DIM),col.regions = heat.colors,main='Levelplot of weighted averages from third scenario',
    xlab='',ylab='')


# Arrange in grids
grid.arrange(Q2,Q1,Q3,layout_matrix = cbind(c(1,3), c(2,2)))


###---------------###
##### TAKE 2
###---------------###

labels <- c('White & Not smoothed')

# First quadrant: coverages
Q1 <- ggplot(CI.coverages, aes(x=factor(Scenario), y=Mean, colour=CI)) +
  geom_point(aes(colour=CI),size=3) +
  scale_x_discrete(name="", labels=labels) +
  scale_colour_manual(values = c('#2b8cbe', '#016c59','#8c510a'), name='CI', labels = c('Normal', 't', 'weighted variance')) +
  ggtitle(label='Coverages of 3 types of CI in univariate setting') +
  theme(plot.title = element_text(lineheight=.4,size=10, face="plain"))


# Second quadrant: distribution of weighted average from the one voxel over all simulations
colours <- c('#00B6EB')
Q2 <- ggplot(mean.wAvg.OneVox, aes(x=WeightAvg,group = Scenario)) + geom_density(aes(fill = Scenario),adjust=1,position='stack', colour='white') +
scale_x_continuous(name="") +
geom_vline(xintercept=0,colour='red') +
scale_fill_manual(values= colours, name="Scenario", labels = labels) +
ggtitle(label='Weighted averages over all simulations') + theme_bw() +
  theme(plot.title = element_text(lineheight=.4,size=10, face="plain"),
    legend.position="left")


# Arrange in grids
grid.arrange(Q2,Q1,nrow=1)




###---------------###
##### TAKE 3
###---------------###
limits <- aes(ymax = Mean + SD, ymin= Mean - SD)
colours <- c('#d53e4f','#f46d43','#fdae61','#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd')
ggplot(CI.coverages, aes(x=SampleSize, y=Mean, colour=Studies)) +
  geom_point(aes(colour=Studies),size=1.3) +
  geom_line(aes(group=Studies), size=1) +
	facet_wrap( ~ CI) +
  scale_x_discrete(name="Sample Size") +
  scale_colour_manual(values = colours, name='Amount of studies', labels = c(2:10)) +
  ggtitle(label='Coverages of 3 types of CI over all voxels and simulations') +
  theme(plot.title = element_text(lineheight=.4,size=13, face="plain"),
    axis.text.x = element_text(angle = 315, hjust = 0, vjust = 0.95),
		legend.key = element_rect(fill='#d9d9d9', colour = '#d9d9d9'),
		legend.background = element_rect(colour = '#d9d9d9', fill = '#d9d9d9'),
		strip.background = element_rect(fill='#d9d9d9'))


# Maybe we should go to 3D scatterplot
	# Surface plot based on linear regression of the sample size and amount of studies on the mean coverage of each coverage
# We create a X and Z value, which are 10 values
x3 <- seq(1,10,1)
y3 <- x3
	# We will use a function that creates for each interesection of x3 an y3 a value, based on the coefficients of the linear regression
	FU <- function(x,y,coef) as.numeric(coef[1] + (x * coef[2]) + (y * coef[3]) + (x * y * coef[4]))
# Here we have al the linear regressions
SUR.lm.norm <- lm(Mean ~ as.numeric(SampleSize) * as.numeric(Studies), data=CI.coverages[CI.coverages$CI=='norm',])
SUR.lm.t <- lm(Mean ~ as.numeric(SampleSize) * as.numeric(Studies), data=CI.coverages[CI.coverages$CI=='t',])
SUR.lm.weightVar <- lm(Mean ~ as.numeric(SampleSize) * as.numeric(Studies), data=CI.coverages[CI.coverages$CI=='weightAvg',])
	# The coefficients from the three models
	coef.lm.norm <- SUR.lm.norm$coef
	coef.lm.t <- SUR.lm.t$coef
	coef.lm.weightVar <- SUR.lm.weightVar$coef
# Now we create the Z values according to the combinations of x3 and y3 (outer function)
zval.norm <- outer(x3,y3,FU, coef=coef.lm.norm)
zval.t <- outer(x3,y3,FU, coef=coef.lm.t)
zval.weightVar <- outer(x3,y3,FU, coef=coef.lm.weightVar)
# Now go to plotting
colours <- c('#1b9e77','#d95f02','#7570b3')
persp(x3,y3,zval.weightVar,theta=125,
	phi=10,xlab='Sample Size (x10)',ylab='Studies',zlab='Emperical Coverage',
	zlim=c(0,1), xlim=c(1,10),ylim=c(1,10), ticktype = 'detailed',main='', border=colours[1],axes=TRUE,box=TRUE,d=2,lwd=2,nticks=10)
	par(new = TRUE)
persp(x3,y3,zval.norm,theta=125,
     phi=10,xlab='',ylab='',zlab='',
		 zlim=c(0,1), xlim=c(1,10),ylim=c(1,10), ticktype = 'detailed',main='', border=colours[2],axes=FALSE,box=FALSE,d=2,lwd=2)
par(new = TRUE)
persp(x3,y3,zval.t,theta=125,
    phi=10,xlab='',ylab='',zlab='',
		 zlim=c(0,1), xlim=c(1,10),ylim=c(1,10), ticktype = 'detailed',main='Average emperical coverage over all voxels and simulations.', border=colours[3],box=FALSE,d=2,lwd=2)
legend(x='topright', legend=CIs[c(3,1,2)], bty='n', text.col = colours,cex=1.5)


# Third quadrant: distribtution of weighted average over all voxels and simulations
ggplot(mean.wAvg.sim.F, aes(x=WeightedAvg,group = Scenario)) + geom_density(aes(fill = Scenario),adjust=1,position='stack', colour='white') +
scale_x_continuous(name="") +
geom_vline(xintercept=0,colour='red') +
scale_fill_manual(values= colours, name="Scenario", labels = labels) +
scale_colour_manual(values= 'white', name="Scenario", labels = labels) +
ggtitle(label='Weighted averages over all simulations in ALL voxels') + theme_bw() +
  theme(plot.title = element_text(lineheight=.4,size=10, face="plain"),
    legend.position="none")

head(mean.wAvg.sim.F)
summary(mean.wAvg.sim.F[mean.wAvg.sim.F$Scenario==45,])

hist(mean.wAvg.sim.F[mean.wAvg.sim.F$SampleSize==10,'WeightedAvg'])
par(new=TRUE)
hist(mean.wAvg.sim.F[mean.wAvg.sim.F$SampleSize==100,'WeightedAvg'],axes=FALSE,xlab='',ylab='',col='red')

plot(density(mean.wAvg.sim.F[mean.wAvg.sim.F$SampleSize==10,'WeightedAvg']),axes=FALSE,xlab='',ylab='',border='red')
plot(density(mean.wAvg.sim.F[mean.wAvg.sim.F$SampleSize==100,'WeightedAvg']),axes=FALSE,xlab='',ylab='',border='blue')

# Compare MPG distributions for cars with
# 4,6, or 8 cylinders
library(sm)
attach(mtcars)

# create value labels
cyl.f <- factor(cyl, levels= c(4,6,8),
  labels = c("4 cylinder", "6 cylinder", "8 cylinder"))

# plot densities
sm.density.compare(mpg, cyl, xlab="Miles Per Gallon")
title(main="MPG Distribution by Car Cylinders")

# add legend via mouse click
colfill<-c(2:(2+length(levels(cyl.f))))
legend(locator(1), levels(cyl.f), fill=colfill)



summary(mean.wAvg.sim.F[mean.wAvg.sim.F$Scenario==3,])
head(mean.wAvg.sim.F)

##
###############
### Plotting: simulations needed
###############
##


###---------------###
##### TAKE 1
###---------------###

# Start with 500 simulations.
  # For each 10 simulations, we calculate the SD.
  # Then we add 10 and calculate the SD of the

# Number of calculations
Calcs <- seq(600,nsim,by=10)
  PreCalcs <- seq(501,(nsim-9),by=10)
NumCalc <- length(Calcs)

# Temporarly change labeling
CIs <- c('norm','t','weightVar')

# MEANS's
MEANvox <- data.frame(
    'Mean' = array(NA,dim=c(NumScen*NumCI*NumCalc)),
    'CI' = array(NA,dim=c(NumScen*NumCI*NumCalc)),
    'Scenario' = array(NA,dim=c(NumScen*NumCI*NumCalc))
    )

for(s in 1:NumScen){
  print(paste('At scenario ', s, sep=''))
  for(j in 1:NumCI){
   # Pre-LOOP
    DataMEAN <- mean.coverages[[CIs[j]]]
     dat.tmp <- DataMEAN[,c(1:500),s]
     dat.mean <- apply(dat.tmp,1,mean,na.rm=TRUE)
      rm(dat.tmp)
    for(i in 1:NumCalc){
      Index <- ((s-1) * (NumCalc*NumCI)) + ((j-1) * NumCalc) + i
      dat.tmp <- DataMEAN[,c(PreCalcs[i]:Calcs[i]),s]
      dat.mean <- cbind(dat.mean,apply(dat.tmp,1,mean,na.rm=TRUE))

      # Calculate MEAN
      MEANsim.tmp <- apply(dat.mean,1,mean,na.rm=TRUE)
      MEANvox[Index,'Mean'] <- mean(MEANsim.tmp,na.rm=TRUE)
      MEANvox[Index, 'CI'] <- j
      MEANvox[Index, 'Scenario'] <- s

      # Remove objects
      rm(MEANsim.tmp,dat.tmp)
    }
    rm(dat.mean)
  }
}

MEANvox$Time <- rep(c(1:NumCalc),times=c(NumCI*NumScen))
  MEANvox$Scenario <- factor(MEANvox$Scenario)
  MEANvox$CI <- factor(MEANvox$CI, labels = CIs)
head(MEANvox)

# As this takes some time, let us save this.
save(MEANvox, file=paste(wd,'/RObjects/',date,'-MEANvox',sep=''))
# Load in object
load(paste(wd,'/RObjects/2016-02-03-MEANvox',sep=''))
	  MEANvox$CI <- factor(MEANvox$CI, labels = CIs)

# Plot
LabelTime <- seq(600,nsim,length.out=6)
ggplot(MEANvox, aes(x = Time, y = Mean, group = CI)) +
  geom_line(aes(colour = CI)) +
  facet_wrap( ~ Scenario)+
	scale_x_continuous(name='Number of simulations', labels = LabelTime) +
  scale_y_continuous(name='Average emperical coverage')



















