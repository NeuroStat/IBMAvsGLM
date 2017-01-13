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
#		* Take 4: increasing sample size and amount of studies (white noise only) for UNIVARIATE approach
#		* Take 5: Same as take 3, but with R file differently sourced.
#   * Take 6: Simple design (only one condition), grid of 2 x 2 x 2 voxels (white noise only)


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
	'Take[1]' = "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/Take1",
	'Take[2]' = "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/Take2",
  'Take[3]' = "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/Take3",
	'Take[4]' = "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/Take4/Results",
	'Take[5]' = "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/Take5",
  'Take[6]' = "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/Take6",
	'Take[7]' = "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/IBMA",
	'Take[8]' = "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/Take8"
	)

# Prefixes
prefix <- list(
	'Take[1]' = "WNSm",
	'Take[2]' = "UNI_",
  'Take[3]' = 'SK_',
	'Take[4]' = "UNI_SK_",
	'Take[5]' = 'SK_',
  'Take[6]' = 'SD_',
	'Take[7]' = 'IBMA_',
	'Take[8]' = 'SDG_'
)

NUMDATAwd <- length(DATAwd)
currentWD <- 8

# Number of scenarios
NumScen.tmp <- matrix(c(
                1,7,
                2,1,
                3,45,
                4,45,
                5,45,
                6,30,
								7,30,
								8,30
              ), ncol=2, byrow=TRUE)
NumScen <- NumScen.tmp[currentWD,2]


# Number of conficence intervals
CIs <- c('norm','t','weightVar')
NumCI <- length(CIs)


# Number of executed simulations
nsim.tmp <- matrix(c(
                1,4176,
                2,3000,
                3,350,
                4,1500,
                5,500,
                6,3000,
								7,3000,
								8,3000
              ), ncol=2, byrow=TRUE)
nsim <- nsim.tmp[currentWD,2]


# Dimension of brain
DIM.tmp <- array(NA, dim=c(NUMDATAwd,3))
	DIM.tmp[c(1,3,5,8),] <- c(16,16,16)
	DIM.tmp[c(2,4),] <- c(1,1,1)
	DIM.tmp[c(6,7),] <- c(2,2,2)
DIM <- DIM.tmp[currentWD,]

# Number of subjects and studies
TablesOverview <- list(
	'[1]' = '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/Take3/OverView.Sel.txt',
	'[2]' = '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/Take6/OverView.txt'
	)
overview.tmp <- matrix(c(			# This takes the element from TablesOverview
                1,NA,
                2,NA,
                3,1,
                4,1,
                5,1,
                6,2,
								7,2,
								8,2
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

###############
# Data has been pre-processed in the PreProcess.R file. All simulations are there being read in and combined.
	# Structue: each scenario = a list
		# Each list: rows = voxels, columns = simulations
# Load in these objects
OBJECTS <- c('CI.upper.norm.All','CI.lower.norm.All','CI.upper.t.All',
	'CI.lower.t.All','CI.upper.weightVar.All','CI.lower.weightVar.All',
	'WeightedAvg.All')
	NumObjects <- length(OBJECTS)
for(i in 1:NumObjects){
	load(paste(DATAwd[[currentWD]], '/Take',currentWD,'-',OBJECTS[i],sep=''))
	print(paste('Loaded object ', OBJECTS[i], sep=''))
}

# True value
trueVal <- 0

# CI coverages over all voxels and simulations
mean.coverage.norm <-
mean.coverage.t <-
mean.coverage.weightVar <-
mean.wAvg <-
		array(NA,dim=c(prod(DIM),NumScen))

# For loop over all scenarios
for(s in 1:NumScen){
	print(s)
	mean.coverage.norm[,s] <- indicating(UPPER = CI.upper.norm.All[[s]], LOWER = CI.lower.norm.All[[s]], trueVal = trueVal)
	mean.coverage.t[,s] <- indicating(UPPER = CI.upper.t.All[[s]], LOWER = CI.lower.t.All[[s]], trueVal = trueVal)
	mean.coverage.weightVar[,s] <- indicating(UPPER = CI.upper.weightVar.All[[s]], LOWER = CI.lower.weightVar.All[[s]], trueVal = trueVal)
}

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
	if(!currentWD %in% c(1,2)){
		CI.coverages$SampleSize <- factor(CI.coverages$SampleSize)
		CI.coverages$Studies <- factor(CI.coverages$Studies)
	}



WeightedAvg.All



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
		if(currentWD %in% c(3,4)){
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
	geom_hline(yintercept=0.95,colour='red') +
	facet_wrap( ~ CI) +
  scale_x_discrete(name="Sample Size") +
  scale_colour_manual(values = colours, name='Amount of studies', labels = c(2:10)) +
  ggtitle(label='Coverages of 3 types of CI over all voxels and simulations') +
  theme(plot.title = element_text(lineheight=.4,size=13, face="plain"),
    axis.text.x = element_text(angle = 315, hjust = 0, vjust = 0.95),
		legend.key = element_rect(fill='#d9d9d9', colour = '#d9d9d9'),
		legend.background = element_rect(colour = '#d9d9d9', fill = '#d9d9d9'),
		strip.background = element_rect(fill='#d9d9d9'),
		legend.position="top")


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
SUR.lm.weightVar <- lm(Mean ~ as.numeric(SampleSize) * as.numeric(Studies), data=CI.coverages[CI.coverages$CI=='weightVar',])
	# The coefficients from the three models
	coef.lm.norm <- SUR.lm.norm$coef
	coef.lm.t <- SUR.lm.t$coef
	coef.lm.weightVar <- SUR.lm.weightVar$coef
# Now we create the Z values according to the combinations of x3 and y3 (outer function)
zval.norm <- outer(x3,y3,FU, coef=coef.lm.norm)
zval.t <- outer(x3,y3,FU, coef=coef.lm.t); zval.t[zval.t > 1] <- 1
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



### Investigate one voxel
coverage.norm <- mean.coverages$norm
mean.coverages$norm
str(mean.coverages)
str(mean.coverages$norm)
str(mean.coverage.norm)
mean.coverage.norm[300,,3]
mean(mean.coverage.norm[300,,3],na.rm=TRUE)


###---------------###
##### TAKE 4
###---------------###

# Let's start with a selection of CI.coverages (2, 5 or 10 studies)
IDStud <- CI.coverages$Studies %in% c(2,5,10)
CI.coverages.Sel <- CI.coverages[IDStud,]
colours <- c('#1b9e77','#d95f02','#7570b3')
ggplot(CI.coverages.Sel, aes(x=SampleSize, y=Mean, colour=Studies)) +
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



# 3D scatterplot
	# Surface plot based on linear regression of the sample size and amount of studies on the mean coverage of each coverage
# We create a X and Z value, which are 10 values
x3 <- seq(1,10,1)
y3 <- seq(2,10,1)
	# We will use a function that creates for each interesection of x3 an y3 a value, based on the coefficients of the linear regression
	FU <- function(x,y,coef) as.numeric(coef[1] + (x * coef[2]) + (y * coef[3]) + (x * y * coef[4]))
# Here we have al the linear regressions
SUR.lm.norm <- lm(Mean ~ as.numeric(SampleSize) * as.numeric(Studies), data=CI.coverages[CI.coverages$CI=='norm',])
SUR.lm.t <- lm(Mean ~ as.numeric(SampleSize) * as.numeric(Studies), data=CI.coverages[CI.coverages$CI=='t',])
SUR.lm.weightVar <- lm(Mean ~ as.numeric(SampleSize) * as.numeric(Studies), data=CI.coverages[CI.coverages$CI=='weightVar',])
	# The coefficients from the three models
	coef.lm.norm <- SUR.lm.norm$coef
	coef.lm.t <- SUR.lm.t$coef
	coef.lm.weightVar <- SUR.lm.weightVar$coef
# Now we create the Z values according to the combinations of x3 and y3 (outer function)
zval.norm <- outer(x3,y3,FU, coef=coef.lm.norm)
zval.t <- outer(x3,y3,FU, coef=coef.lm.t); zval.t[zval.t > 1] <- 1
zval.weightVar <- outer(x3,y3,FU, coef=coef.lm.weightVar)
NOMINAL <- matrix(0.95,ncol=9,nrow=10)
# Now go to plotting
colours <- c('#1b9e77','#d95f02','#7570b3')
persp(x3,y3,NOMINAL,theta=134,
	phi=10,xlab='Sample Size (x10)',ylab='Studies',zlab='Emperical Coverage',
	zlim=c(0.9,1), xlim=c(1,10),ylim=c(2,10), ticktype = 'detailed',main='', border='black',axes=TRUE,box=TRUE,d=2,lwd=2,nticks=10)
	par(new = TRUE)
persp(x3,y3,zval.weightVar,theta=134,
	phi=10,xlab='',ylab='',zlab='',
	zlim=c(0.9,1), xlim=c(1,10),ylim=c(2,10), ticktype = 'detailed',main='', border=colours[1],axes=FALSE,box=FALSE,d=2,lwd=2)
	par(new = TRUE)
persp(x3,y3,zval.norm,theta=134,
     phi=10,xlab='',ylab='',zlab='',
		 zlim=c(0.9,1), xlim=c(1,10),ylim=c(2,10), ticktype = 'detailed',main='', border=colours[2],axes=FALSE,box=FALSE,d=2,lwd=2)
par(new = TRUE)
persp(x3,y3,zval.t,theta=134,
    phi=10,xlab='',ylab='',zlab='',
		 zlim=c(0.9,1), xlim=c(1,10),ylim=c(2,10), ticktype = 'detailed',main='Average emperical coverage over ONE voxels and simulations.', border=colours[3],box=FALSE,d=2,lwd=2)
legend(x=0.26,y=-0.1, legend=c(CIs[c(3,1,2)],'nominal'), bty='n', text.col = c(colours,'black'),cex=1.5)




###---------------###
##### TAKE 6
###---------------###

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








##
###############
### Investigating multiple files in Take 3
###############
##


nsim <- 350
if(currentWD %in% c(1,3)) DIM <- c(16,16,16)

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
    # Load in upper and lower confidence bounds
    try(load(paste(DATAwd[[currentWD]],'/',i,'/','SCEN_',s,'/',prefix[[currentWD]],'CI.upper.norm_',i,sep='')),silent=TRUE)
    try(load(paste(DATAwd[[currentWD]],'/',i,'/','SCEN_',s,'/',prefix[[currentWD]],'CI.lower.norm_',i,sep='')),silent=TRUE)
    try(load(paste(DATAwd[[currentWD]],'/',i,'/','SCEN_',s,'/',prefix[[currentWD]],'CI.upper.t_',i,sep='')),silent=TRUE)
    try(load(paste(DATAwd[[currentWD]],'/',i,'/','SCEN_',s,'/',prefix[[currentWD]],'CI.lower.t_',i,sep='')),silent=TRUE)
    try(load(paste(DATAwd[[currentWD]],'/',i,'/','SCEN_',s,'/',prefix[[currentWD]],'CI.upper.weightVar_',i,sep='')),silent=TRUE)
    try(load(paste(DATAwd[[currentWD]],'/',i,'/','SCEN_',s,'/',prefix[[currentWD]],'CI.lower.weightVar_',i,sep='')),silent=TRUE)
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

# --------------------------
########## Check some values
# Scenario 45 in take 3: n=100, k=10
s <- 45
i <- 350
####### COPE: nothing to weird
load(paste(DATAwd[[currentWD]],'/',i,'/','SCEN_',s,'/',prefix[[currentWD]],'SCOPE_',i,sep=''))
dim(SCOPE)
levelplot(SCOPE[,,,1])
hist(SCOPE[,,,1])
summary(SCOPE[,,,1])
####### TMAP
load(paste(DATAwd[[currentWD]],'/',i,'/','SCEN_',s,'/',prefix[[currentWD]],'STMAP_',i,sep=''))
dim(STMAP)
levelplot(STMAP[,,,1])
hist(STMAP[,,,1])
summary(STMAP[,,,1])
	# Compare with random draws from t-distribution with df=n-1
	check.T <- rt(n=prod(DIM),df=99)
	summary(check.T)
	hist(check.T)

TwoDist.T <- data.frame('Value' = c(STMAP[,,,1],check.T),
											'Source' = factor(rep(c(1,2),each = prod(DIM)),labels=c('tmap','TDist'))
											)

ggplot(TwoDist.T, aes(x=Value, colour=Source)) + geom_density()
	# Seems OK

####### Hedges' G
load(paste(DATAwd[[currentWD]],'/',i,'/','SCEN_',s,'/',prefix[[currentWD]],'SHEDGE_',i,sep=''))
dim(SHEDGE)
summary(SHEDGE[,,,1])
hist(SHEDGE[,,,1])

####### VARCOPE
load(paste(DATAwd[[currentWD]],'/',i,'/','SCEN_',s,'/',prefix[[currentWD]],'SVARCOPE_',i,sep=''))
dim(SVARCOPE)
summary(SVARCOPE[,,,1])
hist(SVARCOPE[,,,1])


####### weighted average with its variance component
load(paste(DATAwd[[currentWD]],'/',i,'/','SCEN_',s,'/',prefix[[currentWD]],'WeightedAvg_',i,sep=''))
length(WeightedAvg)
summary(WeightedAvg)
hist(WeightedAvg)
	# Compare this with a standard normal distribution
	check.wA <- rnorm(n=prod(DIM),sd=1)
	summary(check.wA)
	hist(check.wA)

TwoDist.wA <- data.frame('Value' = c(WeightedAvg,check.wA),
											'Source' = factor(rep(c(1,2),each = prod(DIM)),labels=c('WeightedAverage','NormalDist'))
											)
ggplot(TwoDist.wA, aes(x=Value, colour=Source)) + geom_density()
	# We should compare this with a normal distribution and the variance that is estimated
	load(paste(DATAwd[[currentWD]],'/',i,'/','SCEN_',s,'/',prefix[[currentWD]],'varWeightAvg_',i,sep=''))
	length(varWeightAvg)
	summary(varWeightAvg)
	MeanVar <- mean(varWeightAvg)
	# Take random values
	check.wA2 <- rnorm(n=prod(DIM),sd=sqrt(MeanVar))
	summary(check.wA2)
	hist(check.wA2)
TwoDist.wA2 <- data.frame('Value' = c(WeightedAvg,check.wA2),
											'Source' = factor(rep(c(1,2),each = prod(DIM)),labels=c('WeightedAverage','NormalDistVar'))
											)
ggplot(TwoDist.wA2, aes(x=Value, colour=Source)) + geom_density()







