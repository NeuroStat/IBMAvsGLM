####################
#### TITLE:     Process al the files from simulations so we have one R object
#### Contents:
####
#### Source Files: HPC - Version
#### First Modified: 16/02/2016
#### Notes:
#################



##
###############
### Notes
###############
##


# As it takes huge time to load in all the R files each time we want to look at the data after a simulation, we will combine them in this file in 1 R object.


# Different takes:
#		* Take 1: multivariate, 7 scenario's
#		* Take 2: univariate (1 voxel), 1 scenario (only white noise)
#		* Take 3: increasing sample size and amount of studies (white noise only)
#		* Take 4: increasing sample size and amount of studies (white noise only) for UNIVARIATE approach
#		* Take 5: Same as take 3, but with R file differently sourced.
#   * Take 6: Simple design (only one condition), grid of 2 x 2 x 2 voxels (white noise only)
#   * Take 7: Image based meta-analysis. Simple design (only one condition), grid of 2 x 2 x 2 voxels (white noise only)



##
###############
### Preparation
###############
##


# Reset working directory
rm(list=ls())
gc(verbose = FALSE)
options(warn = -1)

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
	'Take[7]' = "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/IBMA"
	)

# Prefixes
prefix <- list(
	'Take[1]' = "WNSm",
	'Take[2]' = "UNI_",
  'Take[3]' = 'SK_',
	'Take[4]' = "UNI_SK_",
	'Take[5]' = 'SK_',
  'Take[6]' = 'SD_',
	'Take[7]' = 'IBMA_'
)

NUMDATAwd <- length(DATAwd)
currentWD <- 7

# Number of scenarios
NumScen.tmp <- matrix(c(
                1,7,
                2,1,
                3,45,
                4,45,
                5,45,
                6,30,
								7,30
              ), ncol=2, byrow=TRUE)
NumScen <- NumScen.tmp[currentWD,2]
print(paste('We have ', NumScen, ' scenarios on the menu today!', sep=''))

# Number of conficence intervals
PossibleCIs <- c('norm','t','weightVar')
CIs.tmp <- list(
			'1' = c(1,2,3),
			'2' = c(1,2,3),
			'3' = c(1,2,3),
			'4' = c(1,2,3),
			'5' = c(1,2,3),
			'6' = c(1,2,3),
			'7' = c(2)
			)
CIs <- PossibleCIs[c(CIs.tmp[[currentWD]])]
NumCI <- length(CIs)

# Which type of objects are to be loaded in?
ObjType.tmp <- list(
	'1' = c('CIs','WeightedAvg'),
	'2' = c('CIs','WeightedAvg'),
	'3' = c('CIs','WeightedAvg'),
	'4' = c('CIs','WeightedAvg'),
	'5' = c('CIs','WeightedAvg'),
	'6' = c('CIs','WeightedAvg'),
	'7' = c('CIs')
	)
ObjType <- ObjType.tmp[currentWD]

# Number of objects to be loaded in:
	# * Scenario 1-6: CIs + the weighted average
	# * Scenario 7: only t CI
		# * Second column is for initilizing vectors, third for effective loading in objects
NumObjects.tmp <- matrix(c(
	1, NumCI*2 + 1,
	2, NumCI*2 + 1,
	3, NumCI*2 + 1,
	4, NumCI*2 + 1,
	5, NumCI*2 + 1,
	6, NumCI*2 + 1,
	7, NumCI * 2
	),ncol=2,byrow=TRUE)
NumEffObj <- NumObjects.tmp[currentWD,2]

# Number of executed simulations
nsim.tmp <- matrix(c(
                1,4176,
                2,3000,
                3,350,
                4,1500,
                5,500,
                6,3000,
								7,3000
              ), ncol=2, byrow=TRUE)
nsim <- nsim.tmp[currentWD,2]
print(paste('Each scenario takes ', nsim, ' simulations before developing into a flower!', sep=''))

# Dimension of brain
DIM.tmp <- array(NA, dim=c(NUMDATAwd,3))
	DIM.tmp[c(1,3,5),] <- c(16,16,16)
	DIM.tmp[c(2,4),] <- c(1,1,1)
	DIM.tmp[c(6,7),] <- c(2,2,2)
DIM <- DIM.tmp[currentWD,]


##
###############
### Data Wrangling
###############
##

# How are the objects that we want to load in called:
NameObj.All <- c()
if(any(unlist(ObjType) == 'CIs')){
	NameObj.All <- c(NameObj.All,
				paste('CI.upper.',CIs,'.All', sep=''),
				paste('CI.lower.',CIs,'.All', sep='')
		)
	}
if(any(unlist(ObjType) == 'WeightedAvg')){
	NameObj.All <- c(NameObj.All, 'WeightedAvg.All')
}
NameObj.Sim <- gsub('All','Sim', x = NameObj.All)

# Initialize lists according to the scenario.
# Vectors according to the simulations.
for(o in 1:NumEffObj){
	assign(NameObj.Sim[o], c())
}

# Load in R objects and combine them.
NumErr <- c()
print("Let's go!")
for(s in 1:NumScen){
  print(s)
  for(i in 1:nsim){
		for(o in 1:NumEffObj){
			# Load in objects
			OBJ <- try(load(paste(DATAwd[[currentWD]],'/',i,'/','SCEN_',s,'/',prefix[[currentWD]],gsub('.All', '', NameObj.All[o]),'_',i,sep='')),silent=TRUE)
			# Assign this object to a vector
			try(assign(NameObj.Sim[o], cbind(get(NameObj.Sim[o]), get(OBJ))), silent=TRUE)
			rm(OBJ)
		}
			# Check if objects exist
			EVAL <- !exists(x=gsub('.All', '', NameObj.All))
			if(EVAL){
				 NumErr <- rbind(NumErr,c(s,i))
				 next
			}
	}
	# Save all objects in list, according to the scenario
	for(o in 1:NumEffObj){
		if(s == 1){
			try(assign(NameObj.All[o], list(get(NameObj.Sim[o]))), silent=TRUE)
		}else{
			try(assign(NameObj.All[o], c(get(NameObj.All[o]), o = list(get(NameObj.Sim[o])))), silent=TRUE)
		}
		assign(NameObj.Sim[o], c())
	}
}


# Now to save objects!
try(save(CI.upper.norm.All, file = paste(DATAwd[[currentWD]], '/Take',currentWD, '-CI.upper.norm.All', sep='')),silent=TRUE)
try(save(CI.lower.norm.All, file = paste(DATAwd[[currentWD]], '/Take',currentWD, '-CI.lower.norm.All', sep='')),silent=TRUE)
try(save(CI.upper.t.All, file = paste(DATAwd[[currentWD]], '/Take',currentWD, '-CI.upper.t.All', sep='')),silent=TRUE)
try(save(CI.lower.t.All, file = paste(DATAwd[[currentWD]], '/Take',currentWD, '-CI.lower.t.All', sep='')),silent=TRUE)
try(save(CI.upper.weightVar.All, file = paste(DATAwd[[currentWD]], '/Take',currentWD, '-CI.upper.weightVar.All', sep='')),silent=TRUE)
try(save(CI.lower.weightVar.All, file = paste(DATAwd[[currentWD]], '/Take',currentWD, '-CI.lower.weightVar.All', sep='')),silent=TRUE)
try(save(WeightedAvg.All, file = paste(DATAwd[[currentWD]], '/Take',currentWD, '-WeightedAvg.All', sep='')),silent=TRUE)


	# Also save the missing simulations
	if(is.null(NumErr)){
		cat('NULL', file = paste(DATAwd[[currentWD]], '/Take',currentWD, '-NumErr.txt', sep=''))
	}else{
		colnames(NumErr) <- c('Scenario', 'Simulation')
		write.table(NumErr, file = paste(DATAwd[[currentWD]], '/Take',currentWD, '-NumErr.txt', sep=''),sep='\t',row.names=FALSE,quote=FALSE)
	}




