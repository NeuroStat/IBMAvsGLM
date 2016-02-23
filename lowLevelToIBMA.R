####################
#### TITLE:     IBMA using FSL's FLAME on the COPE and VARCOPES.
#### Contents:
####
#### Source Files:
#### First Modified: 17/02/2016
#### Notes:
#################



##
###############
### Notes
###############
##

# Try using FSL's FLAME directly from R.



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
      SCEN <- 10
    }
  # DataWrite directory: where all files are written to
  DataWrite <- try(as.character(input)[4],silent=TRUE)

# Set starting seed: it is the product of the amount of voxels, the number of studies and the number of subjects!
starting.seed <- 36865*K
set.seed(starting.seed)

# Set WD
if(MACHINE=='HPC'){
  wd <- '/user/scratch/gent/gvo000/gvo00022/vsc40728/Simulation'
}
if(MACHINE=='MAC'){
  wd <- '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/Take7'
}
setwd(wd)

# Give path to FSL
if(MACHINE=='HPC'){
  fslpath <- ''
}
if(MACHINE=='MAC'){
  fslpath <- '/usr/local/fsl/bin/'
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
on1 <- seq(1,total,40)
onsets <- list(on1)
duration <- list(20)
effect.null <- list(0)                              ## No effect
effect <- list(1) 			                            ## Effect of 1 for designmatrix
DIM <- c(2,2,2)


####************####
#### Scenario specific simulation details
####************####


# First make the data frame with the combinations of subjects and studies
OverView <- data.frame('Subjects' = rep(seq(10,100,by=10), 3),
                      'Studies' = rep(c(2,5,10), each= 10))

# Now take the correct amount of subjects and studies (according to SCEN)
nsub <- OverView[SCEN,'Subjects']
nstud <- OverView[SCEN,'Studies']


# At the moment, only consider one noise structure:
	Noise <- list(
	  'S1' = c(1,0,0,0,0,0)
	  )

####************####
#### Subject/Study specific simulation details
####************####
# Subject parameters
TrueLocations <- c(4,4,4)
TrueWhiteNoise <- Noise[1]						# MIND THE INDEX HERE!
TrueRadius <- 1
COPE <- VARCOPE <- array(NA,dim=c(prod(DIM),nsub))


####************####
#### Study parameters
####************####
STCOPE <- STVARCOPE <- array(NA,dim=c(prod(DIM),nstud))


####************####
#### Design matrices
####************####
# Design Matrices via neuRosim:
#     * We need two design vectors:
#     * The first one have an intercept (needed for analysis).
#        * This will be the column of the design matrix in the analysis.
#     * The second one is used to generate data with a NULL effect.
design.Cond1 <- simprepTemporal(onsets = list(on1), durations = list(duration[[1]]),
                       hrf = "double-gamma", TR = TR, totaltime = total,
                       effectsize = list(effect[[1]]))

design.null <- simprepTemporal(regions = 1, onsets = onsets, durations = duration,
                       hrf = "double-gamma", TR = TR, totaltime = total,
                       effectsize = effect.null)

# X-matrix in order to fit the model later on.
x <- matrix(c(simTSfmri(design.Cond1, nscan=nscan, TR=TR, noise="none")),ncol=1)



####************####
#### GENERATE DATA
####************####
# For loop over studies
for(t in 1:nstud){
  print(paste('At study ', t, ', scenario ',SCEN, ' in simulation ', K, sep=''))
  # For loop over nsub
  for(s in 1:nsub){
    # Define two regions (which does nothing as there is no effect, )
    regions <- simprepSpatial(regions = 1, coord = TrueLocations, radius = list(TrueRadius), form ="cube", fading = 0)

    # Weighting structure.                                                                    MIND THE INDEX OF NOISE HERE!!!!!
    #   * Order = white, temporal, low-frequency, physyiological, task related and spatial.
    w <- Noise[[1]]

    # Base value
    base <- 5

    # Actual simulated data
    sim.data <- simVOLfmri(design=design.null, image=regions, base=base, dim=DIM, SNR=0.5,
                 type ="gaussian", noise= "mixture", spat="gaussRF", FWHM=2, weights=w, verbose = TRUE)
      # Transform it to correct dimension (Y = t x V)
      Y.data <- t(matrix(sim.data,ncol=nscan))

      rm(w, sim.data)
    ####************####
    #### ANALYZE DATA: 1e level GLM
    ####************####

    # Fitting GLM model.
    model.lm <- lm(Y.data ~ x)
    b1 <- coef(model.lm)['x',]
    COPE[,s] <- b1

      # Estimate residual (we need to extend the design matrix with an intercept)
      xIN <- cbind(1,x)
      BETA <- coef(model.lm)
      res <- (t(Y.data - xIN %*% BETA) %*% (Y.data - xIN %*% BETA))/nscan - 2
      res <- diag(res)
      # Contrast: not interested in intercept
      CONTRAST <- matrix(c(0,1),nrow=1)
    # Calculate varcope
    VARCOPE[,s] <- CONTRAST %*% (solve(t(xIN) %*% xIN )) %*% t(CONTRAST) %*% res

    # Clean objects
    rm(model.lm, b1, xIN, BETA,res,CONTRAST)

  }
  ####************####
  #### GROUP ANALYSIS: 2e level using FLAME
  ####************####

  # Write auxiliarly files to DataWrite. We need:
    # GRCOPE in nifti
    # GRVARCOPE in nifti
    # 4D mask
    # design.mat file
    # design.grp file
    # design.con file

    #----- 1 ----#
    ### Design.mat
    fileCon <- paste(DataWrite,"/design.mat",sep="")
    # Text to be written to the file
    cat('/NumWaves\t1
    /NumPoints\t',paste(nsub,sep=''),'
    /PPheights\t\t1.000000e+00

    /Matrix
    ',rep("1.000000e+00\n",nsub),file=fileCon)

    #----- 2 ----#
    ### Design.con
    fileCon <- file(paste(DataWrite,"/design.con", sep=""))
    	writeLines('/ContrastName1	Group Average
    /NumWaves	1
    /NumContrasts	1
    /PPheights		1.000000e+00
    /RequiredEffect		5.034

    /Matrix
    1.000000e+00
    ',fileCon)
    close(fileCon)

      #----- 3 ----#
      ### Design.grp
    fileCon <- paste(DataWrite,"/design.grp",sep="")
    # Text to be written to the file
    cat('/NumWaves\t1
    /NumPoints\t',paste(nsub,sep=''),'

    /Matrix
    ',rep("1\n",nsub),file=fileCon)

    #----- 4 ----#
    ### COPE.nii
    GRCOPE4D <- nifti(img=array(COPE,dim=c(DIM,nsub)),dim=c(DIM,nsub),datatype = 16)
    writeNIfTI(GRCOPE4D, filename = paste(DataWrite,'/GRCOPE',sep=''),gzipped=FALSE)

    #----- 5 ----#
    ### VARCOPE.nii
    GRVARCOPE4D <- nifti(img=array(VARCOPE,dim=c(DIM,nsub)),dim=c(DIM,nsub),datatype = 16)
    writeNIfTI(GRVARCOPE4D, filename = paste(DataWrite,'/GRVARCOPE',sep=''),gzipped=FALSE)

    #----- 6 ----#
    ### mask.nii
    mask <- nifti(img=array(1, dim=c(DIM,nsub)), dim=c(DIM,nsub), datatype=2)
    writeNIfTI(mask, filename = paste(DataWrite,'/mask',sep=''),gzipped=FALSE)

    # FSL TIME!
    setwd(DataWrite)
    command <- paste(fslpath, 'flameo --cope=GRCOPE --vc=GRVARCOPE --mask=mask --ld=study',t,'_stats --dm=design.mat --cs=design.grp --tc=design.con --runmode=flame1', sep='')
    Sys.setenv(FSLOUTPUTTYPE="NIFTI")
    system(command)

    # Put the result of pooling subjects in a vector for the COPE and VARCOPE
    STCOPE[,t] <- readNIfTI(paste(DataWrite,"/study",t,"_stats/cope1.nii",sep=""), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
    STVARCOPE[,t] <- readNIfTI(paste(DataWrite,"/study",t,"_stats/varcope1.nii",sep=""), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]

    # Clean up objects
    rm(GRCOPE4D,GRVARCOPE4D,command)
}


####************####
#### META-ANALYSIS: 3e level using FLAME
####************####

# Write auxiliarly files to DataWrite. We need:
  # STCOPE in nifti
  # STVARCOPE in nifti
  # 4D mask
  # design.mat file
  # design.grp file
  # design.con file

  #----- 1 ----#
  ### Design.mat
  fileCon <- paste(DataWrite,"/STdesign.mat",sep="")
  # Text to be written to the file
  cat('/NumWaves\t1
  /NumPoints\t',paste(nstud,sep=''),'
  /PPheights\t\t1.000000e+00

  /Matrix
  ',rep("1.000000e+00\n",nstud),file=fileCon)

  #----- 2 ----#
  ### Design.con
  fileCon <- file(paste(DataWrite,"/STdesign.con", sep=""))
  	writeLines('/ContrastName1	Group Average
  /NumWaves	1
  /NumContrasts	1
  /PPheights		1.000000e+00
  /RequiredEffect		5.034

  /Matrix
  1.000000e+00
  ',fileCon)
  close(fileCon)

  #----- 3 ----#
  ### Design.grp
  fileCon <- paste(DataWrite,"/STdesign.grp",sep="")
  # Text to be written to the file
  cat('/NumWaves\t1
  /NumPoints\t',paste(nstud,sep=''),'

  /Matrix
  ',rep("1\n",nstud),file=fileCon)

  #----- 4 ----#
  ### STCOPE.nii
  STCOPE4D <- nifti(img=array(STCOPE,dim=c(DIM,nstud)),dim=c(DIM,nstud),datatype = 16)
  writeNIfTI(STCOPE4D, filename = paste(DataWrite,'/STCOPE',sep=''),gzipped=FALSE)

  #----- 5 ----#
  ### VARCOPE.nii
  STVARCOPE4D <- nifti(img=array(STVARCOPE,dim=c(DIM,nstud)),dim=c(DIM,nstud),datatype = 16)
  writeNIfTI(STVARCOPE4D, filename = paste(DataWrite,'/STVARCOPE',sep=''),gzipped=FALSE)

  #----- 6 ----#
  ### mask.nii
  mask <- nifti(img=array(1, dim=c(DIM,nstud)), dim=c(DIM,nstud), datatype=2)
  writeNIfTI(mask, filename = paste(DataWrite,'/mask',sep=''),gzipped=FALSE)


# FSL TIME!
setwd(DataWrite)
command <- paste(fslpath, 'flameo --cope=STCOPE --vc=STVARCOPE --mask=mask --ld=MA_stats --dm=STdesign.mat --cs=STdesign.grp --tc=STdesign.con --runmode=flame1', sep='')
Sys.setenv(FSLOUTPUTTYPE="NIFTI")
system(command)


### Now CI around T-value
TVAL <- matrix(readNIfTI(paste(DataWrite,"/MA_stats/tstat1.nii",sep=""), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,],ncol=1)
SE <- sqrt(matrix(readNIfTI(paste(DataWrite,"/MA_stats/varcope1.nii",sep=""), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,],ncol=1))
  # Degrees of freedom:
  tdof_t1 <- readNIfTI(paste(DataWrite,"/MA_stats/tdof_t1.nii",sep=""), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[1,1,1]

CI.upper.t <- TVAL +  (qt(0.975,df=tdof_t1) * SE)
CI.lower.t <- TVAL -  (qt(0.975,df=tdof_t1) * SE)



#### NOTE:
  # Not sure what the difference is with:
    # zflame1lowertstat1.nii
    # zflame1uppertstat1.nii
    # which is in the output after running flame


##
###############
### Save objects
###############
##


save(CI.upper.t, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/IBMA_CI.upper.t_',K,sep=''))
save(CI.lower.t, file=paste(wd,'/Results/',K,'/SCEN_',SCEN,'/IBMA_CI.lower.t_',K,sep=''))








