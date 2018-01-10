####################
#### TITLE:     Define true values for simulations
#### Contents:
####
#### Source Files:
#### First Modified: 09/01/2018
#### Notes:
#################


##
###############
### Notes
###############
##

# Script to calculate the true values for the noise, effect sizes and weighted averages.

# Note on calculating the noise (sigma squared of the model) of the time series:
# --- Let us use 3 %BOLD change as a fixed quantity.
# --- We know values for Cohen's d from the Poldrack et al. (2017) study.
# --- These are calculated using 186 unrelated subjects.
# --- Hence we have: d = COPE/[sqrt(VARCOPE) * sqrt(186)]
# --- This is to obtain: d = average effect / sigma.
# --- Sigma in this case corresponds to between subject variability, not the
# ---   noise in a within-subject time series.
# --- To get the latter, note that between subject variability corresponds to 
# ---   the variance of the first level beta.
# --- In matrix notation, we have for the GLM:
# ---   var(beta) = sigma^2 (X'X)^(-1)
# --- Thus, we also need to calculate (X'X)^(-1), which depends on the design.
# --- For this reason, we generate a design here. 
# --- This is a blocked design, in unity (peak = 1).


##
###############
### Preparation
###############
##

# Load in libraries
library(AnalyzeFMRI)
library(lattice)
library(gridExtra)
library(oro.nifti)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(Hmisc)
library(devtools)
library(neuRosim)



##
###############
### Obtain true noise parameter
###############
##

###################
#### Generate design 
###################

# Signal characteristics
TR <- 2
nscan <- 200
total <- TR*nscan
on1 <- seq(1,total,40)
onsets <- list(on1)
duration <- list(20)

# True % BOLD CHANGE
BOLDC <- 3

# Generating a design matrix
X <- simprepTemporal(total,1,onsets=onsets,effectsize = 1, durations=duration,
                     TR = TR, acc=0.1, hrf="double-gamma")

# Generate time series for ONE active voxel: predicted signal from design matrix
pred <- simTSfmri(design=X, base=100, SNR=1, noise="none", verbose=FALSE)
# plot(pred, type = 'l')

# Extend the design matrix with an intercept
xIN <- cbind(1,pred)

# Contrast: not interested in intercept
CONTRAST <- matrix(c(0,1),nrow=1)

# Calculate (X'X)^(-1) with contrast
design_factor <- CONTRAST %*% (solve(t(xIN) %*% xIN )) %*% t(CONTRAST)

# Now we need sensible values for Cohen's d
# We look at Poldrack et al. (2017).
# Here, several contrasts are analyzed into a group analysis of N = 186.
# Then within ROIs, the values for Cohen's d are recorded.
# -------
# See the file Values_d_sigma.R in (https://github.com/NeuroStat/MultivarCBMA/tree/master/1_Scripts/05_Values_d_sigma)
# --- we take the median Cohen's d over several contrasts in fMRI as a medium effect (i.e. 0.55)
# --- the 90% quantile as a high effect (1.02).
# --- and 10% quantile for a low effect (0.14)
TrueD <- c(0.14, 0.55, 1.02)

# Calculate values for sigma
TrueSigma <- BOLDC/(TrueD * as.vector(sqrt(design_factor)))


##
###############
### Obtain true Hedges' g
###############
##

# Hedges' g can be obtained by multiplying Cohen's d with the correction factor.
# This depends on amount of subjects. I don't know the amount of subjects yet,
# so I need to wait with this. 



##
###############
### Save objects
###############
##

saveRDS(TrueSigma, file = 'TrueSigma.rda')







