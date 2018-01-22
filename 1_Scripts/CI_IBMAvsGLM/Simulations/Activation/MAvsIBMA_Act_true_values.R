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
# Also save spatial truth info (smoothed and masked GT area).

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


# Note: all actual parameters are generated in the fMRIGI package.
# See:

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
library(NeuRRoStat)
library(fMRIGI)

# Number of subject: median sample size at 2018 = 28.5 (Poldrack et al., 2017)
nsub <- trueMCvalues('sim_act', 'nsub')

##
###############
### Obtain true noise parameter
###############
##

###################
#### Generate design
###################

# Signal characteristics
TR <- trueMCvalues(ID = 'sim_act', keyword = 'TR')
nscan <- trueMCvalues('sim_act', 'nscan')
total <- trueMCvalues('sim_act', 'total')
on1 <- trueMCvalues('sim_act', 'on1')
onsets <- trueMCvalues('sim_act', 'onsets')
duration <- trueMCvalues('sim_act', 'duration')


# Image characteristics
DIM <- trueMCvalues('sim_act', 'DIM')
voxdim <- trueMCvalues('sim_act', 'voxdim') # Voxelsize
ext <- trueMCvalues('sim_act', 'ext')       #  Extend
nregio <- trueMCvalues('sim_act', 'nregio')

# True % BOLD CHANGE
BOLDC <- trueMCvalues('sim_act', 'BOLDC')

# Generating a design matrix
X <- trueMCvalues('sim_act', 'X')

# Generate time series for ONE active voxel: predicted signal from design matrix
pred <- trueMCvalues('sim_act', 'pred')

# Extend the design matrix with an intercept
xIN <- trueMCvalues('sim_act', 'xIN')

# Contrast: not interested in intercept
CONTRAST <- trueMCvalues('sim_act', 'CONTRAST')

# Calculate (X'X)^(-1) with contrast
design_factor <- trueMCvalues('sim_act', 'design_factor')

# Now we need sensible values for Cohen's d
# We look at Poldrack et al. (2017).
# Here, several contrasts are analyzed into a group analysis of N = 186.
# Then within ROIs, the values for Cohen's d are recorded.
# -------
# See the file Values_d_sigma.R in (https://github.com/NeuroStat/MultivarCBMA/tree/master/1_Scripts/05_Values_d_sigma)
# --- we take the median Cohen's d over several contrasts in fMRI as a medium effect (i.e. 0.55)
# --- the 90% quantile as a high effect (1.02).
# --- and 10% quantile for a low effect (0.14)
TrueD <- trueMCvalues('sim_act', 'TrueD')

# Calculate values for sigma
TrueSigma <- trueMCvalues('sim_act', 'TrueSigma')

# Tau: values come from estimateBSvar.R, no covariate
  # 0th, 50th and 100th percentile of observed between-study variability
Tau <- trueMCvalues('sim_act', 'Tau')

##
###############
### Obtain true Hedges' g
###############
##

# Hedges' g can be obtained by multiplying Cohen's d with the correction factor.
TrueG <- trueMCvalues('sim_act', 'TrueG')


##
###############
### Smoothed and binary GT mask (spatial truth)
###############
##


# True center of activation
TrueLocations <- trueMCvalues('sim_act', 'TrueLocations')

# Spatial smoothing of signal
fwhm <- trueMCvalues('sim_act', 'fwhm')
sigma <- trueMCvalues('sim_act', 'sigma')
width <- trueMCvalues('sim_act', 'width')

# We generate a temporary design for getting a true signal
truthdesign <- trueMCvalues('sim_act', 'truthdesign')

# Now use this to get a sphere shaped area
area <- trueMCvalues('sim_act', 'area')
truth <- trueMCvalues('sim_act', 'truth')

# Unsmoothed ground truth
GroundTruth <- trueMCvalues('sim_act', 'GroundTruth')

# Smooth the GT and put it into the map
SmGT <- trueMCvalues('sim_act', 'SmGT')

# Create the smoothed ground truth mask (where is the true signal)
MaskGT <- trueMCvalues('sim_act', 'MaskGT')

##
###############
### Save objects
###############
##

# List with true parameter values
SavedParam <- list(Nsub = nsub, TrueD = TrueD, TrueSigma = TrueSigma,
               TrueG = TrueG, Tau =  Tau)

# Saving this list, along with spatial truth (smoothed and mask)
saveRDS(SavedParam, file = 'TrueValues.rda')
saveRDS(SmGT, file = 'TrueSmoothedArea.rda')
saveRDS(MaskGT, file = 'TrueMaskedArea.rda')






