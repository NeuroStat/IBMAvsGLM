####################
#### TITLE:     Estimate between-study variability in fMRI database about 'pain'
#### Contents:
####
#### Source Files:  //BSVar
#### First Modified: 12/01/2018
#### Notes:
#################

# tmp wd
setwd("~/Dropbox/PhD/PhDWork/Meta Analysis/R Code/Studie_Simulation/SimulationGit/1_Scripts/CI_IBMAvsGLM/Simulations/Activation/BSVar")

##
###############
### Notes
###############
##

#   In this script, we estimate the observed amount of between-study heterogeneity.
#   We have a database containing 32 studies measuring PAIN > NO PAIN.
#   In total 631 subjects.
#   ---
#   Data/* folder contains the processed (often upsampled to MNI 2 mm standard template)
#         t-maps of the individual studies.
#   In the database.csv folder, we have the name of the nifti files with the sample size.
#   ---
#   First we load in the data, then convert to Hedges'g
#   Next we mask the data using only regions of interest determined using Neurosynth (reverse inference)
#   Within these regions, we estimate between-study heterogeneity.


##
###############
### Preparation
###############
##

# Libraries
library(tidyverse)
library(purrr)
library(oro.nifti)
library(mvmeta)
library(metafor)
library(ggridges)
library(NeuRRoStat)

# Custom function
ExtractTau <- function(Y, V){
  fit <- rma(yi = Y, vi = V, method = "DL")
  return(fit$tau2)
}

# Number of studies
nstud <- 32

##
###############
### Read in data
###############
##

# Start with database file
database <- read.csv2('database.csv', header = TRUE, stringsAsFactors = FALSE)

# Read in ROI mask 
ROI <- readNIfTI('pain_pFgA_z_FDR_0.01_mask')[,,]
# Dimension in 3D
DIM3D <- dim(ROI)
# Switch to array
ROI <- array(ROI, dim = prod(DIM3D))

# Number of masked voxels
table(ROI)

# Check if number of studies in database corresponds to nstud
if(dim(database)[1] != nstud) stop('ERROR: number of studies does not correspond to number found in database')

# Empty data matrix
allStud <- matrix(NA, nrow = prod(DIM3D), ncol = nstud)

# For loop over the studies
for(i in 1:nstud){
  # Progress
  print(i)
  # Name of dataset is first column of database
  studDat <- readNIfTI(paste('Data/', database[i,'img'], '.nii', sep = ''))[,,]
    # Check if dimensions match
    if(all(dim(studDat) != DIM3D)) stop(paste0('Dimension of image ', i,
                                               'does not match MNI space'))
  # Switch to vector
  studDat <- array(studDat, dim = prod(DIM3D))
  # Values outside ROI to NA
  studDat[ROI == 0] <- NA
  # Bind to data matrix
  allStud[,i] <- studDat
  
  # Reset
  rm(studDat)
}

# Summary
summary(allStud)

# Distribution of all t-values
distrAllT <- data.frame(allStud) %>% as.tibble()
  names(distrAllT) <- paste('S', 1:nstud, sep = '')
distrAllT <- gather(distrAllT, key = 'study', value = 'Tvalue')
# Raw violin plot
ggplot(distrAllT, aes(x = study, y = Tvalue)) +
  geom_violin()
# Raw ridge plot
ggplot(distrAllT, aes(y = study, x = Tvalue)) +
  geom_density_ridges2()

# Now only for masked values
distrAllT <- data.frame(distrAllT, ROImask = rep(ROI, nstud)) %>% as.tibble()
maskedVox <- filter(distrAllT, ROImask == 1)
  # Create factor of study
  maskedVox$study <- factor(maskedVox$study, levels = paste('S', 1:nstud, sep = ''))
  # Add variable for Kragel, Maumet and other study
  IDvar <- c(rep('Other',4), rep('Kragel', 6), 'Other', rep('Maumet', 21))
  maskedVox$IDvar <- factor(rep(IDvar, each = sum(ROI)))
# Violin plot
ggplot(maskedVox, aes(x = study, y = Tvalue)) +
  geom_violin(aes(fill = IDvar)) +
  scale_y_continuous('T-value') +
  scale_x_discrete('Study') +
  scale_fill_brewer('ID study*', type = 'qual', palette = 3) +
  ggtitle("T-values from voxels in ROI") + 
  labs(caption = "*t-values from Kragel are obtained through OLS. \n Studies from Maumet are from same site.")

# Ridge plot
ggplot(maskedVox, aes(y = study, x = Tvalue)) +
  geom_density_ridges2(aes(fill = IDvar)) +
  scale_y_discrete('Study', expand = c(0.01, 0)) +
  scale_x_continuous('T-value', expand = c(0.01, 0)) +
  scale_fill_brewer('ID study*', type = 'qual', palette = 3) +
  ggtitle("T-values from voxels in ROI") +
  labs(caption = "*t-values from Kragel are obtained through OLS. \n Studies from Maumet are from same site.") +
  theme_ridges()

  
##
###############
### Hedges' g
###############
##

# Transform the t-values to Hedges' g and calcualte variance
  # First add sample size to dataset 
dataG <- database %>% as.tibble() %>%
  mutate(study = factor(paste('S', 1:nstud, sep = ''),
                        levels = paste('S', 1:nstud, sep = ''))) %>%
  select(-img) %>% right_join(., maskedVox, by = 'study') %>%
  arrange(samplesize) %>% 
  # Now transform to hedges' g and calculate variance
  mutate(hedgeG = hedgeG(t = Tvalue, N = samplesize),
         varG = varHedge(g = hedgeG, N = samplesize))

# Ridge plot
# First re order study factor according to sample size
OrderFactorSS <- database %>% as.tibble() %>%
  mutate(study = factor(paste('S', 1:nstud, sep = ''),
            levels = paste('S', 1:nstud, sep = ''))) %>%
            select(-img) %>% arrange(samplesize)
dataG$study <- factor(dataG$study, levels = rev(OrderFactorSS$study))

ggplot(dataG, aes(y = study, x = hedgeG)) +
  geom_density_ridges2(aes(fill = samplesize), scale = 2) +
  scale_y_discrete('Study', expand = c(0.01, 0)) +
  scale_x_continuous("Hedges' g", expand = c(0.01, 0)) +
  scale_fill_gradient2('Sample size', 
                       low = '#3288bd',
                       mid = '#ffffbf',
                       high = '#d53e4f') +
  ggtitle("Hedges' g from voxels within ROI", 
          subtitle = 'Sorted according to study sample size') +
  theme_ridges()


##
###############
### Meta-analysis
###############
##

# First add voxel ID to dataframe
voxID <- dataG %>% mutate(voxel = rep(1:sum(ROI), nstud))

# # Now estimate tau using the function ExtractTau
# t1 <- Sys.time()
# estTau <- voxID %>% select(-IDvar, -Tvalue) %>% group_by(voxel) %>% 
#   summarise(tau2 = ExtractTau(Y = hedgeG, V = varG))
# Sys.time() - t1       # Time difference of 59.65802 secs

# # Or using purrr
# t2 <- Sys.time()
# estTauP <- voxID %>% split(.$voxel) %>% 
#   map(~ rma(yi = hedgeG, vi = varG, method = 'DL', data = .)) %>%
#   map_dbl("tau2") %>% data.frame(tau2 = .) %>% as.tibble()
# Sys.time() - t2         # Time difference of 1.118482 mins

# Although purrr is a bit slower, it is possible to conveniently extract multiple coefficients.
# We also extract I2, H2 and the beta estimate (which is the weighted average).
estTau <- voxID %>% split(.$voxel) %>% 
  map(~ rma(yi = hedgeG, vi = varG, method = 'DL', data = .)) %>%
  map_df(., magrittr::extract, c("tau2", "I2", "H2", "beta"))

# Plot the observed heterogeneity
summary(estTau)
estTau %>% gather(key = 'Parameter', value = 'Value') %>%
  ggplot(., aes(x = Parameter, y = Value)) +
  geom_violin(aes(fill = Parameter), alpha = 0.5) +
  scale_fill_brewer('', type = 'qual', palette = 5) +
  facet_wrap(~ Parameter, scales = 'free', ncol = 1) +
  coord_flip() +
  scale_x_discrete('Parameter', labels = c('beta' = expression(beta),
        'H2'   = expression(H^2),
        'I2'   = expression(I^2),
        'tau2'   = expression(tau^2))) +
  guides(fill = FALSE)  + 
  ggtitle('Observed distributions for various parameters inside ROI') +
  theme_bw() +
  theme(strip.background = element_blank(),
    strip.text.x = element_blank())







# Older code
testyi <- filter(voxID, voxel == 10000)$hedgeG
testvi <- filter(voxID, voxel == 10000)$varG
rma(yi = testyi, vi = testvi, method = "DL")$tau2
str(rma(yi = testyi, vi = testvi, method = "DL"))
rma(yi = testyi, vi = testvi, method = "DL")[['tau2']]
rma(yi = testyi, vi = testvi, method = "DL")$beta


mtcars %>%
  split(.$cyl) %>% # from base R
  map(~ lm(mpg ~ wt, data = .)) %>%
  map(summary) %>%
  map_dbl("r.squared")

test <- voxID %>% select(-IDvar, -Tvalue) %>% group_by(study) %>% 
  sample_n(size = 10)
test %>% split(.$voxel) %>% 
  map(~ rma(yi = hedgeG, vi = varG, data = .)) %>%
  map_dbl("tau2")
  
  group_by(voxel) %>% 
  summarise(tau = ExtractTau(Y = hedgeG, V = varG))

  
  test <- voxID %>% select(-IDvar, -Tvalue) %>%
    filter(voxel %in% sample(1:sum(ROI), size = 10))
  test %>% split(.$voxel) %>% 
    map(~ rma(yi = hedgeG, vi = varG, method = 'DL', data = .)) %>%
    map_dbl(c("tau2", "I2", "H2"))
  test %>% group_by(voxel) %>% 
    summarise(tau2 = ExtractTau(Y = hedgeG, V = varG))
  
  test %>% split(.$voxel) %>% 
    map(~ rma(yi = hedgeG, vi = varG, method = 'DL', data = .)) %>%
    map_df(., magrittr::extract, c("tau2", "I2", "H2"))
      
  