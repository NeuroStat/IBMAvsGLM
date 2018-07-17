####################
#### TITLE:   Simulation: generate and fit general linear mixed model
#### Contents:
####
#### Source Files:
#### First Modified: 13/07/2018
#### Notes:
#################

##
###############
### Notes
###############
##


# Generate data of N subjects each in M classes.
# Fit a linear mixed model and save the parameters.
# Are we able to extract the true values?


##
###############
### Preparation
###############
##

# Seed
set.seed(1473)

# Load in libraries
library(AnalyzeFMRI)
library(lattice)
library(gridExtra)
library(oro.nifti)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(reshape2)
library(lme4)
library(MASS)
library(RColorBrewer)
library(mvmeta)
library(metafor)
library(devtools)
library(Matrix)
library(neuRosim)
library(NeuRRoStat)
library(fMRIGI)


##
###############
### Functions
###############
##

##
###############
### Simulation parameters
###############
##

# Number of simulations
nsim <- 500

# True values
## Class and sample size
numclass <- 40
NinClass <- 50

## Fixed parameters
beta0 <- 5
beta1 <- 2

## Random parameters
sigma_e <- 1
sigmab0 <- 1
sigmab1 <- 2

# Generate the variance covariance matrix of the random effects
var_cov_U <- rbind(c(sigmab0**2, 0), c(0, sigmab1**2))
# Generate the values for b0 and b1
B_matrix <- MASS::mvrnorm(numclass, mu = c(0,0), Sigma = var_cov_U)


# Generate the predictors of all participants
X <- round(runif(n = NinClass * numclass, min = 1, max = 20))
X <- rep(list(cbind(1, round(runif(n = NinClass, min = 1, max = 20)))), numclass)

# Z-matrix: block diagonal of individual X's
ComplZ <- as.matrix(bdiag(X))

# The V matrix equals: ZDZ' + R
ComplV <- ComplZ %*% var_cov_U %*% t(ComplZ)


ComplX <- ComplV <- matrix(NA,nrow = 1, ncol = 1)
VarCovBeta_raw <- matrix(0, ncol = 2, nrow = 2)
Xlist <- Zlist <- list()
# True variance-covariance matrix for fixed effects parameters
for(i in 1:numclass){
  # Predictors for this class
  X <- cbind(1, round(runif(n = NinClass, min = 1, max = 20)))
  
  # Z-matrix for this class
  Z <- X
  
  # V-matrix
  V <- Z %*% var_cov_U %*% t(Z) +
    diag(sigma_e**2, NinClass)
  
  # Part of var-covar-beta matrix
  VarCovBeta_raw <- VarCovBeta_raw + t(X) %*% solve(V) %*% X
  
  # Add X and V to supermatrix
  #ComplX <- as.matrix(bdiag(ComplX,X))
  #ComplV <- as.matrix(bdiag(ComplV,V))
  
  # Save X and Z
  Xlist[[i]] <- X
  Zlist[[i]] <- Z
}
# Remove NA row in ComplX and ComplV
#ComplX <- ComplX[-1,-1]
#ComplV <- ComplV[-1,-1]

# Now calculate true variance-covariance matrix
#VarCovBeta <- solve(t(ComplX) %*% solve(ComplV) %*% ComplX)
VarCovBeta <- solve(VarCovBeta_raw)

# The true variance-covariance matrix of the parameter estimates 
# of the fixed effects then becomes:
#TrueV <- cbind(1,X) %*% var_cov_U %*% t(cbind(1,X)) +
  diag(sigma_e**2, NinClass)
#VarCovBeta <- solve(numclass * t(cbind(1,X)) %*% solve(TrueV) %*% cbind(1,X))

# Standard error of beta = sqrt(var(beta))
SEBeta <- data.frame('term' = c('(Intercept)', 'X'),
                     'TrueSE' = sqrt(diag(VarCovBeta)), 
                     stringsAsFactors = FALSE)


##
###############
### Simulation
###############
##

# Empty data frame with simulation results
FitTotDat <- data.frame() %>% as_tibble()

# For loop over the simulations
for(r in 1:nsim){
  # Print simulation
  print(r)
  
  # Empty data frame
  TotDat <- data.frame()
  
  # Empty true var-covar matrix of fixed parameter estimates
  # VarCovBeta <- diag(0,2)
  
  # Loop over the classes
  for(i in 1:numclass){
    
      # Generate predictor for subjects in this class
      # X <- runif(n = NinClass, min = 1, max = 20)
      
      # Calculate non-inverted true variance-covariance matrix for this class
      # TrueV <- cbind(1,X) %*% var_cov_U %*% t(cbind(1,X)) +
      #   diag(sigma_e**2, NinClass)
      # # Also: add it to all subjects (later on invert it)
      # VarCovBeta <- t(cbind(1,X)) %*% solve(TrueV) %*% cbind(1,X) +
      #   VarCovBeta

      # Loop over the subjects
      for(j in 1:NinClass){
      # Generate data using: X*beta + Z*u + e
      # dat <- cbind(1,X) %*% matrix(c(beta0, beta1), ncol = 1) + 
      #   cbind(1,X) %*% matrix(c(B_matrix[i,1], B_matrix[i,2]), ncol = 1) +
      #   rnorm(n = NinClass, mean = 0, sd = sigma_e)
        
        dat <- Xlist[[i]] %*% matrix(c(beta0, beta1), ncol = 1) + 
          Zlist[[i]] %*% matrix(c(B_matrix[i,1], B_matrix[i,2]), ncol = 1) +
          rnorm(n = NinClass, mean = 0, sd = sigma_e)
      
      # Add to data frame
      TotDat <- data.frame(Y = dat, X = Xlist[[i]][,2], class = i) %>% as_tibble() %>%
        bind_rows(TotDat,.)
    }
  }
  
  # Invert the true variance-covariance matrix
  # VarCovBetaT <- solve(VarCovBeta)
  
  # Standard error of beta = sqrt(var(beta))
  # SEBeta <- data.frame('term' = c('(Intercept)', 'X'),
  #            'TrueSE' = sqrt(diag(VarCovBetaT)), stringsAsFactors = FALSE)
  
  # Analysis
  fit <- lmer(Y ~ 1 + X + (1 + X|class), data = TotDat, REML = TRUE)
  FitTotDat <- broom::tidy(fit) %>% 
    # Add true SE
    left_join(.,SEBeta, by = 'term') %>%
    mutate(sim = r) %>%
    bind_rows(FitTotDat,.)
   
}


##
###############
### Analysis
###############
##

# FitTotDat %>% 
#   mutate(filEst = ifelse(term %in% c('(Intercept)', 'X'),
#                 estimate, NA)) %>%
#   group_by(term) %>%
#   summarise(AvgEst = mean(estimate),
#             AvgSE = mean(std.error),
#             EmpSE = sd(filEst),
#             AvgTrueSE = mean(TrueSE))

# Fixed part
FitTotDat %>% 
  filter(term %in% c('(Intercept)', 'X')) %>%
  # Add true estimate
  left_join(.,data.frame(term = c('(Intercept)', 'X'),
                         TrueEst = c(beta0, beta1), stringsAsFactors = FALSE),
            by = 'term') %>%
  group_by(term) %>%
  summarise(AvgEst = mean(estimate),
            TrueEst = mean(TrueEst),
            AvgSE = mean(std.error),
            TrueSE = mean(TrueSE))

# Random part
FitTotDat %>% 
  filter(term %in% c('sd_Observation.Residual',
                     'sd_(Intercept).class',
                     'sd_X.class')) %>%
  dplyr::select(term, estimate, sim) %>%
  left_join(.,data.frame('term' = c('sd_Observation.Residual',
                                    'sd_(Intercept).class',
                                    'sd_X.class'),
                         'TrueSD_ran' = c(sigma_e, sigmab0, sigmab1),
                         stringsAsFactors = FALSE),
            by = 'term') %>%
  group_by(term) %>%
  summarise(AvgEst = mean(estimate),
            AvgTrueSD = mean(TrueSD_ran))




