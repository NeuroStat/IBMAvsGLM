####################
#### TITLE:     Generate ES's and fit random effects model: consistent estimator?
#### Contents:
####
#### Source Files:
#### First Modified: 11/03/2019
#### Notes:
#################



##
###############
### Notes
###############
##

# Let us generate ES's and then fit a random effects model.
# Do we get the expected values?

# Using model: Y = x*beta + u + e


##
###############
### Preparation
###############
##

# Libraries
library(metafor)
library(dplyr)
library(tidyr)
library(ggplot2)
library(NeuRRoStat)

# Seed
set.seed(245)

# Number of simulations
nsim <- 1000

##
###############
### Custom functions
###############
##

# Estimator for B
EstB <- function(X, W, Y){
  estimate <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% Y
  return(estimate)
}

##
###############
### Simulation parameters
###############
##

# True ES
Tes <- 0.9

# Number of studies
nstud <- 200

# X vector of 1's
X <- matrix(1, nrow = nstud)

# Between study variability
sigm2U <- c(0, 1, 4, 9)

# Within study-variability: for each study a value. 
# Based on number of subjects in each study. 
# Number of subjects in each study
nsubS <- sample(x = 50:100, size = nstud, replace = TRUE)
#nsubS <- rep(30, nstud)
sigm2E <- NeuRRoStat::varHedgeT(g = Tes, N = nsubS)

##
###############
### Generate data
###############
##

# Empty data frame
results <- data.frame() %>% as_tibble()

# Progress
PROG <- seq(0, nsim, length.out = 10)

# For loop over the simulations
for(i in 1:nsim){
  # Progress
  if(i %in% PROG) print(paste0('At ', i/nsim * 100, '%'))
  # For loop over the amount of between-study heterogeneity
  for(j in 1:length(sigm2U)){
    # Generate data
    Y <- X %*% Tes + rnorm(n = nstud, mean = 0, sd = sqrt(sigm2U[j])) +
      rnorm(n = nstud, mean = 0, sd = sqrt(sigm2E))
    
    # Weights of each study
    vi <- NeuRRoStat::varHedgeT(g = Y, N = nsubS)
    
    # Fit the model using rma + HE or REML estimator for tau2
    fitHE <- rma(yi = Y, vi = vi, method = 'HE')
    fitREML <- rma(yi = Y, vi = vi, method = 'REML')
    
    # Between study variability: DL estimator
    tau2_DL <- NeuRRoStat::tau(Y = Y, W = (1/vi), k = length(Y))
    # HE estimator
    tau2_HE <- c(fitHE$tau2)
    # REML estimator
    tau2_REML <- c(fitREML$tau2)

    # Weights
    W <- diag(c(1/(vi + tau2_DL)))
    
    # Estimate for B using various estimators for tau
    # First one = DerSimonian and Laird
    B_DL <- EstB(X = X, W = W, Y = Y)
      # controle: rma(yi = Y, vi = vi, method = 'DL') --> klopt
    # One using Hedges
    B_HE <- c(fitHE$beta)
    # One using REML
    B_REML <- c(fitREML$beta)
    
    # Save in data frame
    results <- data.frame('B' = c(B_DL, B_HE, B_REML),
              'TrueB' = Tes,
              'tau2' = c(tau2_DL, tau2_HE, tau2_REML),
              'TrueTau2' = sigm2U[j],
              'estimatorBS' = c('DL', 'HE', 'REML')
               ) %>% as_tibble() %>%
      bind_rows(results, .)
  }
}

# Average over simulations
results %>%
  group_by(estimatorBS, TrueTau2) %>%
  summarise(AvgEstB = mean(B),
            TrueB = mean(TrueB),
            AvgEstTau2 = mean(tau2))

# True value using the estimator ==> same as expected value
TrueW <- diag(1/(sigm2E + sigm2U[2]))
TrueY <- matrix(Tes, nrow = nstud)
EstB(X = X, W = TrueW, Y = TrueY)

# Make a plot for B
results %>%
  group_by(estimatorBS, TrueTau2) %>%
  summarise(AvgEstB = mean(B),
            TrueB = mean(TrueB),
            AvgEstTau2 = mean(tau2)) %>%
  group_by(estimatorBS, TrueTau2) %>%
  ggplot(aes(x = estimatorBS, y = AvgEstB)) + 
  geom_bar(stat = 'identity',
                 position = 'dodge2',
                 aes(fill = factor(TrueTau2))) +
  geom_hline(aes(yintercept = Tes), colour = 'black', linetype = 'dashed') +
  labs(caption = paste0('Averaged over ', nsim, ' simulations. Dashed line = true value.')) +
  scale_fill_brewer('True between-study variance', type = 'qual', palette = 5) +
  scale_y_continuous('Average estimate of B') +
  scale_x_discrete('Estimator of between-study variability') +
  theme_bw() +
  theme(panel.grid.major = element_line(size = 0.8),
        panel.grid.minor = element_line(size = 0.8),
        panel.spacing = unit(1, "lines"),
        axis.title.x = element_text(face = 'plain'),
        axis.title.y = element_text(face = 'plain'),
        axis.text = element_text(size = 9, face = 'bold'),
        axis.ticks = element_line(size = 0.9),
        axis.ticks.length=unit(.20, "cm"),
        axis.line = element_line(size = .75),
        title = element_text(face = 'plain'),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(face = 'bold', size = 12),
        legend.position = 'bottom',
        legend.text = element_text(face = 'bold'))


# Make a plot for between-study heterogeneity
results %>%
  dplyr::select(-B,-TrueB) %>%
  group_by(TrueTau2, estimatorBS) %>%
  summarise(AvgEstTau2 = mean(tau2)) %>%
  gather(value = 'value', key = 'type', -2) %>%
  mutate(type_rename = ifelse(type == 'TrueTau2', ' True value',
                              'Estimated value')) %>%
  mutate(type_factor = factor(type_rename)) %>%
  mutate(TrueValue = rep(rep(sigm2U, each = 3), 2)) %>%
  group_by(estimatorBS, type_factor) %>%
  ggplot(aes(x = factor(TrueValue), y = value)) + 
  geom_bar(stat = 'identity',
           width = .9,
           position = 'dodge',
           aes(fill = type_factor)) +
  facet_wrap(~estimatorBS) + 
  labs(caption = paste0('Averaged over ', nsim, ' simulations.')) +
  scale_fill_brewer('True between-study variance', type = 'qual', palette = 3) +
  scale_y_continuous('Average estimate of between-study variability') +
  scale_x_discrete('Amount of between-study variability') +
  theme_bw() +
  theme(panel.grid.major = element_line(size = 0.8),
        panel.grid.minor = element_line(size = 0.8),
        panel.spacing = unit(1, "lines"),
        axis.title.x = element_text(face = 'plain'),
        axis.title.y = element_text(face = 'plain'),
        axis.text = element_text(size = 9, face = 'bold'),
        axis.ticks = element_line(size = 0.9),
        axis.ticks.length=unit(.20, "cm"),
        axis.line = element_line(size = .75),
        title = element_text(face = 'plain'),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(face = 'bold', size = 9),
        legend.position = 'bottom',
        legend.text = element_text(face = 'plain'))






































