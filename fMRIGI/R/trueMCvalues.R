#' @title
#' True parameter values.
#'
#' @description
#' This function returns true parameter value for the Monte-Carlo simulations
#' It takes an ID for the scenario (e.g. data with activation) and keyword as argument.
#' It is used to define global variables once (here). The scripts running the
#' simulations, analyses, visualizations etc., rely on this function.
#' The idea is to reduce the amount of copying global variables between scripts.
#'
#' @param ID An identifier to select the appropriate simulation scenario.
#' @param keyword The parameter that needs to be selected.
#'
#' @details For ID, we have the following possibilities at the moment:
#' * act_sim
#' @md
#'
#' @export
trueMCvalues <- function(ID = c('sim_act'), keyword){
  # Check ID
  ID <- match.arg(ID)

  # Objects for different simulation settings
  if(ID == 'sim_act'){
    ###################
    #### Global variables
    ###################

    # Number of subject: median sample size at 2018 = 28.5 (Poldrack et al., 2017)
    nsub <- 29

    # Number of studies: not backed by empirical data.
    nstud <- seq(5, 50, by = 5)

    ###################
    #### Data characteristics
    ###################

    # Signal characteristics
    TR <- 2
    nscan <- 200
    total <- TR*nscan
    on1 <- seq(1,total,40)
    onsets <- list(on1)
    duration <- list(20)

    # Image characteristics
    DIM <- c(9,9,9)
    voxdim <- c(1,1,1)          # Voxelsize
    ext <- 1                    #  Extend
    nregio <- 1

    ###################
    #### Generate a design: GROUND TRUTH DESIGN
    ###################

    # %BOLD change => fixed quantity
    #   We will change the amount of noise to change effect size, Cohen's d
    # See MAvsIBMA_Act_true_values.R on how we obtained values for Cohen's d
    #   and the amount of noise within subjects to achieve these ES.
    # We have a conditition with a null effect, and one with an effect.
    BOLDC <- c(0,3)

    # Base of signal
    base <- 100

    # Generating a design matrix
    X <- neuRosim::simprepTemporal(total,1,onsets=onsets,effectsize = 1, durations=duration,
                         TR = TR, acc=0.1, hrf="double-gamma")

    # Generate time series for ONE active voxel: predicted signal from design matrix
    pred <- neuRosim::simTSfmri(design=X, base=0, SNR=1, noise="none", verbose=FALSE)
    # plot(pred, type = 'l')

    # Extend the design matrix with an intercept
    xIN <- cbind(base,pred)

    # Contrast: not interested in intercept
    CONTRAST <- matrix(c(0,1),nrow=1)

    # Calculate (X'X)^(-1) with contrast
    design_factor <- CONTRAST %*% (solve(t(xIN) %*% xIN )) %*% t(CONTRAST)

    ###################
    #### True effect sizes + variability (within- and between-subject)
    ###################

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

    # The true hedges g is obtained by multiplying TrueD with the h correction factor.
    TrueG <- TrueD * NeuRRoStat::corrH(Ne = nsub, type = 'one.sample')

    # Our goal is now to define sensible vaules for within- and between-subject variability.
    # --------------------------------------------------------------------------------
    # Denote within-subject variance as: TrueSigma2W
    # And between-subject variance as: TrueSigma2B
    # Consider the following simplified GLM in fMRI:
    # Y = beta0 + beta1X + E (--> first level)
    # Y_G = beta1*X_G + E* (--> second level)
    # Hence Var(Y_G) = Var(E*) = sigma^2G = sigma^2_B + sigma^2*design_factor
    # Hence, the variance at second level is sum of between-subject variability (sigma^2_B)
        # and imperfect first level intra-subject estimation variability (sigma^2 = white noise)
    # Furthermore we have: d = beta/sigma_G
    # So we get for sigma_G:
    TrueSigmaG <- BOLDC[2]/TrueD
    # Now we want that sigma^2_B / sigma^2_W = 0.5.
      # Work this out and you get:
    # ---------------------------------------------------------------------------------
    TrueSigma2B <- TrueSigmaG^2/3
    TrueSigma2W <- (2*TrueSigmaG^2 / 3) / c(design_factor)
    # Written otherwise: solve(design_factor) %*% (TrueSigmaG^2 - (TrueSigmaG^2/3))
    # CHECK: TrueSigma2B / (design_factor %*% TrueSigma2W)
    # CHECK: (design_factor %*% TrueSigma2W + TrueSigma2B) == TrueSigmaG^2
    # CHECK: TrueSigma2B/TrueSigma2W
    # Furthermore, the true value for Cohen's d should be:
    # ----> BODLC/(sqrt(design_factor %*% TrueSigma2W + TrueSigma2B))
    # CHECK: BOLDC[2]/sqrt((design_factor %*% TrueSigma2W + TrueSigma2B))

    ###################
    #### True variability between-studies
    ###################

    # Tau: values come from estimateBSvar.R, no covariate
    # 0th, 50th and 100th percentile of observed between-study variability
    Tau <- sqrt(c(0,0.10,0.495))
    # I^2 is the observed excessive dispersion over the total amount of variability
    # Thus we have: TrueSigma2M / (TrueSigma2M + TrueSigma2B + TrueSigma2W) == I2
    I2 <- c(0, 62.61, 87.28) / 100
    # Work this out and we get:
    TrueSigma2M <- (I2*TrueSigma2W + I2*TrueSigma2B) / (1 - I2)
    #(I2*TrueSigma2W*c(design_factor) + I2*TrueSigma2B*c(design_lvl2)) / (1 - I2)
    # CHECK: TrueSigma2M / (TrueSigma2M + TrueSigma2B + TrueSigma2W)
    #Xg <- matrix(1, nrow = 10)
    #design_lvl2 <- solve(t(Xg) %*% Xg)
    #(I2*TrueSigma2W + I2*TrueSigma2B) / (1 - I2)
    #test <- (I2*TrueSigma2W*c(design_factor) + I2*TrueSigma2B*c(design_lvl2)) / (1 - I2)
    #test <- (I2*TrueSigma2W*c(design_factor) + I2*TrueSigma2B) / (1 - I2)
    #sqrt(test)

    # Group design matrix
    Xg <- matrix(1, nrow = nsub)
    design_lvl2 <- solve(t(Xg) %*% Xg)
    var2LVLS <- (TrueSigma2B + TrueSigma2W * c(design_factor)) * c(design_lvl2)
    #var2LVLS <- (TrueSigma2B * c(design_factor_lvl2) + TrueSigma2W * c(design_factor))
    COMBs <- expand.grid(var2LVLS, TrueSigma2M)
    COMBs$TrueDs <- BOLDC[2] / sqrt(COMBs$Var1 + COMBs$Var2)
    COMBs
  }

  # Check whether keyword matches one of the objects
  objCheck <- try(get(keyword), silent = TRUE)
  if(grepl(pattern = 'not found', x = objCheck[1])) {
    stop('Keyword not found in settings (see file trueMCvalues.R)')
  }

  # Return the object
  return(get(keyword))
}
