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
    voxdim <- c(3.5, 3.5, 3.51) # Voxelsize
    ext <- 1                    #  Extend
    nregio <- 1

    ###################
    #### Generate a design: GROUND TRUTH DESIGN
    ###################

    # %BOLD change => fixed quantity
    #   We will change the amount of noise to change effect size, Cohen's d
    # See MAvsIBMA_Act_true_values.R on how we obtained values for Cohen's d
    #   and the amount of noise within subjects to achieve these ES.
    BOLDC <- 3

    # Base of signal
    base <- 100

    # Spatial smoothing of signal
    fwhm <- 8
    sigma <- fwhm/sqrt(8*log(2))
    width <- 5

    # Generating a design matrix
    X <- neuRosim::simprepTemporal(total,1,onsets=onsets,effectsize = 1, durations=duration,
                         TR = TR, acc=0.1, hrf="double-gamma")

    # Generate time series for ONE active voxel: predicted signal from design matrix
    pred <- neuRosim::simTSfmri(design=X, base=100, SNR=1, noise="none", verbose=FALSE)
    # plot(pred, type = 'l')

    # Extend the design matrix with an intercept
    xIN <- cbind(1,pred)

    # Contrast: not interested in intercept
    CONTRAST <- matrix(c(0,1),nrow=1)

    # Calculate (X'X)^(-1) with contrast
    design_factor <- CONTRAST %*% (solve(t(xIN) %*% xIN )) %*% t(CONTRAST)

    ###################
    #### True effect sizes + variability (within and between)
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

    # Now we have in fMRI:
    # Y = beta0 + beta1X + E (first level)
    # Y_G = beta1*X_G + E* (second level)
    # Hence Var(Y_G) = Var(E*) = sigma^2G = sigma^2_B + sigma^2(X'X)^(-1)
    # Hence, the variance at second level is sum of between-subject variability (sigma^2_B)
        # and imperfect first level intra-subject estimation variability (sigma^2 = white noise)
    # Furthermore we have: d = beta/sigma_G
    # So we get for sigma_G:
    TrueSigmaG <- BOLDC/TrueD
    # Now we want that sigma^2_B / sigma^2_W = 0.5.
      # Work this out and you get:
    TrueSigma2B <- TrueSigmaG^2/3
    TrueSigma2W <- (2*TrueSigmaG^2 / 3) / c(design_factor)
    # NOTE: TrueSigma2B corresponds to variance of the random slopes between-subjects
    #   We also want a random intercept (B0), variance = 1
    Trueb0_Bsub <- 1

    # Calculate values for sigma
    TrueSigma <- BOLDC/(TrueD * as.vector(sqrt(design_factor)))

    # Tau: values come from estimateBSvar.R, no covariate
    # 0th, 50th and 100th percentile of observed between-study variability
    Tau <- sqrt(c(0,0.10,0.495))
    I2 <- c(0, 62.61, 87.28)
    # NOTE: tau^2 corresponds to variance of random slopes between-studies.
      # Also add random intercept
    Trueb0_Bstud <- 1

    # Hedges' g can be obtained by multiplying Cohen's d with the correction factor.
    TrueG <- TrueD * NeuRRoStat::corrJ(N = nsub)


    ###################
    #### True spatial location
    ###################

    # True center of activation
    TrueLocations <- c(5,5,5)

    # We generate a temporary design for getting a true signal
    truthdesign <- neuRosim::simprepTemporal(1, 1, onsets = 1, effectsize = 1,
                                   durations = 1, TR = 1, acc = 0.1)
    # Now use this to get a sphere shaped area
    area <- neuRosim::simprepSpatial(regions = 1, coord = list(TrueLocations),
                           radius = ext, form = "sphere", fading = 0)
    truth <- neuRosim::simVOLfmri(design = truthdesign, image = area,
                        dim = DIM, SNR = 1, noise = "none")[,,,1]
    # Unsmoothed ground truth
    GroundTruth <- ifelse(truth > 0, 1, 0)

    # Smooth the GT and put it into the map
    SmGT <- AnalyzeFMRI::GaussSmoothArray(GroundTruth, voxdim = voxdim,
                                          ksize = width, sigma = diag(sigma,3))

    # Create the smoothed ground truth mask (where is the true signal)
    MaskGT <- SmGT
    MaskGT[SmGT == 0] <- 0
    MaskGT[SmGT != 0] <- 1
  }

  # Check whether keyword matches one of the objects
  objCheck <- try(get(keyword), silent = TRUE)
  if(grepl(pattern = 'not found', x = objCheck[1])) {
    stop('Keyword not found in settings (see file trueMCvalues.R)')
  }

  # Return the object
  return(get(keyword))
}
