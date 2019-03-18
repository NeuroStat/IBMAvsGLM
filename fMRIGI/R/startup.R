.onAttach <- function(...) {
  messageBlock <- "fMRIGI is the accompanying R package for the fMRI GLM versus MA project."
  packageStartupMessage(messageBlock)
}
