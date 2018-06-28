

# Parameters
delta <- 0.8
N <- 20


# Function to use in second moment of non-central t-distribution
ratioGamma <- function(N){
  num <- gamma((N - 2)/2)
  denom <- gamma((N - 1) / 2)
  return(num/denom)
}

# Function corresponding to true correction factor, h in paper
TrueCorr <- function(N){
  num <- gamma((N - 1)/2)
  denom <- (sqrt((N - 1)/2) * gamma((N - 2)/2))
  return(num/denom)
}

# Check difference between h and J
NeuRRoStat::corrJ(10) - TrueCorr(10)
NeuRRoStat::corrJ(N) - TrueCorr(N)

# Plug in values of second moment of non-central t-distribution
varA <- ((N-1) * (1 + N*delta**2))/(N - 3)
varB <- ((N*delta**2*(N - 1))/2) * ratioGamma(N)**2
1/N*(varA - varB)

# Compare with version written in paper
PaperVarA <- 1/N*((N**2*delta**2 - N*delta**2+N-1)/(N-3))
PaperVarB <- (delta**2)/(TrueCorr(N)**2)
PaperVarA - PaperVarB

# Compare with Hedges 1981: variance of biased g
VarBiasedG <- (N - 1)/((N - 3)*N)*(1 + (delta**2*N)) - (delta**2/TrueCorr(N)**2)
VarBiasedG

# Now variance of unbiased g
VarBiasedG * (NeuRRoStat::corrJ(N)**2)
VarBiasedG * (TrueCorr(N)**2)

# Variance based on correction J
((N - 1)/((N - 3)*N)*(1 + (delta**2*N)) - (delta**2/NeuRRoStat::corrJ(N)**2)) * (NeuRRoStat::corrJ(N)**2)
# Or written in full
(((NeuRRoStat::corrJ(N)**2) * N**2 * delta**2 - (NeuRRoStat::corrJ(N)**2)*N*delta**2 + (NeuRRoStat::corrJ(N)**2)*N - (NeuRRoStat::corrJ(N)**2)) / (N**2 - 3*N)) - delta**2

# Earlier paper
NeuRRoStat::varHedge(g = delta, N)
# Written in other form
((1 + delta**2*N)/N) - (((delta**2*(N - 3))/2) * ratioGamma(N)**2)

# Difference
VarBiasedG * (TrueCorr(N)**2) - NeuRRoStat::varHedge(g = delta, N)


((N**2*delta**2 - N*delta**2+N-1)/(N**2 - 3*N))


























# OLDER CODE


TrueCorr <- function(N){
  num <- gamma((N - 1)/2)
  denom <- (sqrt((N - 1)/2) * gamma((N - 2)/2))
  return(num/denom)
}
TrueCorr(10)


1 * NeuRRoStat::corrJ(20)


Pa <- 1/N*((N**2*delta**2 - N*delta**2+N-1)/(N-3))
Pb <- (delta**2)/(NeuRRoStat::corrJ(N)**2)
Pb <- (delta**2)/(TrueCorr(N)**2)
Vard <- Pa-Pb

Varg <- (NeuRRoStat::corrJ(N))**2 * Vard

NeuRRoStat::varHedge(g = 0.8 * NeuRRoStat::corrJ(N), N)


(N - 1)/((N - 3)*N)*(1 + (delta**2*N)) - (delta**2/TrueCorr(N)**2)

varA <- ((N-1) * (1 + N*delta**2))/(N - 3)
varB <- ((N*delta**2*(N - 1))/2)*((
  (gamma((N - 2)/2))/(gamma((N - 1)/2)))**2)
1/N*(varA - varB)





