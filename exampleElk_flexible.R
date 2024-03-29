library(IPMbook)
library(nimble)

# SETUP #
#--------#

## Set seed
mySeed <- 0
set.seed(mySeed)

## Define constant variables
Amax <- 17 # Total number of age classes
Tmax <- 6 # Total number of years

UIpriorN.max <- 1000 # Upper bound of non-informative prior on initial population size

## Set prior toggle
#inform.prior <- TRUE # Use informative prior for initial population size
inform.prior <- FALSE # Use non-informative prior for initial population size

## Set data inclusion toggles
# Age-at-harvest data
useData.AaH <- TRUE

# Harvest reporting data
useData.Hreports <- FALSE

# Telemetry data
useData.telemetry <- FALSE

## Set test-run toggle
#testRun <- TRUE # Tests setup with only 10 iterations / chain
testRun <- FALSE # Full run


# DATA PREPARATION #
#------------------#

## Load data
data(elk)

## Set parameter "guesstimates" for making informative priors
s <- 0.8 # annual survival
f <- 0.25 # recruitment
r <- 0.5 # reporting rate
h <- 0.1 # hunting mortality probability

## Create matrix population model (transition matrix A)
A <- matrix(0, ncol = Amax, nrow = Amax)
A[1, 2:Amax] <- f
for (a in 2:Amax){
  A[a, a-1] <- s
}

## Compute stable age distribution (right eigenvector of A)
z <- which.max(Re(eigen(A)$values))
revec <- Re(eigen(A)$vectors[, z])

## Population size in first age class of first year
N1 <- elk$C[1, 1] / (h * r)

## Compute age-specific population sizes in first year
n <- N1 * revec / revec[1]

## Bundle data and constants and produce data overview
nim.data <- with(elk, list(C = C, a = H[1,], b = H[2,], R = R))

nim.constants <- with(elk, list(total = colSums(R), n.years = Tmax,
                                n.age = Amax, 
                                N1min = C[,1], upper = UIpriorN.max))

## Add informative prior information if necessary
if(inform.prior){
  nim.data$n <- n
}


# FLEXIBLE MODEL CODE #
#---------------------#

elk.IPMcode <- nimbleCode({
  
  ## Priors and linear models
  for (t in 1:(n.years-1)){
    f[t] ~ dunif(0, 1)
  }
  for (t in 1:n.years){
    s[t] ~ dunif(0, 1)
    h[t] ~ dunif(0, 1)
    r[t] ~ dunif(0, 1)
  }
  
  ## Population process model
  
  # Initial population size
  for (a in 1:n.age){
    
    if(inform.prior){
      # Informative Poisson priors
      N[a, 1] ~ dpois(n[a])
    }else{
      # Non-informative uniform priors
      #n[a] ~ dunif(0.5, upper+0.4999)    # Ensures that 1 and upper are chosen with the same probability as values 2 to 999
      n[a] ~ dunif(N1min[a], upper)
      N[a, 1] <- round(n[a])
    }
  }
  
  # Process model over time: our model of population dynamics
  for (t in 1:(n.years-1)){
    N[1, t+1] ~ dpois((Ntot[t] - N[1, t]) * f[t])
    for (a in 2:n.age){
      N[a, t+1] ~ dbin((1-h[t]) * s[t], N[a-1, t])
    } #a
  } #t
  
  # Derived quantity: total year-specific population size
  for (t in 1:n.years){
    Ntot[t] <- sum(N[1:n.age, t])
  }
  
  
  ## Age-at-Harvest data likelihood
  if(useData.AaH){
    for (t in 1:n.years){
      for (a in 1:n.age){
        C[a, t] ~ dbin(h[t] * r[t], N[a, t])
      } #a
    } #t
  }

  
  ## Hunter survey data likelihood (logistic regression model)
  if(useData.Hreports){
    for (t in 1:n.years){
      b[t] ~ dbin(r[t], a[t])
    }
  }

  
  ## Radio tracking data likelihood (multinomial model)
  if(useData.telemetry){
    for (t in 1:n.years){
      R[1:3, t] ~ dmulti(prt[1:3, t], total[t])
      prt[1, t] <- h[t]
      prt[2, t] <- (1-h[t]) * (1-s[t])
      prt[3, t] <- (1-h[t]) * s[t]
    }
  }

})


# MODEL SETUP #
#-------------#

## Set MCMC parameters
if(testRun){
  ni <- 10
  nb <- 0
  nc <- 3
  nt <- 1
}else{
  ni <- 250000
  nb <- 50000
  nc <- 3
  nt <- 20
}

## Initial value simulation function
inits <- function(){
  
  InitList <- list(s = runif(nim.constants$n.years, 0.8, 1), 
                   h = runif(nim.constants$n.years, 0.05, 0.15))
  
  if(!inform.prior){
    InitList$n <- round(n)
  }
  
  return(InitList)
}


inits_full <- function(){
  
  # Simulate rates
  s <- runif(nim.constants$n.years, 0.8, 1)
  h <- runif(nim.constants$n.years, 0.05, 0.15)
  f <- runif(nim.constants$n.years, 0.25, 0.5)
  r <- runif(nim.constants$n.years, 0.6, 0.9)
  
  # Set initial population sizes
  n <- round(n)
  N <- matrix(NA, nrow = nim.constants$n.age, ncol = nim.constants$n.years)
  
  N[,1] <- n
  
  # Simulate population sizes
  for (t in 1:(nim.constants$n.years-1)){
    
    N[1, t+1] <- rpois(1, sum(N[2:nim.constants$n.age, t]) * f[t])
    
    for (a in 2:(nim.constants$n.age)){
      
      N[a, t+1] <- rbinom(1, size = N[a-1, t], prob = (1-h[t]) * s[t])
    }
  }
  
  # Check for discrepancies between simulated N and AaH data
  if(any(nim.data$C > N)){
    stop("Simulated population size lower than number reported as harvested.")
  }
  
  # Assemble and return initial values
  InitList <- list(s = s, h = h, f = f, r = r,
                   N = N)
  if(!inform.prior){
    InitList$n <- n
  }
  
  return(InitList)
}


## Parameters monitored
parameters <- c("h", "s", "f", "r", "Ntot", "N")

## Pre-sample initial values
Inits <- list()

for(i in 1:nc){
  if(inform.prior){
    Inits[[i]] <- inits()
  }else{
    Inits[[i]] <- inits_full()
  }
}

# MODEL RUN #
#-----------#

## Run model with NIMBLE
out <- nimbleMCMC(code = elk.IPMcode,
                  data = nim.data,
                  constants = nim.constants,
                  inits = Inits,
                  monitors = parameters,
                  samplesAsCodaMCMC = TRUE,
                  nburnin = nb,
                  niter = ni,
                  thin = nt,
                  nchains = nc,
                  setSeed = mySeed)

# NOTE (if using the reduced initial value simulation from the original example, function inits(), only):
# Nimble will point out that initialization is incomplete and that there
# are NA in some model parameters. The former means that not all
# nodes have been initialized, while the latter means that some lowest-level
# parameters (i.e. those that have priors in the model) do not have initial 
# values. Ideally, we want to at the very least initialize all lowest-level 
# parameters. 

# PRIOR-POSTERIOR OVERLAP #
#-------------------------#

library(MCMCvis)

## List relevant parameters
plotParams1 <- c(paste0("f[", 1:(Tmax-1), "]"),
                paste0("s[", 1:Tmax, "]"),
                paste0("h[", 1:Tmax, "]"),
                paste0("r[", 1:Tmax, "]"))

plotParams2 <- paste0("N[", 1:Amax, ", 1]")

## Specify priors for relevant parameters
simNo <- nc*(ni-nb)/nt

# Fecundity, survival, harvest, and reporting priors
priors1 <- matrix(NA, nrow = simNo, ncol = length(plotParams1))

for(i in 1:length(plotParams1)){
  priors1[ ,i] <- runif(simNo, 0, 1)
}

# Initial population sizes
priors2 <- matrix(NA, nrow = simNo, ncol = length(plotParams2))

for(i in 1:length(plotParams2)){
  
  if(inform.prior){
    priors2[ ,i] <- rpois(simNo, n[i])
  }else{
    priors2[ ,i] <- round(runif(simNo, elk$C[i, 1], UIpriorN.max))
  }
}

# All parameters
priorsAll <- cbind(priors1, priors2)
  
## Plot prior-posterior overlap using MCMCvis package

# All parameters, print to pdf
MCMCtrace(out, 
          params = c(plotParams1, plotParams2),
          ISB = FALSE,
          exact = TRUE,
          priors = priorsAll,
          pdf = TRUE)

# # Zoom on posterior - vital rate parameters
# MCMCtrace(out, 
#           params = plotParams1,
#           ISB = FALSE,
#           exact = TRUE,
#           priors = priors1,
#           pdf = FALSE)
# 
# # Zoom on prior - vital rate parameters
# MCMCtrace(out, 
#           params = plotParams1,
#           ISB = FALSE,
#           exact = TRUE,
#           priors = priors1,
#           pdf = FALSE,
#           post_zm = FALSE)
# 
# # Zoom on posterior - initial population sizes
# MCMCtrace(out, 
#           params = plotParams2,
#           ISB = FALSE,
#           exact = TRUE,
#           priors = priors2,
#           pdf = FALSE)
# 
# # Zoom on prior - initial population sizes
# MCMCtrace(out, 
#           params = plotParams2,
#           ISB = FALSE,
#           exact = TRUE,
#           priors = priors2,
#           pdf = FALSE,
#           post_zm = FALSE)
