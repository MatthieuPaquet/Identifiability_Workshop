#code adapted from:
# Schaub & Kéry (2022) Integrated Population Models
# Chapter 6 : Benefits of integrated population modeling
# ------------------------------------------------------

#to be run ith Nimble and implement prior-posterior overlaps
#and include/exclude the different datasets
library(IPMbook) 
library(lattice)
library(coda)
library(nimble)
#Set which datasets are included
CRDATA <- TRUE
#CRDATA <- FALSE
if (CRDATA) { crname <- "CR_"
} else { crname <- ""}
REPRODATA <- TRUE
#REPRODATA <- FALSE
if (REPRODATA) { reproname <- "P_"
} else {reproname <- ""}
COUNTDATA <- TRUE
#COUNTDATA <- FALSE
if (COUNTDATA) { countname <- "C_"
} else {countname <- ""}
#to name files later
datascenario <- paste(crname,reproname,countname,sep="")
#--------------------------------------------
# NIMBLE Code
#---------------------------------------------
IPMcode <- nimbleCode(	{
  # Priors and linear models
  sj ~ dunif(0, 1)
  sa ~ dunif(0, 1)
  p ~ dunif(0, 1)
  f ~ dunif(0, 10)

  sigma ~ dunif(0.5, 100)
  tau <- pow(sigma, -2)

  # Population count data (state-space model)
  # Model for the initial population size: uniform priors
  N[1,1] ~ dunif(1, 300)
  N[2,1] ~ dunif(1, 300)

  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1)){
    N[1,t+1] <- f / 2 * sj * (N[1,t] + N[2,t])
    N[2,t+1] <- sa * (N[1,t] + N[2,t])
  }
if (COUNTDATA) {
  # Observation model
  for (t in 1:n.occasions){
    C[t] ~ dnorm(N[1,t] + N[2,t], tau)
  }
}#COUNTDATA
  if (REPRODATA) {
  # Productivity data (Poisson regression model)
  nJ ~ dpois(n.rep * f)
  }#REPRODATA
  if (CRDATA) {
  # Capture-recapture data (CJS model with multinomial likelihood)
  # Define the multinomial likelihood
  for (t in 1:(n.occasions-1)){
    marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,1:n.occasions], rel.j[t])
    marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,1:n.occasions], rel.a[t])
  }
  }#CRDATA
  # Define the cell probabilities of the m-arrays
  q <- 1 - p                   # Probability of non-recapture
  for (t in 1:(n.occasions-1)){
    # Main diagonal
    pr.j[t,t] <- sj * p
    pr.a[t,t] <- sa * p
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
      pr.j[t,j] <- sj * sa^(j-t) *q^(j-t)*p
      pr.a[t,j] <- sa^(j-t+1)*q^(j-t)*p 
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
      pr.j[t,j] <- 0
      pr.a[t,j] <- 0
    } #j
  } #t
  # Last column: probability of non-recapture
  for (t in 1:(n.occasions-1)){
    pr.j[t,n.occasions] <- 1-sum(pr.j[t,1:(n.occasions-1)])
    pr.a[t,n.occasions] <- 1-sum(pr.a[t,1:(n.occasions-1)])
  }

  # Derived parameters
  # Annual population growth rate
  for (t in 1:(n.occasions-1)){
    ann.growth.rate[t] <- (N[1,t+1] + N[2,t+1]) / (N[1,t] + N[2,t])
  }
  # Total population size
  for (t in 1:n.occasions){
    Ntot[t] <- N[1,t] + N[2,t]
  }
}
)

data(woodchat6)
str(woodchat6)
# List of 6
# $ ch   : num [1:947, 1:10] 1 1 1 1 1 0 1 1 1 0 ...
# $ age  : num [1:947] 2 2 2 2 2 2 2 2 2 2 ...
# $ count: num [1:10] 110 104 100 85 85 71 118 112 91 104
# $ J    : num [1:10] 147 144 132 131 178 178 235 169 177 186
# $ B    : num [1:10] 48 48 45 37 53 56 74 59 55 60
# $ f    : num [1:535] 4 7 5 5 1 5 4 1 3 0 ...

marr <- marrayAge(woodchat6$ch, woodchat6$age)

# Bundle data
cr.data <- repro.data <- count.data <- list()
cr.constants <- repro.constants <- count.constants <- list()
base.constants <- list(n.occasions=length(woodchat6$count))
if (CRDATA) {
  cr.data <- list(marr.j=marr[,,1], marr.a=marr[,,2])
  cr.constants <- list(rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]))
}
if (COUNTDATA) {
  count.data <- list(C=woodchat6$count)
}
if (REPRODATA) {
  repro.data <- list(nJ=sum(woodchat6$J))
  repro.constants <- list(n.rep=sum(woodchat6$B))
}

constants <- c(cr.constants,repro.constants,base.constants)
nim.data <- c(cr.data,count.data,repro.data)

# Initial values
inits <- function(){list(sj=runif(1, 0, 0.5))}
set.seed(1)
# MCMC settings
ni <- 40000; nb <- 10000; nc <- 3; nt <- 3;
## Pre-sample initial values
Inits <- list()
for(i in 1:nc){
  Inits[[i]] <- inits()
}
# Parameters monitored
parameters <- c("sj", "sa", "p", "f", "sigma", "N", "ann.growth.rate", "Ntot")



out <- nimbleMCMC(code = IPMcode,
                  data = nim.data,
                  constants = constants,
                  inits = Inits,
                  monitors = parameters,
                  samplesAsCodaMCMC = TRUE,
                  nburnin = nb,
                  niter = ni,
                  thin = nt,
                  nchains = nc,
                  setSeed = T)

# PRIOR-POSTERIOR OVERLAP #
#-------------------------#

library(MCMCvis)

## List relevant parameters
plotParams1 <- c("f",
                 "p","sj","sa")
plotParams2 <- c("sigma",paste0("N[", 1:2, ", 1]"))

## Specify priors for relevant parameters
simNo <- nc*(ni-nb)/nt

# Fecundity, survival, harvest, and reporting priors
priors1 <- matrix(NA, nrow = simNo, ncol = length(plotParams1))
priors1[ ,1] <- runif(simNo, 0, 10)
priors1[ ,2:length(plotParams1)] <- runif(simNo, 0, 1)
# Initial population sizes
priors2 <- matrix(NA, nrow = simNo, ncol = length(plotParams2))
priors2[ ,1] <- runif(simNo, 0.5, 100)
priors2[ ,2:length(plotParams2)]  <- runif(simNo, 1, 300)
## Plot prior-posterior overlap using MCMCvis package

# Zoom on prior - vital rate parameters
MCMCtrace(out, 
          params = plotParams1,
          ISB = FALSE,
          exact = TRUE,
          priors = priors1,
          pdf = TRUE,
          post_zm = FALSE,
          iter=simNo/nc,
          Rhat = T,
          filename=paste0("woodchatshrike",datascenario,"mcmctraceparam1"))
MCMCtrace(out, 
          params = plotParams2,
          ISB = FALSE,
          exact = TRUE,
          priors = priors2,
          pdf = TRUE,
          post_zm = FALSE,
          iter=simNo/nc,
          Rhat = T,
          filename=paste0("woodchatshrike",datascenario,"mcmctraceparam2"))
