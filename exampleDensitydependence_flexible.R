#code adapted from:
# Schaub & Kéry (2022) Integrated Population Models
# Chapter 8 : Integrated population models with density-dependence

#to be run ith Nimble and implement prior-posterior overlaps
#and include/exclude the different datasets
#here we deleted random time variation of the demographic rates
library(lattice)
library(coda)
library(nimble)
# To load the red-backed shrike data 
library(IPMbook)

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
  # Models for demographic rates
  for (t in 1:(n.occasions-1)){
    logit(phij[t]) <- alpha[1] + beta[1] * (N[t] - mean.C)
    logit(phia[t]) <- alpha[2] + beta[2] * (N[t] - mean.C)
  }
  for (t in 1:n.occasions){
    log(f[t]) <- alpha[3] + beta[3] * (N[t] - mean.C)
  }

  # Priors for variance parameters (hyperparameters)
  tau ~ dgamma(0.001, 0.001)
  sigma2 <- 1 / tau

  # Priors for immigration and resighting rates
  omega ~ dunif(0, 75)
  pj ~ dunif(0, 1)
  pa ~ dunif(0, 1)

  # Priors for regression parameters of the demographic rates
  for (j in 1:3){
    alpha[j] ~ dnorm(0, 0.001)
    beta[j] ~ dunif(-5, 5)
  }
  # Population count data (state-space model)
  # Model for initial stage-spec. population sizes: discrete uniform priors
  R[1] ~ dcat(pNinit[1:50])                             # Local recruits
  S[1] ~ dcat(pNinit[1:50])                             # Surviving adults
  
  for (t in 1:n.occasions){
  I[t] ~ dpois(omega)                          # No. Immigrants
    }
  # Process model over time: our model of population dynamics
  for (t in 2:n.occasions){
    R[t] ~ dpois(f[t-1] / 2 * phij[t-1] * N[t-1]) # No. local recruits
    S[t] ~ dbin(phia[t-1], N[t-1])                # No. surviving adults
  }

  # Observation model
  for (t in 1:n.occasions){
    N[t] <- S[t] + R[t] + I[t]                    # Total number of breeding females
    logN[t] <- log(N[t])
    if (COUNTDATA) {
    C[t] ~ dlnorm(logN[t], tau)
    }#COUNTDATA
  }
  if (CRDATA) {
  # Capture-recapture data (CJS model with multinomial likelihood)
  # Define the multinomial likelihood
  for (t in 1:(n.occasions-1)){
    marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,1:n.occasions], rel.j[t])
    marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,1:n.occasions], rel.a[t])
  }
  }#CRDATA
  # Define the cell probabilities of the m-arrays
  for (t in 1:(n.occasions-1)){
    # Main diagonal
    pr.j[t,t] <- phij[t] * pj
    pr.a[t,t] <- phia[t] * pa
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
      pr.j[t,j] <- phij[t] * prod(phia[(t+1):j]) * (1-pj) * (1-pa)^(j-t-1) * pa
      pr.a[t,j] <- prod(phia[t:j]) * (1-pa)^(j-t) * pa
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

  # Productivity data (Poisson regression model)
  if (REPRODATA) {
  for (t in 1:n.occasions){
    J[t] ~ dpois(B[t] * f[t])
  }
  }#REPRODATA
}
)
data(redbacked)
str(redbacked) # Not shown
# Bundle data
cr.data <- repro.data <- count.data <- list()
cr.constants <- repro.constants <- count.constants <- list()
base.constants <- with(redbacked,list(n.occasions=length(count),mean.C=mean(count),pNinit=dUnif(1, 50)))
if (CRDATA) {
  cr.data <- with(redbacked,list(marr.j=marr.j, marr.a=marr.a))
  cr.constants <- with(redbacked,list(rel.j=rowSums(marr.j), rel.a=rowSums(marr.a)))
}
if (COUNTDATA) {
  count.data <- list(C=redbacked$count)
}
if (REPRODATA) {
  repro.data <- list(J=redbacked$J)
  repro.constants <- list(B=redbacked$B)
}

constants <- c(cr.constants,repro.constants,base.constants)
nim.data <- c(cr.data,count.data,repro.data)
set.seed(1)
# Define initial values
inits <- function() {list(beta=c(0,0,0),alpha=c(qlogis(c(0.05,0.35)),log(3)),pj=runif(1,0.2,0.7), pa=runif(1,0.1,0.6))}
# MCMC settings
ni <- 40000; nb <- 10000; nc <- 3; nt <- 4; 
## Pre-sample initial values
Inits <- list()
for(i in 1:nc){
  Inits[[i]] <- inits()
}
# Define parameters to be monitored
parameters <- c("phij", "phia", "f", "omega",
                "R", "S", "I", "N", "pa", "pj",
                "sigma2","tau", "alpha", "beta")


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
plotParams<- c(paste0("alpha[",1:3,"]"),
               paste0("beta[",1:3,"]"),"R[1]",
               "S[1]","tau","omega","pj","pa")
## Specify priors for relevant parameters
simNo <- nc*(ni-nb)/nt

priors <- matrix(NA, nrow = simNo, ncol = length(plotParams))
priors[ ,1:3] <- rnorm(simNo, 0, sd=sqrt(1/0.001))
priors[ ,4:6] <- runif(simNo, -5, 5)
priors[ ,7:8] <- runif(simNo, 0, 50) #continuous instead of discrete prior
priors[ ,9] <- rgamma(simNo,shape=0.001,rate=0.001)
priors[ ,10] <- runif(simNo,0,75)
priors[ ,11:length(plotParams)] <- runif(simNo, 0, 1)

## Plot prior-posterior overlap using MCMCvis package

# Zoom on prior - all parameters
MCMCtrace(out, 
          params = plotParams,
          ISB = FALSE,
          exact = TRUE,
          priors = priors,
          pdf = TRUE,
          post_zm = FALSE,
          iter=simNo/nc,
          Rhat = T,
          filename=paste0("redbackedshrikeDD",datascenario,"mcmctrace"))
save(out,file=paste0("output_redbackedshrikeDD",datascenario,"mcmc.Rdata"))
