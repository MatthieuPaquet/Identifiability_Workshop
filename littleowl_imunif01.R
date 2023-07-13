#==========================================================
#==========================================================
# R code for Book Parameter Redundancy and Identifiability
# by Diana J. Cole
# R code is for Immigration Integrated Model
# Section 9.2.2
#==========================================================
#==========================================================
setwd("~/Documents/Identifiability_Workshop")
###we converted the Winbugs code in NIMBLE cade
###and added convergence diagnostics and traceplots
### AND used a more reasonable prior for immigration rate: ~Unif(0,1)
### AND added productivity data (REPRODATA <- TRUE/FALSE)
library(lattice)
library(coda)
library(nimble)
library(MCMCvis)
library(truncnorm)
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
  # priors 
  phij ~ dunif(0,1)
  phia ~ dunif(0,1)
  im~dunif(0,1)
  rho ~ dunif(0,30)
  p ~ dunif(0,1)
  # Census
  N1[1] <- round(n1)
  #note that these means are unrealistically high.
  n1 ~ T(dnorm(100,0.0001),0,)  
  NadSurv[1] <- round(nadsurv)
  nadsurv~ T(dnorm(100,0.0001),0,)  	 
  Nadimm[1] <- round(nadim)
  nadim ~ T(dnorm(100,0.0001),0,) 	
  Ntot[1] <- NadSurv[1] + Nadimm[1] + N1[1]  
  
  for(t in 2:T){ 
    mean1[t] <- round(0.5*rho*phij*Ntot[t-1])
    N1[t] ~ dpois(mean1[t])
    mpo[t] <- round(Ntot[t-1]*im)
    NadSurv[t] ~ dbin(phia,Ntot[t-1])
    Nadimm[t] ~ dpois(mpo[t])
    Ntot[t] <- NadSurv[t] + Nadimm[t] + N1[t]
  }
  if (COUNTDATA) {
    # Observation model
  for(t in 1:T){ 
  #In Cole 2020 this loop started from t=2 only but not in Abadi et al.
  y[t] ~ dpois(Ntot[t])
  }#t
  }#COUNTDATA
  if (REPRODATA) {
#productivity
    f ~ dpois(rho * s)
      }
  # CJS 
  if (CRDATA) {
    # Capture-recapture data (CJS model with multinomial likelihood)
  for( i in 1:2*(T-1)) {
    m[i,1:T] ~ dmulti(pr[i,1:T],r[i])
  }#i
  }#CRDATA
  q <- 1-p
  # m-array cell probabilities for juveniles
  for(i in 1:(T-1)) {
    pr[i,i]<-phij*p
    for(j in (i+1):(T-1)) {
      pr[i,j] <- phij*phia^(j-i)*q^(j-i)*p
    } #j
    for( j in 1:(i-1)) {
      pr[i,j] <- 0
    }#j
    pr[i,T] <- 1-sum(pr[i,1:(T-1)])
  }#i
  # m-array cell probabilities for adults
  for(i in 1:(T-1))            {
    pr[i+T-1,i] <- phia*p
    for(j in (i+1):(T-1))            {
      pr[i+T-1,j] <- phia^(j-i+1)*q^(j-i)*p 
    } #j
    for( j in 1:(i-1))                {
      pr[i+T-1,j] <- 0
    } #j
    pr[i+T-1,T] <- 1-sum(pr[i+T-1,1:(T-1)])
  } #i
}

)
# Population count data
y <- c(14,9,8,15,17,15,13,9,8,9,9,10,12,11,10,11,10,7,5,5,8,12,14,15,18,21)
#number of ocasions
nyears=length(y)
# Capture recapture data for females
mfem <- matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,13,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,
                 0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,11,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,12,
                 0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,22,
                 0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,22,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,12,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,11,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,
                 0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,9,
                 0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,15,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,14,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,17,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,17,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,3,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,4,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,11,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,14,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,13,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,15,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,15,
                 1,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,
                 0,3,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
                 0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7,
                 0,0,0,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,
                 0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,
                 0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,
                 0,0,0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,
                 0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,
                 0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,
                 0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,
                 0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,2,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,3,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,2,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,2,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,1,0,0,1,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,4,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,3,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3),nrow=50,ncol=26,byrow=T)
r <- rep(NA,2*(nyears-1))
for(i in 1:(2*(nyears-1))) {
  r[i] <- sum(mfem[i,])
} 
  #number of nestlings
  f <- sum(c(27,19,25,25,47,46,26,29,23,24,20,21,33,32,
             35,35,8,7,17,10,24,31,28,30,33,25)) # number of offspring produced
  s <- sum(c(15,9,8,17,18,16,13,9,8,9,9,11,13,11,
             11,11,10,7,5,5,8,12,14,15,20,24)) # number of breeding females counted
  # Bundle data
  cr.data <- repro.data <- count.data <- list()
  cr.constants <- repro.constants <- count.constants <- list()
  base.constants <- list(T=nyears)
if (CRDATA) {
  cr.data <- list(m=mfem)
  cr.constants <- list(r=r)
}
if (COUNTDATA) {
  count.data <- list(y=y)
}
if (REPRODATA) {
  repro.data <- list(f=f)
  repro.constants <- list(s=s)
}
constants <- c(cr.constants,repro.constants,base.constants)
data <- c(cr.data,count.data,repro.data)

nadsurv <- round(y[1]/2)
nadim <- 1
n1 <- round(y[1]/2)
inits <-  function(){list(phij=runif(1,0.3,0.5),phia=runif(1,0.5,0.7),rho=2.5,im=0.1,p=0.4,
              n1=n1,nadsurv=nadsurv,nadim=nadim,NadSurv=rep(nadsurv,nyears),N1=rep(n1,nyears),Nadimm=rep(nadim,nyears))}
# MCMC settings
ni <- 20000; nb <- 10000; nc <- 3; nt <- 1;
## Pre-sample initial values
#(here they are the same for all chains but I set seed in case we'd change that,
#and the seed is useful later)
set.seed(1)
Inits <- list()
for(i in 1:nc){
  Inits[[i]] <- inits()
}
parameters <- c("phij","phia","rho","im","p","n1","nadsurv","nadim")
#Build and run the model for the IPM
out <- nimbleMCMC(code = IPMcode,
                  data = data,
                  constants = constants,
                  inits = Inits,
                  monitors = parameters,
                  samplesAsCodaMCMC = TRUE,
                  nburnin = nb,
                  niter = ni,
                  thin = nt,
                  nchains = nc,
                  setSeed = T)
## Specify priors for relevant parameters
simNo <- nc*(ni-nb)/nt
priors <- matrix(NA, nrow = simNo, ncol = length(parameters))
priors[ ,1:2] <- runif(simNo, 0, 1)
priors[ ,3] <- runif(simNo, 0, 30)
priors[ ,4:5] <- runif(simNo, 0, 1)
priors[ ,6:8] <- rtruncnorm(simNo, mean=100, a=0, b= Inf, sd=sqrt(1/0.0001))

# Zoom on prior - all parameters
MCMCtrace(out, 
          params = parameters,
          ISB = FALSE,
          exact = TRUE,
          priors = priors,
          pdf = TRUE,
          post_zm = FALSE,
          iter=simNo/nc,
          Rhat = TRUE,
          filename=paste0("plots/littleowl",datascenario,"mcmctrace"))
save(out,file=paste0("data/output_littleowl",datascenario,"mcmc.Rdata"))
