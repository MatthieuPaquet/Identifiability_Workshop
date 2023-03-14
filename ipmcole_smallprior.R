#==========================================================
#==========================================================
# R code for Book Parameter Redundancy and Identifiability
# by Diana J. Cole
# R code is for Immigration Integrated Model
# Section 9.2.2
#==========================================================
#==========================================================

###we converted the Winbugs code in NIMBLE cade
###and added convergence diagnostics and traceplots
### AND used a more reasonable prior for immigration rate: ~Unif(0,1)

library(lattice)
library(coda)
library(nimble)

#------------------------------
# Function for Prior Overlap
#-------------------------------

overlap <- function(data,prior,minv,maxv,freqv,xlabel) {
  # OVERLAP calculates the proportion overlap between the
  # prior and posterior using a kernel density to approximate
  # posterior
  # Also plots a graph of prior and posterior
  # 'data' contains the posterior chain
  # 'prior' contains a vector of prior values evaluated at same interval
  # as 'minv', 'maxv' and 'freqv' values given
  
  k1 <- 0.9 # Controls the smoothness of the kernel
  
  x <- seq(minv,maxv,freqv)
  nn <- length(x)
  fK <- c(rep(0,nn))
  
  overlap<-0
  for (i in 1:nn) {
    fK[i]<-kernel(x[i],data,k1)
    if (fK[i]<prior[i]){
      overlap<-overlap+fK[i]*freqv
    }
    else {
      overlap=overlap+prior[i]*freqv
    }
  }
  
  plot(x,fK,type = "l",ylab="f",xlab=xlabel)
  lines(x,prior,lty=2)  
  return(overlap)
}


kernel <- function(y,data,k1) {
  # KERNEL Calculates a kernel density estimate for a sample in 'data'.
  #   'y' is the value at which we want the kernel estimate.
  #   'k1' can be chosen to vary the amount of smoothing;
  #   Calls the function DELTA.
  
  n <- length(data)
  h <- k1*min(sd(data),IQR(data)/1.34)/n^0.2	
  
  z <- 0
  for (i in 1:n ) {
    z<-z+delta((y-data[i])/h)
  }				            
  z<-z/(n*h);
}

delta <- function(u) {
  # DELTA calculates a normal kernel
  
  y <-(exp(-u*u/2))/sqrt(2*pi);
}

#--------------------------------------------
# NIMBLE Code
#---------------------------------------------


# Population count data
y <- c(14,9,8,15,17,15,13,9,8,9,9,10,12,11,10,11,10,7,5,5,8,12,14,15,18,21)
T <- length(y)
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

r <- rep(NA,2*(T-1))
for(i in 1:(2*(T-1))) {
  r[i] <- sum(mfem[i,])
} 
# Winbugs model:
IPMcode <- nimbleCode(	{
  # priors 
  phij ~ dunif(0,1)
  phia ~ dunif(0,1)
  im~dunif(0,1)
  rho ~ dunif(0,30)
  p ~ dunif(0,1)
  
  # Census
  ##had to change these to avoid slice samplers issues (values are unrealistically high)
  N1[1] <- round(n1)
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
    y[t] ~ dpois(Ntot[t])
  }
  
  
  # CJS 
  q <- 1-p
  for( i in 1:2*(T-1)) {
    m[i,1:T] ~ dmulti(pr[i,1:T],r[i])
  } 	
  # m-array cell probabilities for juveniles
  for(i in 1:(T-1)) {
    pr[i,i]<-phij*p
    for(j in (i+1):(T-1)) {
      pr[i,j] <- phij*phia^(j-i)*q^(j-i)*p
    } 
    for( j in 1:(i-1)) {
      pr[i,j] <- 0
    }
    pr[i,T] <- 1-sum(pr[i,1:(T-1)])
  } 
  # m-array cell probabilities for adults
  for(i in 1:(T-1))            {
    pr[i+T-1,i] <- phia*p
    for(j in (i+1):(T-1))            {
      pr[i+T-1,j] <- phia^(j-i+1)*q^(j-i)*p 
    } 
    for( j in 1:(i-1))                {
      pr[i+T-1,j] <- 0
    } 
    pr[i+T-1,T] <- 1-sum(pr[i+T-1,1:(T-1)])
  } 
}

)

y2 <- y
mK <- mfem

nadsurv <- round(y2[1]/2)
nadim <- 0
n1 <- round(y2[1]/2)
constants <- list(T=T,r=r)#note that for nimble r has to be calculated outside the model
data.or <-list(y=y2,m=mK)
inits <- list(phij=0.4,phia=0.6,rho=2.5,im=0.1,p=0.4,
              n1=n1,nadsurv=nadsurv,nadim=nadim)
#Build the model for the IPM
modelIPM  <-  nimbleModel(IPMcode,
                          constants = constants,data=data.or,inits = inits)
#compile model for IPM
cmodelIPM  <-  compileNimble(modelIPM)
#configure the MCMC
mcmcConf  <-  configureMCMC(cmodelIPM,monitors=c("phij","phia","rho","im","p"))
#Build the MCMC
mcmc  <-  buildMCMC(mcmcConf)
#Compile the  MCMC
cmcmc  <-  compileNimble(mcmc, project = cmodelIPM)
#Run the MCMC
#ADDED traceplots for fecundity and immigration and their posterior correlation
list.samples <- runMCMC(cmcmc,niter=20000,nburnin=10000,thin=1,nchains=3)
#here I voluntarily didn't set the seed so that we can see how r hat varies
#for convergence diagnostic
gelman.diag(list(as.mcmc(list.samples$chain3),as.mcmc(list.samples$chain2),
                 as.mcmc(list.samples$chain1)))
plot(list.samples$chain1[,"im"],type="l",ylim=c(0,1))
lines(list.samples$chain2[,"im"],col="red")
lines(list.samples$chain3[,"im"],col="blue")
plot(list.samples$chain1[,"rho"],type="l",ylim=c(0,30))
lines(list.samples$chain2[,"rho"],col="red")
lines(list.samples$chain3[,"rho"],col="blue")

plot(list.samples$chain1[,"rho"],list.samples$chain1[,"im"],ylim=c(0,1),xlim=c(0,30))
points(list.samples$chain2[,"rho"],list.samples$chain3[,"im"],col=rgb(red=1,green=0,blue=0,alpha=0.5))
points(list.samples$chain3[,"rho"],list.samples$chain3[,"im"],col=rgb(red=0,green=0,blue=1,alpha=0.5))

#change niter and burn in for final versions if necessary
list.samples <- list(runMCMC(cmcmc,niter=20000,nburnin=10000,thin=1,nchains=1,setSeed=T))


results <- do.call(rbind.data.frame, list.samples)

post<-results$rho
minv<-0
maxv<-30
freqv<-0.1
xx <- seq(minv,maxv,freqv)
xlabel<-expression(rho)
prior <-  dunif(xx,0,30)
rhooverlap<-overlap(post,prior,minv,maxv,freqv,xlabel)
rhooverlap

post<-results$im
minv<-0
maxv<-1
freqv<-0.1
xx <- seq(minv,maxv,freqv)
xlabel<-expression(im)
prior <-  dunif(xx,0,1)
imoverlap<-overlap(post,prior,minv,maxv,freqv,xlabel)
imoverlap

post<-results$phia
minv<-0
maxv<-1
freqv<-0.01
xx <- seq(minv,maxv,freqv)
xlabel<-expression(phi[a])
prior <-  dunif(xx,0,1)
phiaoverlap<-overlap(post,prior,minv,maxv,freqv,xlabel)
phiaoverlap

post<-results$phij
minv<-0
maxv<-1
freqv<-0.01
xx <- seq(minv,maxv,freqv)
xlabel<-expression(phi[j])
prior <-  dunif(xx,0,1)
phijoverlap<-overlap(post,prior,minv,maxv,freqv,xlabel)
phijoverlap

post<-results$p
minv<-0
maxv<-1
freqv<-0.01
xx <- seq(minv,maxv,freqv)
xlabel<-expression(p)
prior <-  dunif(xx,0,1)
poverlap<-overlap(post,prior,minv,maxv,freqv,xlabel)
poverlap
