# Schaub & Kéry (2022) Integrated Population Models
# Chapter 17 : Elk
# ----------------

# Run time for test script xx mins, full run xx mins

# 17.4 Component data likelihoods
# =============================================

library(IPMbook); library(nimble)
data(elk)
str(elk)
# List of 3
# $ C: num [1:17, 1:6] 21 36 15 8 7 8 6 5 3 3 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:17] "1" "2" "3" "4" ...
# .. ..$ : chr [1:6] "1988" "1989" "1990" "1991" ...
# $ H: num [1:2, 1:6] 275 143 290 154 211 211 360 272 201 201 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:2] "Hunters surveyed" "Registered harvest"
# .. ..$ : chr [1:6] "1988" "1989" "1990" "1991" ...
# $ R: num [1:3, 1:6] 1 2 10 4 1 20 0 1 24 0 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:3] "hunted" "died naturally" "survived"
# .. ..$ : chr [1:6] "1988" "1989" "1990" "1991" ...

# Set seed
mySeed <- 0
set.seed(mySeed)

# 17.4.1 Age-at-harvest data (no code)
# 17.4.2 Hunter survey data (no code)
# 17.4.3 Radio tracking data (no code)


# 17.5 The integrated population model
# ====================================

s <- 0.8                                          # Guess of annual survival
f <- 0.25                                         # Guess of recruitment

# Create matrix population model (transition matrix A)
A <- matrix(0, ncol=17, nrow=17)
A[1,2:17] <- f
for (a in 2:17){
  A[a,a-1] <- s
}

# Compute stable age distribution (right eigenvector of A)
z <- which.max(Re(eigen(A)$values))
revec <- Re(eigen(A)$vectors[,z])

# Population size in first age class of first year
r <- 0.5                                          # Guess of reporting rate
h <- 0.1                                          # Guess of hunting mortality
N1 <- elk$C[1,1] / (h * r)

# Compute age-specific population sizes in first year
n <- N1 * revec / revec[1]

# Bundle data and constants and produce data overview
nim.data <- with(elk, list(C=C, a=H[1,], b=H[2,], R=R, n=n))
nim.constants <- with(elk, list(total=colSums(R), n.years=ncol(C),
                           n.age=nrow(C)))

str(nim.data)
str(nim.constants)


# Write NIMBLE model
elk.IPMcode1 <- nimbleCode({
  
  # Priors and linear models
  for (t in 1:(n.years-1)){
    f[t] ~ dunif(0, 1)
  }
  for (t in 1:n.years){
    s[t] ~ dunif(0, 1)
    h[t] ~ dunif(0, 1)
    r[t] ~ dunif(0, 1)
  }

  # Age-at-harvest data (state-space model)
  # Model for the initial population size: Poisson priors
  for (a in 1:n.age){
    N[a,1] ~ dpois(n[a])
  }

  # Process model over time: our model of population dynamics
  for (t in 1:(n.years-1)){
    N[1,t+1] ~ dpois((Ntot[t] - N[1,t]) * f[t])
    for (a in 2:n.age){
      N[a,t+1] ~ dbin((1-h[t]) * s[t], N[a-1,t])
    } #a
  } #t

  # Derived quantity: total year-specific population size
  for (t in 1:n.years){
    Ntot[t] <- sum(N[1:n.age,t])
  }

  # Observation model
  for (t in 1:n.years){
    for (a in 1:n.age){
      C[a,t] ~ dbin(h[t] * r[t], N[a,t])
    } #a
  } #t

  # Hunter survey data (logistic regression model)
  for (t in 1:n.years){
    b[t] ~ dbin(r[t], a[t])
  }

  # Radio tracking data (multinomial model)
  for (t in 1:n.years){
    R[1:3,t] ~ dmulti(prt[1:3,t], total[t])
    prt[1,t] <- h[t]
    prt[2,t] <- (1-h[t]) * (1-s[t])
    prt[3,t] <- (1-h[t]) * s[t]
  }
})

# Initial value simulation function
inits <- function() {list(s=runif(nim.constants$n.years, 0.8, 1), h=runif(nim.constants$n.years,
                                                                      0.05, 0.15))}

# Parameters monitored
parameters <- c("h", "s", "f", "r", "Ntot", "N")

# MCMC settings
#ni <- 10; nb <- 0; nc <- 3; nt <- 1
ni <- 250000; nb <- 50000; nc <- 3; nt <- 20

# Pre-sample initial values
Inits <- list()
for(i in 1:nc){
  Inits[[i]] <- inits()
}

# Run NIMBLE (ART 8 min) and check convergence
out1 <- nimbleMCMC(code = elk.IPMcode1,
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

# NOTE: Nimble will point out that initialization is incomplete and that there
# are NA in initial values for some parameters. The former means that not all
# nodes have been initialized, while the latter means that some lowest-level
# parameters (i.e. those that have priors in the model) do not have initial 
# values. Ideally, we want to at the very least initialize all lowest-level 
# parameters. 

plot(out1, ask = T)


# 17.6 Results on elk population dynamics
# =======================================

print(summary(out1))

# Re-format output as matrix
out1.mat <- as.matrix(out1)

# Calculation of annual survival (s*)
s.star <- (1-out1.mat[, paste0("h[", 1:6, "]")]) * out1.mat[, paste0("s[", 1:6, "]")]

# ~~~~ Code for Fig. 17.3 ~~~~
library(scales)
cl <- viridis_pal(option='E')(20)[c(18,2)]
qu <- function(x) quantile(x, c(0.025, 0.975))
op <- par(mfrow=c(2,1), mar=c(3.5, 5, 1, 1))
d <- 0.1
plot(y=colMeans(out1.mat[, paste0("h[", 1:6, "]")]), x=(1:6)+d, ylim=c(0,0.4), xlim=c(1-d, 6+d), type="b",
     pch=16, axes=FALSE, ylab="Probability", xlab=NA, col=cl[2])
segments((1:6)+d, matrixStats::colQuantiles(out1.mat[, paste0("h[", 1:6, "]")], probs = 0.025), (1:6)+d, matrixStats::colQuantiles(out1.mat[, paste0("h[", 1:6, "]")], probs = 0.975), col=cl[2])
axis(1, at=1:6, labels = 1988:1993)
axis(2, las=1)

points(y=1-colMeans(out1.mat[, paste0("s[", 1:6, "]")]), x=(1:6)-d, type="b", pch=16, col=cl[1])
segments((1:6)-d, 1-matrixStats::colQuantiles(out1.mat[, paste0("s[", 1:6, "]")], probs = 0.025), (1:6)-d, 1-matrixStats::colQuantiles(out1.mat[, paste0("s[", 1:6, "]")], probs = 0.975), col=cl[1])
legend('topleft', pch=rep(16,2), col=rev(cl),
       legend=c('Hunting mortality', 'Background mortality'), bty='n')

plot(colMeans(out1.mat[, paste0("f[", 1:5, "]")]), x=2:6, ylim=c(0,1.2), xlim=c(1-d, 6+d), type="b",
     pch=16, axes=FALSE, ylab="Probability", xlab=NA, col=cl[1])
segments(2:6, matrixStats::colQuantiles(out1.mat[, paste0("f[", 1:5, "]")], probs = 0.025), 2:6, matrixStats::colQuantiles(out1.mat[, paste0("f[", 1:5, "]")], probs = 0.975), col=cl[1])
axis(1, at=1:6, labels = 1988:1993)
axis(2, las=1, at=c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c(0, 0.2, 0.4, 0.6, 0.8, 1))
points(y=apply(s.star, 2, mean), x=1:6, type="b", pch=16, col=cl[2])
segments(1:6, apply(s.star, 2, qu)[1,], 1:6, apply(s.star, 2, qu)[2,], col=cl[2])
legend('topleft', pch=c(16,16), col=rev(cl),
       legend=c('Annual survival', 'Recruitment'), bty='n')
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ Fig. 17.4 ~~~~
library(scales)
cl <- viridis_pal(option='E')(20)[c(18,11)]
qu <- function(x) quantile(x, c(0.025, 0.975))

# Calculate the number of harvested elk
H <- out1.mat[, paste0("Ntot[", 1:6, "]")] * out1.mat[, paste0("h[", 1:6, "]")]

op <- par(mar=c(3.5, 5, 1, 1), las=1)
z <- cbind(out1.mat[, "Ntot[1]"], H[,1], out1.mat[, "Ntot[2]"], H[,2],  # ~~~ FIXME
           out1.mat[, "Ntot[3]"], H[,3], out1.mat[, "Ntot[4]"], H[,4],
           out1.mat[, "Ntot[5]"], H[,5], out1.mat[, "Ntot[6]"], H[,6])

a <- barplot(apply(z, 2, mean), space = rep(c(0.5,0.1), 6), ylim = c(0, 3000),
             ylab="Numbers", col=rep(cl,6), border=NA)
axis(1, at = c(mean(a[1:2]), mean(a[3:4]), mean(a[5:6]), mean(a[7:8]),
               mean(a[9:10]), mean(a[11:12])), labels = 1988:1993, lwd=0)
segments(a, apply(z, 2, qu)[1,], a, apply(z, 2, qu)[2,], col='black')
legend('topleft', pch=c(15,15), legend=c('Pre-hunting population', 'Harvested'),
       bty='n', col=cl, pt.cex=1.75)
segments(a[1]-0.55, 0, a[12]+0.55, 0)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 17.7 Prior sensitivity analysis
# ===============================

# ~~~~ code for analysis with different priors ~~~~
# Prior Poisson 2 (P2): h = 0.075
s <- 0.8   # guess of annual survival
f <- 0.25  # guess of recruitment

# Create matrix population model
A <- matrix(0, ncol=17, nrow=17)
A[1,2:17] <- f
for (a in 2:17){
  A[a,a-1] <- s
}

# Compute right eigenvector
z <- which.max(Re(eigen(A)$values))
revec <- Re(eigen(A)$vectors[,z])

# Population size in first age class of first year
r <- 0.5     # guess of reporting rate
h <- 0.075   # guess of hunting mortality
N1 <- elk$C[1,1] / (h * r)

# Compute age-dependent population size in first year
n <- N1 * revec / revec[1]

# Bundle data and constants and produce data overview
nim.data <- with(elk, list(C=C, a=H[1,], b=H[2,], R=R, n=n))
nim.constants <- with(elk, list(total=colSums(R), n.years=ncol(C),
                                n.age=nrow(C)))

# Initial value simulation function
inits <- function() {list(s=runif(nim.constants$n.years, 0.8, 1), h=runif(nim.constants$n.years,
                                                                          0.05, 0.15))}

# Parameters monitored
parameters <- c("h", "s", "f", "r", "Ntot", "N")

# MCMC settings
#ni <- 10; nb <- 0; nc <- 3; nt <- 1
ni <- 250000; nb <- 50000; nc <- 3; nt <- 20

# Pre-sample initial values
Inits <- list()
for(i in 1:nc){
  Inits[[i]] <- inits()
}

# Run NIMBLE (ART 8 min) and check convergence
out2 <- nimbleMCMC(code = elk.IPMcode1,
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


# Prior Poisson 3 (P3): h = 0.125
s <- 0.8   # guess of annual survival
f <- 0.25  # guess of recruitment

# Create matrix population model
A <- matrix(0, ncol=17, nrow=17)
A[1,2:17] <- f
for (a in 2:17){
  A[a,a-1] <- s
}

# Compute right eigenvector
z <- which.max(Re(eigen(A)$values))
revec <- Re(eigen(A)$vectors[,z])

# Population size in first age class of first year
r <- 0.5     # guess of reporting rate
h <- 0.125   # guess of hunting mortality
N1 <- elk$C[1,1] / (h * r)

# Compute age-dependent population size in first year
n <- N1 * revec / revec[1]

# Bundle data and constants and produce data overview
nim.data <- with(elk, list(C=C, a=H[1,], b=H[2,], R=R, n=n))
nim.constants <- with(elk, list(total=colSums(R), n.years=ncol(C),
                                n.age=nrow(C)))

# Initial value simulation function
inits <- function() {list(s=runif(nim.constants$n.years, 0.8, 1), h=runif(nim.constants$n.years,
                                                                          0.05, 0.15))}

# Parameters monitored
parameters <- c("h", "s", "f", "r", "Ntot", "N")

# MCMC settings
#ni <- 10; nb <- 0; nc <- 3; nt <- 1
ni <- 250000; nb <- 50000; nc <- 3; nt <- 20

# Pre-sample initial values
Inits <- list()
for(i in 1:nc){
  Inits[[i]] <- inits()
}

# Run NIMBLE (ART 8 min) and check convergence
out3 <- nimbleMCMC(code = elk.IPMcode1,
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


# Discrete uniform prior (U)
s <- 0.8   # guess of annual survival
f <- 0.25  # guess of recruitment

# Create matrix population model
A <- matrix(0, ncol=17, nrow=17)
A[1,2:17] <- f
for (a in 2:17){
  A[a,a-1] <- s
}

# Compute right eigenvector
z <- which.max(Re(eigen(A)$values))
revec <- Re(eigen(A)$vectors[,z])

# Population size in first age class of first year
r <- 0.5     # guess of reporting rate
h <- 0.1     # guess of hunting mortality
N1 <- elk$C[1,1] / (h * r)

# Compute age-dependent population size in first year
n <- N1 * revec / revec[1]
lo <- 0.5 * n
up <- 1.5 * n

# Create a matrix with the values for the categorical distribution
pinit <- dUnif(lower=lo, upper=up)

# Bundle data and constants and produce data overview
nim.data <- with(elk, list(C=C, a=H[1,], b=H[2,], R=R, p=pinit))
nim.constants <- with(elk, list(total=colSums(R), n.years=ncol(C),
                                n.age=nrow(C), up=round(up)))


# Write NIMBLE model
elk.IPMcode2 <- nimbleCode({
  
  # Priors and linear models
  for (t in 1:(n.years-1)){
    f[t] ~ dunif(0, 1)
  }
  for (t in 1:n.years){
    s[t] ~ dunif(0, 1)
    h[t] ~ dunif(0, 1)
    r[t] ~ dunif(0, 1)
  }
  
  # Age-at-harvest data (state-space model)
  # Model for the initial population size: discrete uniform priors
  for (a in 1:n.age){
    N[a,1] ~ dcat(p[a, 1:up[a]])
  }
  
  # Process model over time: our model of population dynamics
  for (t in 1:(n.years-1)){
    N[1,t+1] ~ dpois((Ntot[t] - N[1,t]) * f[t])
    for (a in 2:n.age){
      N[a,t+1] ~ dbin((1-h[t]) * s[t], N[a-1,t])
    } #a
  } #t
  
  # Derived quantity: total year-specific population size
  for (t in 1:n.years){
    Ntot[t] <- sum(N[1:n.age,t])
  }
  
  # Observation model
  for (t in 1:n.years){
    for (a in 1:n.age){
      C[a,t] ~ dbin(h[t] * r[t], N[a,t])
    } #a
  } #t
  
  # Hunter survey data (logistic regression model)
  for (t in 1:n.years){
    b[t] ~ dbin(r[t], a[t])
  }
  
  # Radio tracking data (multinomial model)
  for (t in 1:n.years){
    R[1:3,t] ~ dmulti(prt[1:3,t], total[t])
    prt[1,t] <- h[t]
    prt[2,t] <- (1-h[t]) * (1-s[t])
    prt[3,t] <- (1-h[t]) * s[t]
  }
})

# Initial value simulation function
inits <- function() {list(s=runif(nim.constants$n.years, 0.8, 1), h=runif(nim.constants$n.years,
                                                                          0.05, 0.15))}

# Parameters monitored
parameters <- c("h", "s", "f", "r", "Ntot", "N")

# MCMC settings
#ni <- 10; nb <- 0; nc <- 3; nt <- 1
ni <- 250000; nb <- 50000; nc <- 3; nt <- 20

# Pre-sample initial values
Inits <- list()
for(i in 1:nc){
  Inits[[i]] <- inits()
}

# Run NIMBLE (ART x min) and check convergence
out4 <- nimbleMCMC(code = elk.IPMcode2,
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
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ Fig. 17.5 ~~~~
# library(scales)
# cl <- viridis_pal(option='E')(20)[c(18,11,5,1)]
# op <- par(mfrow=c(2,1), las=1, mar=c(3,4,1,1))
# d <- 0.1
# plot(x=(1:6)-1.5*d, y=out1$mean$Ntot, ylim = c(0, 4000), xlim=c(1-1.5*d, 6+1.5*d), pch=16,
#      axes=FALSE, xlab=NA, ylab=expression(paste("Population size (", italic(N)[tot],")")), col=cl[1])
# segments((1:6)-1.5*d, out1$q2.5$Ntot, (1:6)-1.5*d, out1$q97.5$Ntot, col=cl[1])
# points(x=(1:6)-0.5*d, y=out2$mean$Ntot, pch=16, col=cl[2])
# segments((1:6)-0.5*d, out2$q2.5$Ntot, (1:6)-0.5*d, out2$q97.5$Ntot, col=cl[2])
# points(x=(1:6)+0.5*d, y=out3$mean$Ntot, pch=16, col=cl[3])
# segments((1:6)+0.5*d, out3$q2.5$Ntot, (1:6)+0.5*d, out3$q97.5$Ntot, col=cl[3])
# points(x=(1:6)+1.5*d, y=out4$mean$Ntot, pch=16, col=cl[4])
# segments((1:6)+1.5*d, out4$q2.5$Ntot, (1:6)+1.5*d, out4$q97.5$Ntot, col=cl[4])
# axis(1, at=1:6, labels=1988:1993)
# axis(2)
# 
# plot(x=(1:6)-1.5*d, y=out1$mean$h, ylim = c(0, 0.31), xlim=c(1-1.5*d, 6+1.5*d), pch=16,
#      axes=FALSE, xlab=NA, ylab=expression(paste("Hunting mortality (", italic(h), ")")), col=cl[1])
# segments((1:6)-1.5*d, out1$q2.5$h, (1:6)-1.5*d, out1$q97.5$h, col=cl[1])
# points(x=(1:6)-0.5*d, y=out2$mean$h, pch=16, col=cl[2])
# segments((1:6)-0.5*d, out2$q2.5$h, (1:6)-0.5*d, out2$q97.5$h, col=cl[2])
# points(x=(1:6)+0.5*d, y=out3$mean$h, pch=16, col=cl[3])
# segments((1:6)+0.5*d, out3$q2.5$h, (1:6)+0.5*d, out3$q97.5$h, col=cl[3])
# points(x=(1:6)+1.5*d, y=out4$mean$h, pch=16, col=cl[4])
# segments((1:6)+1.5*d, out4$q2.5$h, (1:6)+1.5*d, out4$q97.5$h, col=cl[4])
# axis(1, at=1:6, labels=1988:1993)
# axis(2)
# legend(x=0.75, y=0.33, pch=rep(16,2), col=cl[1:2],
#        legend=c(expression('P'[1]),expression('P'[2])), bty='n')
# legend(x=1.5, y=0.33, , inset = 0.1, pch=rep(16,2), col=cl[3:4],
#        legend=c(expression('P'[3]), 'U'), bty='n')
# par(op)
# 
# # Fig. 17.6
# library(scales)
# cl <- viridis_pal(option='E')(20)[c(18,11,5,1)]
# qu <- function(x) quantile(x, c(0.025, 0.975))
# 
# # Calculate annual population growth rates
# lam1 <- lam2 <- lam3 <- lam4 <-
#   matrix(NA, ncol=ncol(out1$sims.list$Ntot)-1, nrow=nrow(out1$sims.list$Ntot))
# for (t in 1:5){
#   lam1[,t] <- out1$sims.list$Ntot[,t+1] / out1$sims.list$Ntot[,t]
#   lam2[,t] <- out2$sims.list$Ntot[,t+1] / out2$sims.list$Ntot[,t]
#   lam3[,t] <- out3$sims.list$Ntot[,t+1] / out3$sims.list$Ntot[,t]
#   lam4[,t] <- out4$sims.list$Ntot[,t+1] / out4$sims.list$Ntot[,t]
# }
# 
# op <- par(las = 1, mar = c(3,4,1,1))
# d <- 0.1
# xa <- c(1.5, 2.5, 3.5, 4.5, 5.5)
# plot(x=xa-1.5*d, y=apply(lam1, 2, mean), type='b', pch=16, xlim=c(1,6), ylim=c(0.6, 1.2),
#      col=cl[1], axes=FALSE, xlab=NA, ylab='Population growth rate')
# segments(xa-1.5*d, apply(lam1, 2, qu)[1,], xa-1.5*d, apply(lam1, 2, qu)[2,], col=cl[1])
# points(x=xa-0.5*d, y=apply(lam2, 2, mean), type='b', pch=16, col=cl[2])
# segments(xa-0.5*d, apply(lam2, 2, qu)[1,], xa-0.5*d, apply(lam2, 2, qu)[2,], col=cl[2])
# points(x=xa+0.5*d, y=apply(lam3, 2, mean), type='b', pch=16, col=cl[3])
# segments(xa+0.5*d, apply(lam3, 2, qu)[1,], xa+0.5*d, apply(lam3, 2, qu)[2,], col=cl[3])
# points(x=xa+1.5*d, y=apply(lam4, 2, mean), type='b', pch=16, col=cl[4])
# segments(xa+1.5*d, apply(lam4, 2, qu)[1,], xa+1.5*d, apply(lam4, 2, qu)[2,], col=cl[4])
# axis(2)
# axis(1, at=1:6, labels=1988:1993)
# legend(x=2.5, y=0.75, pch=rep(16,2), col=cl[1:2],
#        legend=c(expression('P'[1]),expression('P'[2])), bty='n')
# legend(x=3.5, y=0.75, pch=rep(16,2), col=cl[3:4], legend=c(expression('P'[3]), 'U'), bty='n')
# par(op)


# Define an IPM that uses a uniform prior for the initial population size. In principle this is the same model as model2.txt, however, instead of using the categorical distribution in JAGS to create a discrete uniform prior, we use here the uniform distribution and round the simulated number to the nearest integer. This is computationally more efficient, and is expected to have no impact on the results.

# Write NIMBLE model
elk.IPMcode3 <- nimbleCode({
  
  # Priors and linear models
  for (t in 1:(n.years-1)){
    f[t] ~ dunif(0, 1)
  }
  for (t in 1:n.years){
    s[t] ~ dunif(0, 1)
    h[t] ~ dunif(0, 1)
    r[t] ~ dunif(0, 1)
  }
  
  # Age-at-harvest data (state-space model)
  # Model for the initial population size: rounded uniform priors
  for (a in 1:n.age){
    n[a] ~ dunif(0.5, upper+0.4999)    # Ensures that 1 and upper are chosen with the same probability as values 2 to 999
    N[a,1] <- round(n[a])
  }
  
  # Process model over time: our model of population dynamics
  for (t in 1:(n.years-1)){
    N[1,t+1] ~ dpois((Ntot[t] - N[1,t]) * f[t])
    for (a in 2:n.age){
      N[a,t+1] ~ dbin((1-h[t]) * s[t], N[a-1,t])
    } #a
  } #t
  
  # Derived quantity: total year-specific population size
  for (t in 1:n.years){
    Ntot[t] <- sum(N[1:n.age,t])
  }
  
  # Observation model
  for (t in 1:n.years){
    for (a in 1:n.age){
      C[a,t] ~ dbin(h[t] * r[t], N[a,t])
    } #a
  } #t
  
  # Hunter survey data (logistic regression model)
  for (t in 1:n.years){
    b[t] ~ dbin(r[t], a[t])
  }
  
  # Radio tracking data (multinomial model)
  for (t in 1:n.years){
    R[1:3,t] ~ dmulti(prt[1:3,t], total[t])
    prt[1,t] <- h[t]
    prt[2,t] <- (1-h[t]) * (1-s[t])
    prt[3,t] <- (1-h[t]) * s[t]
  }
})

# Fit models

# 1. Original data set; uniform prior U(1, 1000)
nim.data <- with(elk, list(C=C, a=H[1,], b=H[2,], R=R))
nim.constants <- with(elk, list(total=colSums(R), n.years=ncol(C),
                                n.age=nrow(C), upper=1000))

# Initial values
# We need good initial values for the population size in the first year
s <- 0.8   # guess of annual survival
f <- 0.25  # guess of recruitment

# Create matrix population model
A <- matrix(0, ncol=17, nrow=17)
A[1,2:17] <- f
for (a in 2:17){
  A[a,a-1] <- s
}
# Compute stable age distribution (right eigenvector)
z <- which.max(Re(eigen(A)$values))
revec <- Re(eigen(A)$vectors[,z])
# Population size in first age class of first year
r <- 0.5   # guess of reporting rate
h <- 0.1   # guess of hunting mortality
N1 <- elk$C[1,1] / (h * r)
# Compute age-specific population sizes in first year
n <- N1 * revec / revec[1]

# Initial value simulation function
inits <- function() {list(s=runif(nim.constants$n.years, 0.8, 1), n=round(n),
                          h=runif(nim.constants$n.years, 0.05, 0.15))}

# Parameters monitored
parameters <- c("h", "s", "f", "r", "Ntot", "N")

# MCMC settings
#ni <- 10; nb <- 0; nc <- 3; nt <- 1
ni <- 250000; nb <- 50000; nc <- 3; nt <- 20

# Pre-sample initial values
Inits <- list()
for(i in 1:nc){
  Inits[[i]] <- inits()
}

# Run NIMBLE (ART x min) and check convergence
out5 <- nimbleMCMC(code = elk.IPMcode3,
                   data = nim.data,
                   constants = nim.constants,
                   inits = Inits,
                   monitors = parameters,
                   samplesAsCodaMCMC = TRUE,
                   nburnin = nb,
                   niter = ni,
                   thin = nt,
                   nchains = nc,
                   setSeed = mySeed + 1)

# NOTE: Here, I added + 1 to the seed since it otherwise picked a "bad" starting
# value for some of the N nodes. This is not a good solution though. It's a 
# "quick" fix and what should be done instead is simulating reasonable initial 
# values for more - if not all - nodes in the model. 
# more

# 2. Original data set; uniform prior U(1, 3000)
nim.data <- with(elk, list(C=C, a=H[1,], b=H[2,], R=R))
nim.constants <- with(elk, list(total=colSums(R), n.years=ncol(C),
                                n.age=nrow(C), upper=3000))

# Initial value simulation function
inits <- function() {list(s=runif(nim.constants$n.years, 0.8, 1), n=round(n),
                          h=runif(nim.constants$n.years, 0.05, 0.15))}

# Parameters monitored
parameters <- c("h", "s", "f", "r", "Ntot", "N")

# MCMC settings
#ni <- 10; nb <- 0; nc <- 3; nt <- 1
ni <- 160000; nt <- 15; nb <- 10000; nc <- 3

# Pre-sample initial values
Inits <- list()
for(i in 1:nc){
  Inits[[i]] <- inits()
}

# Run NIMBLE (ART x min) and check convergence
out6 <- nimbleMCMC(code = elk.IPMcode3,
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

# 3. 20 times more radio tracking data; uniform prior U(1, 1000)
nim.data <- with(elk, list(C=C, a=H[1,], b=H[2,], R=R*20))
nim.constants <- with(elk, list(total=colSums(R*20), n.years=ncol(C),
                                n.age=nrow(C), upper=1000))

# Initial value simulation function
inits <- function() {list(s=runif(nim.constants$n.years, 0.8, 1), n=round(n),
                          h=runif(nim.constants$n.years, 0.05, 0.15))}

# Parameters monitored
parameters <- c("h", "s", "f", "r", "Ntot", "N")

# MCMC settings
#ni <- 10; nb <- 0; nc <- 3; nt <- 1
ni <- 250000; nb <- 50000; nc <- 3; nt <- 20

# Pre-sample initial values
Inits <- list()
for(i in 1:nc){
  Inits[[i]] <- inits()
}

# Run NIMBLE (ART x min) and check convergence
out7 <- nimbleMCMC(code = elk.IPMcode3,
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

# 4. 20 times more radio tracking data; uniform prior U(1, 3000)
nim.data <- with(elk, list(C=C, a=H[1,], b=H[2,], R=R*20))
nim.constants <- with(elk, list(total=colSums(R*20), n.years=ncol(C),
                                n.age=nrow(C), upper=3000))

# Initial value simulation function
inits <- function() {list(s=runif(nim.constants$n.years, 0.8, 1), n=round(n),
                          h=runif(nim.constants$n.years, 0.05, 0.15))}

# Parameters monitored
parameters <- c("h", "s", "f", "r", "Ntot", "N")

# MCMC settings
#ni <- 10; nb <- 0; nc <- 3; nt <- 1
ni <- 250000; nb <- 50000; nc <- 3; nt <- 20

# Pre-sample initial values
Inits <- list()
for(i in 1:nc){
  Inits[[i]] <- inits()
}

# Run NIMBLE (ART x min) and check convergence
out8 <- nimbleMCMC(code = elk.IPMcode3,
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


save(out1, out2, out3, out4, out5, out6, out7, out8, file="ElkResults.Rdata")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~ Fig. 17.7 ~~~~
# load('ElkResults.Rdata')
# library(scales)
# cl <- c(viridis_pal(option='E')(20)[c(18,8)], 'red')
# 
# op <- par("mfrow", "mar")
# layout(matrix(1:16, 4, 4, byrow=TRUE), widths=c(1.1, 1, 1, 1), heights=c(1.1, 1, 1, 1), TRUE)
# 
# up <- 1000
# upy <- 0.006
# 
# par(mar=c(3,3,2,1))
# hist(out5$sims.list$N[,1,1], col=cl[1], border=cl[1], freq=FALSE, xlim=c(0, up),
#      ylim=c(0, upy), main=NA, ylab=NA, xlab=NA, axes=FALSE, breaks=seq(0,1000, by=50))
# abline(h=1/up, col=cl[3])
# text(x=0, y=0.9*upy, labels=bquote(bar(x) == .(round(out5$mean$N[1,1]))), pos=4)
# text(x=0, y=0.8*upy, labels=paste('SD: ', round(out5$sd$N[1,1])), pos=4)
# axis(1)
# mtext('U(1,1000)', side=2, line=1)
# mtext('1y old', side=3, line=0.5)
# 
# par(mar=c(3,2,2,1))
# hist(out5$sims.list$N[,2,1], col=cl[1], border=cl[1], freq=FALSE, xlim=c(0, up),
#      ylim=c(0, upy), main=NA, ylab=NA, xlab=NA, axes=FALSE, breaks=seq(0,1000, by=50))
# abline(h=1/up, col=cl[3])
# text(x=0, y=0.9*upy, labels=bquote(bar(x) == .(round(out5$mean$N[2,1]))), pos=4)
# text(x=0, y=0.8*upy, labels=paste('SD: ', round(out5$sd$N[2,1])), pos=4)
# axis(1)
# mtext('2y old', side=3, line=0.5)
# 
# par(mar=c(3,2,2,1))
# hist(out5$sims.list$N[,3,1], col=cl[1], border=cl[1], freq=FALSE, xlim=c(0, up),
#      ylim=c(0, upy), main=NA, ylab=NA, xlab=NA, axes=FALSE, breaks=seq(0,1000, by=50))
# abline(h=1/up, col=cl[3])
# text(x=0, y=0.9*upy, labels=bquote(bar(x) == .(round(out5$mean$N[3,1]))), pos=4)
# text(x=0, y=0.8*upy, labels=paste('SD: ', round(out5$sd$N[3,1])), pos=4)
# axis(1)
# mtext('3y old', side=3, line=0.5)
# 
# par(mar=c(3,2,2,1))
# hist(out5$sims.list$N[,4,1], col=cl[1], border=cl[1], freq=FALSE, xlim=c(0, up),
#      ylim=c(0, upy), main=NA, ylab=NA, xlab=NA, axes=FALSE, breaks=seq(0,1000, by=50))
# abline(h=1/up, col=cl[3])
# text(x=600, y=0.9*upy, labels=bquote(bar(x) == .(round(out5$mean$N[4,1]))), pos=4)
# text(x=600, y=0.8*upy, labels=paste('SD: ', round(out5$sd$N[4,1])), pos=4)
# axis(1)
# mtext('4y old', side=3, line=0.5)
# 
# par(mar=c(3,3,1,1))
# up <- 3000
# upy <- 0.0015
# hist(out6$sims.list$N[,1,1], col=cl[1], border=cl[1], freq=FALSE, xlim=c(0, up),
#      ylim=c(0, upy), main=NA, ylab=NA, xlab=NA, axes=FALSE, breaks=seq(0,3000, by=150))
# abline(h=1/up, col=cl[3])
# text(x=0, y=0.9*upy, labels=bquote(bar(x) == .(round(out6$mean$N[1,1]))), pos=4)
# text(x=0, y=0.8*upy, labels=paste('SD: ', round(out6$sd$N[1,1])), pos=4)
# axis(1)
# mtext('U(1,3000)', side=2, line=1)
# 
# par(mar=c(3,2,1,1))
# hist(out6$sims.list$N[,2,1], col=cl[1], border=cl[1], freq=FALSE, xlim=c(0, up),
#      ylim=c(0, upy), main=NA, ylab=NA, xlab=NA, axes=FALSE, breaks=seq(0,3000, by=150))
# abline(h=1/up, col=cl[3])
# text(x=0, y=0.9*upy, labels=bquote(bar(x) == .(round(out6$mean$N[2,1]))), pos=4)
# text(x=0, y=0.8*upy, labels=paste('SD: ', round(out6$sd$N[2,1])), pos=4)
# axis(1)
# 
# par(mar=c(3,2,1,1))
# hist(out6$sims.list$N[,3,1], col=cl[1], border=cl[1], freq=FALSE, xlim=c(0, up),
#      ylim=c(0, upy), main=NA, ylab=NA, xlab=NA, axes=FALSE, breaks=seq(0,3000, by=150))
# abline(h=1/up, col=cl[3])
# text(x=0, y=0.9*upy, labels=bquote(bar(x) == .(round(out6$mean$N[3,1]))), pos=4)
# text(x=0, y=0.8*upy, labels=paste('SD: ', round(out6$sd$N[3,1])), pos=4)
# axis(1)
# 
# par(mar=c(3,2,1,1))
# hist(out6$sims.list$N[,4,1], col=cl[1], border=cl[1], freq=FALSE, xlim=c(0, up),
#      ylim=c(0, upy), main=NA, ylab=NA, xlab=NA, axes=FALSE, breaks=seq(0,3000, by=150))
# abline(h=1/up, col=cl[3])
# text(x=1700, y=0.9*upy, labels=bquote(bar(x) == .(round(out6$mean$N[4,1]))), pos=4)
# text(x=1700, y=0.8*upy, labels=paste('SD: ', round(out6$sd$N[4,1])), pos=4)
# axis(1)
# 
# par(mar=c(3,3,1,1))
# up <- 1000
# upy <- 0.008
# hist(out7$sims.list$N[,1,1], col=cl[2], border=cl[2], freq=FALSE, xlim=c(0, up),
#      ylim=c(0, upy), main=NA, ylab=NA, xlab=NA, axes=FALSE, breaks=seq(0,1000, by=50))
# abline(h=1/up, col=cl[3])
# text(x=0, y=0.9*upy, labels=bquote(bar(x) == .(round(out7$mean$N[1,1]))), pos=4)
# text(x=0, y=0.8*upy, labels=paste('SD: ', round(out7$sd$N[1,1])), pos=4)
# axis(1)
# mtext('U(1,1000)', side=2, line=1)
# 
# par(mar=c(3,2,1,1))
# hist(out7$sims.list$N[,2,1], col=cl[2], border=cl[2], freq=FALSE, xlim=c(0, up),
#      ylim=c(0, upy), main=NA, ylab=NA, xlab=NA, axes=FALSE, breaks=seq(0,1000, by=50))
# abline(h=1/up, col=cl[3])
# text(x=0, y=0.9*upy, labels=bquote(bar(x) == .(round(out7$mean$N[2,1]))), pos=4)
# text(x=0, y=0.8*upy, labels=paste('SD: ', round(out7$sd$N[2,1])), pos=4)
# axis(1)
# 
# par(mar=c(3,2,1,1))
# hist(out7$sims.list$N[,3,1], col=cl[2], border=cl[2], freq=FALSE, xlim=c(0, up),
#      ylim=c(0, upy), main=NA, ylab=NA, xlab=NA, axes=FALSE, breaks=seq(0,1000, by=50))
# abline(h=1/up, col=cl[3])
# text(x=500, y=0.9*upy, labels=bquote(bar(x) == .(round(out7$mean$N[3,1]))), pos=4)
# text(x=500, y=0.8*upy, labels=paste('SD: ', round(out7$sd$N[3,1])), pos=4)
# axis(1)
# 
# par(mar=c(3,2,1,1))
# hist(out7$sims.list$N[,4,1], col=cl[2], border=cl[2], freq=FALSE, xlim=c(0, up),
#      ylim=c(0, upy), main=NA, ylab=NA, xlab=NA, axes=FALSE, breaks=seq(0,1000, by=50))
# abline(h=1/up, col=cl[3])
# text(x=500, y=0.9*upy, labels=bquote(bar(x) == .(round(out7$mean$N[4,1]))), pos=4)
# text(x=500, y=0.8*upy, labels=paste('SD: ', round(out7$sd$N[4,1])), pos=4)
# axis(1)
# 
# par(mar=c(3,3,1,1))
# up <- 3000
# upy <- 0.006
# hist(out8$sims.list$N[,1,1], col=cl[2], border=cl[2], freq=FALSE, xlim=c(0, up),
#      ylim=c(0, upy), main=NA, ylab=NA, xlab=NA, axes=FALSE, breaks=seq(0,3000, by=150))
# abline(h=1/up, col=cl[3])
# text(x=1000, y=0.9*upy, labels=bquote(bar(x) == .(round(out8$mean$N[1,1]))), pos=4)
# text(x=1000, y=0.8*upy, labels=paste('SD: ', round(out8$sd$N[1,1])), pos=4)
# axis(1)
# mtext('U(1,3000)', side=2, line=1)
# 
# par(mar=c(3,2,1,1))
# hist(out8$sims.list$N[,2,1], col=cl[2], border=cl[2], freq=FALSE, xlim=c(0, up),
#      ylim=c(0, upy), main=NA, ylab=NA, xlab=NA, axes=FALSE, breaks=seq(0,3000, by=150))
# abline(h=1/up, col=cl[3])
# text(x=1000, y=0.9*upy, labels=bquote(bar(x) == .(round(out8$mean$N[2,1]))), pos=4)
# text(x=1000, y=0.8*upy, labels=paste('SD: ', round(out8$sd$N[2,1])), pos=4)
# axis(1)
# 
# par(mar=c(3,2,1,1))
# hist(out8$sims.list$N[,3,1], col=cl[2], border=cl[2], freq=FALSE, xlim=c(0, up),
#      ylim=c(0, upy), main=NA, ylab=NA, xlab=NA, axes=FALSE, breaks=seq(0,3000, by=150))
# abline(h=1/up, col=cl[3])
# text(x=1500, y=0.9*upy, labels=bquote(bar(x) == .(round(out8$mean$N[3,1]))), pos=4)
# text(x=1500, y=0.8*upy, labels=paste('SD: ', round(out8$sd$N[3,1])), pos=4)
# axis(1)
# 
# par(mar=c(3,2,1,1))
# hist(out8$sims.list$N[,4,1], col=cl[2], border=cl[2], freq=FALSE, xlim=c(0, up),
#      ylim=c(0, upy), main=NA, ylab=NA, xlab=NA, axes=FALSE, breaks=seq(0,3000, by=150))
# abline(h=1/up, col=cl[3])
# text(x=1500, y=0.9*upy, labels=bquote(bar(x) == .(round(out8$mean$N[4,1]))), pos=4)
# text(x=1500, y=0.8*upy, labels=paste('SD: ', round(out8$sd$N[4,1])), pos=4)
# axis(1)
# par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 17.8 Specification of the survival process with hazard rates
# ============================================================

# ~~~~ Recreate the priors used for the first IPM (lines 36-56 above ~~~~
s <- 0.8                    # Guess of annual survival
f <- 0.25                   # Guess of recruitment

# Create matrix population model
A <- matrix(0, ncol=17, nrow=17)
A[1,2:17] <- f
for (a in 2:17){
  A[a,a-1] <- s
}

# Compute right eigenvector
z <- which.max(Re(eigen(A)$values))
revec <- Re(eigen(A)$vectors[,z])

# Population size in first age class of first year
r <- 0.5                    # Guess of reporting rate
h <- 0.1                    # Guess of hunting mortality
N1 <- elk$C[1,1] / (h * r)

# Compute age-dependent population size in first year
n <- N1 * revec / revec[1]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Bundle data and constants and produce data overview
nim.data <- with(elk, list(C=C, a=H[1,], b=H[2,], R=R, n=n))
nim.constants <- with(elk, list(total=colSums(R), n.years=ncol(C),
                                n.age=nrow(C)))


# Write NIMBLE model
elk.IPMcode4 <- nimbleCode({
  
  # Priors and linear models
  for (t in 1:(n.years-1)){
    f[t] ~ dunif(0, 1)
  }
  for (t in 1:n.years){
    # Priors for hazard rates
    mh[t] ~ dgamma(0.1, 0.1)                              # Hunting mortality hazard rate
    mo[t] ~ dgamma(0.1, 0.1)                              # Background mortality hazard rate
    # Calculate probabilities from hazard rates
    s[t] <- exp(-(mh[t] + mo[t]))                         # Overall survival
    h[t] <- (1 - s[t]) * (mh[t] / (mh[t] + mo[t]))        # Hunting mortality
    o[t] <- (1 - s[t]) * (mo[t] / (mh[t] + mo[t]))        # Background mortality
    # Prior for reporting rate
    r[t] ~ dunif(0, 1)
  }
  
  # Age-at-harvest data (state-space model)
  # Model for the initial population size: Poisson priors
  for (a in 1:n.age){
    N[a,1] ~ dpois(n[a])
  }
  
  # Process model over time: our model of population dynamics
  for (t in 1:(n.years-1)){
    N[1,t+1] ~ dpois((Ntot[t] - N[1,t]) * f[t])
    for (a in 2:n.age){
      N[a,t+1] ~ dbin((1-h[t]) * s[t], N[a-1,t])
    } #a
  } #t
  
  # Derived quantity: total year-specific population size
  for (t in 1:n.years){
    Ntot[t] <- sum(N[1:n.age,t])
  }
  
  # Observation model
  for (t in 1:n.years){
    for (a in 1:n.age){
      C[a,t] ~ dbin(h[t] * r[t], N[a,t])
    } #a
  } #t
  
  # Hunter survey data (logistic regression model)
  for (t in 1:n.years){
    b[t] ~ dbin(r[t], a[t])
  }
  
  # Radio tracking data (multinomial model)
  for (t in 1:n.years){
    R[1:3,t] ~ dmulti(prt[1:3,t], total[t])
    prt[1,t] <- h[t]
    prt[2,t] <- (1-h[t]) * (1-s[t])
    prt[3,t] <- (1-h[t]) * s[t]
  }
})


# Initial value simulation function
inits <- function() {list(mh=runif(nim.constants$n.years, 0.01, 0.1), mo=runif(nim.constants$n.years,
                                                                           0.01, 0.1))}
# Parameters monitored
parameters <- c("mh", "mo", "h", "s", "o", "f", "r", "Ntot", "N")

# MCMC settings
#ni <- 10; nb <- 0; nc <- 3; nt <- 1
ni <- 450000; nb <- 50000; nc <- 3; nt <- 40

# Pre-sample initial values
Inits <- list()
for(i in 1:nc){
  Inits[[i]] <- inits()
}

# Run NIMBLE (ART x min) and check convergence
out9 <- nimbleMCMC(code = elk.IPMcode4,
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

plot(out9, ask = TRUE)
print(summary(out9))
