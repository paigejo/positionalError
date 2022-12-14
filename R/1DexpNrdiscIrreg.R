# simple one dimensional example, using the exponential covariance, nugget, and a 
# discrete (0 or 1) covariate in space. Include uneven grid spacing

require(spam)
require(fields)

# set working directory to folder containing this script
code.path <- rstudioapi::getActiveDocumentContext()$path
code.path.splitted <- strsplit(code.path, "/")[[1]]
home.dir <- paste(code.path.splitted[1: (length(code.path.splitted)-2)], collapse = "/")
setwd(home.dir)

theta0 <- c(.1, 1,.1)  # range, sill, nugget
beta0 <- c(0,1)

n <- 150
R <- 25
DeltaFact <- 0.5
perDelta <- DeltaFact * 1/(2*(n+1))  # should be between 0 and 1/(2n+2)

nKnots=10
alpha = 1
covScale = 1/(2*(n+1)) * 1
jitterDist = 1/(2*(n+1))
relUrbProb = 1 # relative probability of selecting points when covariate=1
filename <- paste0("R",R,"_n",n,"_th",theta0[1]*10,theta0[2],
                   theta0[3]*10,"_D",
                   DeltaFact*10,"_b",beta0[2], 
                   "_discK", nKnots, "_a", alpha, "_sc", round(10^3*covScale, 1), 
                   "_jitt", round(jitterDist*10^3, 1), "_rPr", relUrbProb)
set.seed(12)

# covScale is inverse range, alpha is amplitude
fun = function(x, covSc=covScale, a=alpha, 
               nK=nKnots, knots=seq(0, 1, l=nK)) {
  dists = rdist(cbind(x, 0), cbind(knots, 0))
  minDists = apply(dists, 1, min)
  
  a * as.numeric(minDists < covScale)
}
m = 100001
stemp = seq(0, 1, l=m)
stemp = stemp[-m]
delta = 1/(m-1)
plot(stemp, fun(stemp), type="l")

# calculate the probability of sampling each of stemp such that 
# probability of sampling point with fun(stemp) == 1 is relUrbProb 
# larger than with fun(stemp) == 0:

# (nKnots-1)*covScale
# baseProb = 1/((nKnots + 1) * (1 - (nKnots- 1)*covScale))
# urbProb = nKnots / (covScale * (nKnots - 1) * (nKnots + 1))
baseProb = as.numeric(fun(stemp) == 0) / sum(fun(stemp) == 0)
urbProb = relUrbProb * as.numeric(fun(stemp) == 1) / sum(fun(stemp) == 1)
probs = baseProb + urbProb
probs = probs/sum(probs)

probs = 1 + fun(stemp)
probs = probs/sum(probs)

# function for sampling locations randomly based on the probabilties above
sampleLocs = function(thisJitterDist=jitterDist) {
  thisSper = sample(stemp, n, replace=TRUE, prob=probs)
  thisSper = thisSper + runif(n, max=delta)
  
  thisSapp = thisSper + runif(n, min=-thisJitterDist, max=thisJitterDist)
  
  list(sper=thisSper, sapp=thisSapp)
}

# plot an example set of sampling locations
out = sampleLocs()
sapp = out$sapp
sper = out$sper
plot(stemp, fun(stemp), type="l", xlim=c(-jitterDist, 1+jitterDist), ylim=c(0, 1))
arrows(sper, rep(0, n), sper, rep(0.2, n), col="blue", length=0)
arrows(sapp, rep(0, n), sapp, rep(0.2, n), col="red", length=0)

# approximate and true (per for perturbed is the truth) spatial coordinates

# plot(sper, fun(sper), type="l")
# plot(sapp, fun(sper), type="l")
fun(sapp) - fun(sper)
hist(fun(sapp) - fun(sper), breaks=30)
mean(fun(sapp))
mean(fun(sper))
mean(fun(sapp) != fun(sper))
mean(abs(mean(fun(sapp) - fun(sper))))
varxstar = var(fun(sper))
varx = var(fun(sapp))
print(varxstar/varx) # beta1 should converge to this factor of the truth
print(varxstar/varx)

misclassRate = mean(fun(sper)[fun(sapp) == 1]) - mean(fun(sper)[fun(sapp) == 0])
print(misclassRate)

happ <- as.matrix( dist( sapp))
hper <- as.matrix( dist( sper))

happ <- as.matrix( dist( sapp))
Xapp <- cbind(1, fun(sapp))

# yR <- rmvnorm(R, mu=Xper%*%beta0, Sigma=cov.exp(hper, theta0))

Rsel <- 3

parResult <- matrix(0, R, 15)
meanXper = numeric(R)
meanXapp = numeric(R)
meanAbsDiff = numeric(R)
varxstars = numeric(R)
varxs = numeric(R)
misclassRates = numeric(R)

startTime = proc.time()[3]
for (Rsel in 1:R) {
  # simulate perturbed locations, covariates, and responses
  out = sampleLocs()
  sapp = out$sapp
  sper = out$sper
  
  Xapp <- cbind(1, fun(sapp))
  happ <- as.matrix( dist( sapp))
  Xper <- cbind(1, fun(sper))
  hper <- as.matrix( dist( sper))
  yR <- rmvnorm(1, Sigma=cov.exp(hper, theta0))
  
  y <- c(yR) + Xper%*%beta0
  
  meanXper[Rsel] = mean(fun(sper))
  meanXapp[Rsel] = mean(fun(sapp))
  meanAbsDiff[Rsel] = mean(abs(fun(sper) - fun(sapp)))
  varxstars[Rsel] = var(fun(sper))
  varxs[Rsel] = var(fun(sapp))
  misclassRates[Rsel] = mean(fun(sper)[fun(sapp) == 1]) - mean(fun(sper)[fun(sapp) == 0])
  
  # (betaMLE, thetaMLE, 2*nll, convergence)
  (parsper <- mle( y, Xper, hper, cov.exp, beta0, theta0, thetalower=c(0.01,0.01,0.01),
                   thetaupper=Inf)[c(1,2,4)]) 
  (parsapp <- mle( y, Xapp, happ, cov.exp, beta0, theta0, thetalower=c(0.01,.01,.01),
                   thetaupper=Inf)[c(1,2,4)]) 
  
  
  parResult[Rsel,] <- c(unlist(parsper),unlist(parsapp),
                        neg2loglikelihood(y, Xper, hper, cov.exp, parsapp$par[1:2], parsapp$par[3:5]))
  
  currentTime = proc.time()[3]
  propDone = Rsel/R
  propLeft = 1-propDone
  timeTaken = (currentTime - startTime) / 60
  totalTimeEst = timeTaken / propDone
  timeLeftEst = totalTimeEst * propLeft
  print(paste0("realization ", Rsel, "/", R, " complete"))
  print(paste0("Time taken (minutes): ", timeTaken))
  print(paste0("Estimated time left (minutes): ", timeLeftEst))
}

parResult = cbind(parResult, meanAbsXDiff=meanAbsDiff, meanXDiff=meanXper-meanXapp, meanXper=meanXper, meanXapp=meanXapp, 
                  varxstar=varxstars, varx=varxs, estRelBiasBeta1=varxstars/varxs, misclassRate=misclassRates)
colnames(parResult) = c("beta0per", "beta1per", "rangeper", "sillper", "nuggetper", "2nllper", "convergenceper", 
                        "beta0app", "beta1app", "rangeapp", "sillapp", "nuggetapp", "2nllapp", "convergenceapp", 
                        "2nllperAtApp", "meanAbsXDiff", "meanXDiff", "meanXper", "meanXapp", 
                        "varxstar", "varx", "estRelBiasBeta1", "misclassRate")
print(parResult)

save.image(file = paste0(filename,".RData"))
if(FALSE) {
  code.path <- rstudioapi::getActiveDocumentContext()$path
  code.path.splitted <- strsplit(code.path, "/")[[1]]
  home.dir <- paste(code.path.splitted[1: (length(code.path.splitted)-2)], collapse = "/")
  setwd(home.dir)
  
  theta0 <- c(.1, 1,.1)  # range, sill, nugget
  beta0 <- c(0,1)
  
  n <- 150
  R <- 25
  DeltaFact <- 0.5
  perDelta <- DeltaFact * 1/(2*(n+1))  # should be between 0 and 1/(2n+2)
  
  nKnots=10
  alpha = 1
  covScale = 1/(2*(n+1)) * 1
  jitterDist = 1/(2*(n+1))
  relUrbProb = 1 # relative probability of selecting points when covariate=1
  filename <- paste0("R",R,"_n",n,"_th",theta0[1]*10,theta0[2],
                     theta0[3]*10,"_D",
                     DeltaFact*10,"_b",beta0[2], 
                     "_discK", nKnots, "_a", alpha, "_sc", round(10^3*covScale, 1), 
                     "_jitt", round(jitterDist*10^3, 1), "_rPr", relUrbProb)
  
  load(file = paste0(filename,".RData"))
}

diff3 <- parResult[,1:6]- parResult[,8:13] # good pars and 2nll - bad pars and 2nll
truePar = c(beta0, theta0)
biasper = sweep(parResult[,1:5], 2, truePar)
biasapp = sweep(parResult[,8:12], 2, truePar)

# pair plot of good model estimates
pairs( parResult[,c(1:6, 16, 22)], gap=0, pch=19, cex=.4)
# pair plot of difference between good and approximate model estimates
pairs( cbind(diff3, meanXDiff=parResult[,22]),gap=0, pch=19,cex=.4)

# calculate bias relative to the truth
relBiasper = as.data.frame(sweep(parResult[,2:5], 2, truePar[-1], FUN=function(x, y) {(x - y)/y})*100)
relBiasapp = as.data.frame(sweep(parResult[,9:12], 2, truePar[-1], FUN=function(x, y) {(x - y)/y})*100)
relBiasper = cbind(relBiasper, meanAbsXDiff=parResult[,16], estPctBiasBeta1=100*(parResult[,22]-1), misclassPct=100*(parResult[,23]-1))
relBiasapp = cbind(relBiasapp, meanAbsXDiff=parResult[,16], estPctBiasBeta1=100*(parResult[,22]-1), misclassPct=100*(parResult[,23]-1))
colMeans(parResult[,2:5])
relBiasRange = range(c(relBiasper, relBiasapp))
boxplot(as.list(relBiasper), ylim=relBiasRange)
boxplot(as.list(relBiasapp), ylim=relBiasRange)
colMeans(relBiasper)
colMeans(relBiasapp)
pairs(relBiasapp, gap=0, pch=19,cex=.4)

pairs(relBias, gap=0, pch=19,cex=.4)

pairs(biasper, gap=0, pch=19, cex=.4)
pairs(biasapp, gap=0, pch=19, cex=.4)

summary(diff3)


hist(parResult[,15]-parResult[,6],breaks = 25)


