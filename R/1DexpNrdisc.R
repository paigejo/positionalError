# simple one dimensional example, using the exponential covariance, nugget, and a 
# discrete (0 or 1) covariate in space.

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
filename <- paste0("R",R,"_n",n,"_th",theta0[1]*10,theta0[2],
                   theta0[3]*10,"_D",
                   DeltaFact*10,"_b",beta0[2], 
                   "_discnK", nKnots, "_a", alpha, "_sc", covScale)
set.seed(12)

# approximate and true (per for perturbed is the truth) spatial coordinates
sapp <- c(1:n)/(1+n)
sper <- sapp + runif(n, -perDelta, perDelta) # just a test set of perturbed locations

# covScale is inverse range, alpha is amplitude
fun = function(x, covSc=covScale, a=alpha, 
               nK=nKnots, knots=seq(0, 1, l=nK)) {
  dists = rdist(cbind(x, 0), cbind(knots, 0))
  minDists = apply(dists, 1, min)
  
  a * as.numeric(minDists < covScale)
}
stemp = seq(0, 1, l=1001)
plot(stemp, fun(stemp), type="l")
plot(sper, fun(sper), type="l")
plot(sapp, fun(sper), type="l")
fun(sapp) - fun(sper)
hist(fun(sapp) - fun(sper), breaks=30)
mean(fun(sapp))
mean(fun(sper))

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
meanXapp = mean(fun(sapp))
meanAbsDiff = numeric(R)
varxstars = numeric(R)
varxs = var(fun(sapp))
misclassRates = numeric(R)

startTime = proc.time()[3]
for (Rsel in 1:R) {
  # simulate perturbed locations, covariates, and responses
  sper <- sapp + runif(n, -perDelta, perDelta)
  Xper <- cbind(1, fun(sper))
  hper <- as.matrix( dist( sper))
  yR <- rmvnorm(1, Sigma=cov.exp(hper, theta0))
  
  y <- c(yR) + Xper%*%beta0
  
  meanXper[Rsel] = mean(fun(sper))
  meanAbsDiff[Rsel] = mean(abs(fun(sper) - fun(sapp)))
  varxstars[Rsel] = var(fun(sper))
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

parResult = cbind(parResult, meanAbsXDiff=meanAbsDiff, meanXDiff=meanXper-meanXapp, 
                  varxstar=varxstars, varx=varx, estRelBiasBeta1=varxstars/varxs, misclassRate=misclassRates)
colnames(parResult) = c("beta0per", "beta1per", "rangeper", "sillper", "nuggetper", "2nllper", "convergenceper", 
                        "beta0app", "beta1app", "rangeapp", "sillapp", "nuggetapp", "2nllapp", "convergenceapp", 
                        "2nllperAtApp", "meanAbsXDiff", "meanXDiff", 
                        "varxstar", "varx", "estRelBiasBeta1", "misclassRate")
print(parResult)

save.image(file = paste0(filename,".RData"))

diff3 <- parResult[,1:6]- parResult[,8:13] # good pars and 2nll - bad pars and 2nll
truePar = c(beta0, theta0)
biasper = sweep(parResult[,1:5], 2, truePar)
biasapp = sweep(parResult[,8:12], 2, truePar)

# pair plot of good model estimates
pairs( parResult[,c(1:6, 16)], gap=0, pch=19, cex=.4)
# pair plot of difference between good and approximate model estimates
pairs( cbind(diff3, meanXDiff=parResult[,16]),gap=0, pch=19,cex=.4)

# calculate bias relative to the truth
relBiasper = as.data.frame(sweep(parResult[,2:5], 2, truePar[-1], FUN=function(x, y) {(x - y)/y})*100)
relBiasapp = as.data.frame(sweep(parResult[,9:12], 2, truePar[-1], FUN=function(x, y) {(x - y)/y})*100)
relBiasper = cbind(relBiasper, misclassPct=100*(parResult[,21]-1))
relBiasapp = cbind(relBiasapp, misclassPct=100*(parResult[,21]-1))
colMeans(parResult[,2:5])
relBiasRange = range(c(relBiasper, relBiasapp))
boxplot(as.list(relBiasper), ylim=relBiasRange)
boxplot(as.list(relBiasapp), ylim=relBiasRange)
colMeans(cbind(relBiasper, m2ll=parResult[,6]))
colMeans(cbind(relBiasapp, m2ll=parResult[,15]))
pairs(relBiasapp, gap=0, pch=19,cex=.4)

pairs(relBias, gap=0, pch=19,cex=.4)

pairs(biasper, gap=0, pch=19, cex=.4)
pairs(biasapp, gap=0, pch=19, cex=.4)

summary(diff3)


hist(parResult[,15]-parResult[,6],breaks = 25)


