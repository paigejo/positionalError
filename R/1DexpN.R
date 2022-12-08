# simple onedimesional example, using the exponential covariance.

require(spam)
require(fields)

theta0 <- c(.3, 1,.1)  # range, sill 

n <- 150
R <- 100
DeltaFact <- 0.5
perDelta <- DeltaFact * 1/(2*(n+1))  # should be between 0 and 1/(2n+2)


filename <- paste0("R",R,"_n",n,"_th",theta0[1]*10,theta0[2],
                   theta0[3]*10,"_D",
                   DeltaFact*10)
set.seed(12)

xapp <- c(1:n)/(1+n)
xper <- xapp + runif(n, -perDelta, perDelta)

happ <- as.matrix( dist( xapp))
hper <- as.matrix( dist( xper))

yR <- rmvnorm(R, Sigma=cov.exp(hper, theta0))
yR <- yR - rowMeans(yR)

Rsel <- 3

parResult <- matrix(0, R, 11)


for (Rsel in 1:R) {

y <- c(yR[Rsel,])

(parsper <- mle.nomean( y, hper, cov.exp, theta0, thetalower=c(.01, .01,.01),
             thetaupper=Inf)[c(1,2,4)]) 
(parsapp <- mle.nomean( y, happ, cov.exp, theta0, thetalower=c(.01,.01, .01),
             thetaupper=Inf)[c(1,2,4)]) 


parResult[Rsel,] <- c(unlist(parsper),unlist(parsapp),neg2loglikelihood.nomean(y, hper, cov.exp, parsapp$par))

}


print(parResult)

save.image(file = paste0(filename,".RData"))

diff3 <- parResult[,1:4]- parResult[,6:9]

plot( diff3[,1:2], xlab='Range', ylab='sill')
abline(v=0,h=0, col='gray') 
summary(diff3)


hist(parResult[,11]-parResult[,4],breaks = 25)
abline(v=0)

