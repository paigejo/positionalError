# simple onedimesional example, using the exponential covariance.

require(spam)
require(fields)

theta0 <- c(.3, 1,.1)  # range, sill 
beta0 <- c(0,1)

n <- 150
R <- 100
DeltaFact <- 0.5
perDelta <- DeltaFact * 1/(2*(n+1))  # should be between 0 and 1/(2n+2)


filename <- paste0("R",R,"_n",n,"_th",theta0[1]*10,theta0[2],
                   theta0[3]*10,"_D",
                   DeltaFact*10,"_b",beta0[2])
set.seed(12)

xapp <- c(1:n)/(1+n)
xper <- xapp + runif(n, -perDelta, perDelta)

X <- cbind(1, xper)

happ <- as.matrix( dist( xapp))
hper <- as.matrix( dist( xper))

yR <- rmvnorm(R, mu=X%*%beta0, Sigma=cov.exp(hper, theta0))

Rsel <- 3

parResult <- matrix(0, R, 15)


for (Rsel in 1:R) {

y <- c(yR[Rsel,])

(parsper <- mle( y, X, hper, cov.exp, beta0, theta0, thetalower=c(0.01,0.01,0.01),
             thetaupper=Inf)[c(1,2,4)]) 
(parsapp <- mle( y, X, happ, cov.exp, beta0, theta0, thetalower=c(0.01,.01,.01),
             thetaupper=Inf)[c(1,2,4)]) 


parResult[Rsel,] <- c(unlist(parsper),unlist(parsapp),
                      neg2loglikelihood(y, X, hper, cov.exp, parsapp$par[1:2], parsapp$par[3:5]))

}


print(parResult)

save.image(file = paste0(filename,".RData"))

diff3 <- parResult[,1:6]- parResult[,9:14]


pairs( parResult[,1:6], gap=0, pch=19, cex=.4)
pairs( diff3,gap=0, pch=19,cex=.4)

summary(diff3)


hist(parResult[,15]-parResult[,6],breaks = 25)


