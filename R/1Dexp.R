# simple onedimesional example, using the exponential covariance.

require(spam)
require(fields)

theta0 <- c(.3, 1)  # range, sill 

n <- 150
R <- 100
DeltaFact <- 0.5
perDelta <- DeltaFact * 1/(2*(n+1))  # should be between 0 and 1/(2n+2)


filename <- paste0("R",R,"_n",n,"_th",theta0[1]*10,theta0[2],"_D",
                   DeltaFact*10)
set.seed(12)

xapp <- c(1:n)/(1+n)
xper <- xapp + runif(n, -perDelta, perDelta)

happ <- as.matrix( dist( xapp))
hper <- as.matrix( dist( xper))

yR <- rmvnorm(R, Sigma=cov.exp(hper, theta0))
yR <- yR - rowMeans(yR)

Rsel <- 3

parResult <- matrix(0, R, 8)

pdf(paste0(filename,".pdf"))

for (Rsel in 1:R) {

y <- c(yR[Rsel,])

par(mfrow=c(2,2))
plot(xper, y, pch=19, cex=.2)

(parsper <- mle.nomean( y, hper, cov.exp, theta0, thetalower=c(.01,.01),
             thetaupper=Inf)[c(1,2,4)]) 
(parsapp <- mle.nomean( y, happ, cov.exp, theta0, thetalower=c(.01,.01),
             thetaupper=Inf)[c(1,2,4)]) 
parResult[Rsel,] <- c(unlist(parsper),unlist(parsapp))

nr <- 39
ns <- 41
range <- seq(.05, to=.5, length=nr)
sill <- seq(.3, 2.5, length=ns)
grid <- expand.grid( range, sill)


n2llper <- apply( grid, 1, function(theta)
  neg2loglikelihood.nomean( y, hper, cov.exp, theta))
reln2llper <- matrix( n2llper - parsper$value, nr, ns) 

n2llapp <- apply( grid, 1, function(theta)
  neg2loglikelihood.nomean( y, happ, cov.exp, theta))
reln2llapp <- matrix( n2llapp - parsapp$value, nr, ns) 

image.plot( range, sill, reln2llper, zlim=c(0,12), main="exact (perturbed)")
contour( range, sill, reln2llper, levels=qchisq(c(.7,.9), 2),
         labels=c('70%','90%'), labcex=1, add=TRUE)
points( parsper$par[1], parsper$par[2], pch=20, col='white')
points( parsapp$par[1], parsapp$par[2], pch=20, col='yellow')

image.plot( range, sill, reln2llapp, zlim=c(0,12), main="approx (gridded)")
contour( range, sill, reln2llapp, levels=qchisq(c(.7,.9), 2),
         labels=c('70%','90%'), labcex=1, add=TRUE)
points( parsper$par[1], parsper$par[2], pch=20, col='white')
points( parsapp$par[1], parsapp$par[2], pch=20, col='yellow')

image.plot( range, sill, reln2llapp-reln2llper, zlim=c(-10,10), 
            main="rel. neg-2-loglik: approx - exact")
points( parsper$par[1], parsper$par[2], pch=20, col='white')
points( parsapp$par[1], parsapp$par[2], pch=20, col='yellow')


cat('.')
}


dev.off()

print(parResult)

save.image(file = paste0(filename,".RData"))

}
