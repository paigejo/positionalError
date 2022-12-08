

load("R25_n300_th31_D5.RData")
parResult5 <- parResult

load("R25_n300_th31_D9.RData")
parResult9 <- parResult

load("R100_n150_th31_D5.RData")
parResult5s <- parResult

diff5s <- parResult5s[,1:3]- parResult5s[,5:7]

diff5 <- parResult5[,1:3]- parResult5[,5:7]
diff9 <- parResult9[,1:3]- parResult9[,5:7]

plot( diff5s[,1:2], xlab='Range', ylab='sill')
abline(v=0,h=0, col='gray') 
summary(diff5s)



plot( diff5[,1:2], xlab='Range', ylab='sill')
abline(v=0,h=0, col='gray') 
summary(diff5)

plot( diff9[,1:2], xlab='Range', ylab='sill')
abline(v=0,h=0, col='gray') 
summary(diff9)



load("R100_n150_th31_D5.RData")
p1 <- parResult


load("R100_n150_th311_D5_b1.RData")
p1 <- parResult




