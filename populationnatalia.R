# population plot

N = 3
T = 1000
b = 1.5
c = 0.9
d = 0.8
e = 0.5
eps = 0.3
pop = 100
J = 1000
v_0 <- c(0.5,0.3,0.2)

sample_L <- function(b, c, d, e, eps){
  L <- matrix(0,3,3)
  rdm <- runif(4, -eps, eps)
  L[1,2] = b + rdm[1]
  L[1,3] = e + rdm[4]
  L[2,1] = c + rdm[2]
  L[3,2] = d + rdm[3]
  return(L)
}

population <- numeric(T+1)
U <- v_0
distrib <- matrix(0,T+1,N)
population[1] <- pop
distrib[1] <- v_0

for (t in seq_len(T)){
  L <- sample_L(b,c,d,e,eps)
  U=L%*%U
  population[t+1]=sum(U)
  distrib[t+1,]=U/population[t+1]
}
plot(population, type="l", ylab="Population size", xlab="Time")


plot1=plot(distrib[,1], col=1, type="l", ylim=c(0,1))+lines(distrib[,2], col=2, type="l", lty=1)+lines(distrib[,3], col=3, type="l", lty=1)
legend(legend=c("Proportion Age 1", "Proportion Age 2", "Proportion Age 3"), col=1:3, "topright", lty=1)



