# population plot

N = 3
T = 200
b = 1.5
c = 0.9
d = 0.8
e = 0.5
eps = 0.3
pop = 100
J = 1000
v_0 <- c(0.2,0.2,0.6)

sample_L <- function(b, c, d, e, eps){
  L <- matrix(0,3,3)
  rdm <- runif(4, -eps, eps)
  L[1,2] = b + rdm[1]
  L[1,3] = e + rdm[4]
  L[2,1] = c + rdm[2]
  L[3,2] = d + rdm[3]
  return(L)
}

population <- numeric(T)
U <- v_0
distrib <- matrix(0,T,N)

for (t in seq_len(T)){
  L <- sample_L(b,c,d,e,eps)
  U=L%*%U
  population[t]=sum(U)
  distrib[t,]=U/population[t]
}
plot(c(pop, population), type="l", ylab="Population size", xlab="Time")

distrib

plot1=plot(distrib[,1], col=1, type="l")+lines(distrib[,2], col=2, type="l", lty=2)+lines(distrib[,3], col=3, type="l", lty=3, ylim=c(0,1))
legend(legend=c("Age 1", "Age 2", "Age 3"), col=1:3, "topright", lty=1:3)




v_0 <- c(0.6,0.2,0.2)
ratio = numeric(J)

for (j in seq_len(J)){
  
  v <- v_0
  for (t in seq_len(T)){
    
    L <- sample_L(b,c,d,e,eps)
    
    if (t == T){
      
      ratio[j] = log(sum(L %*% v)) - log(sum(v))
    }
    
    else v <- L %*% v
  }
  
}
