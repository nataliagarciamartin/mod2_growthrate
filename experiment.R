library(rlist)
library(shape)
source("aux_functions.R")
source("run_markov.R")

################################################################################
# Parameters

N <- 3 # number og age classes
M <- 5 # number of Leslie matrices

b_1 <- 1.5 # birth rate of class 1
b_2 <- 0.8 # birth rate of class 2
lambda_1 <- 0.9 # survival rate of class 0 (juvenile class)
lambda_2 <- 0.7 # survival rate of class 1

eps <- 0.2 # semi-range of parameters

T <- 500 # number of time steps
pop <- 100 # initial number of individuals  
K <- 10^6 # maximal number of individuals  

################################################################################
# Simulation: 5 chains

ens_M <- create_M(b_1, b_2, lambda_1, lambda_2, eps, M)
p <- rep(1,M) / M # probability distribution on M

w_01 <- c(0.2,0.4,0.4) # initial state 1
results_1 <- run_markov(w_01, pop, ens_M, p, T)

w_02 <- c(0.6,0.2,0.2) # initial state 2
results_2 <- run_markov(w_02, pop, ens_M, p, T)

w_03 <- c(0.4,0.3,0.3) # initial state 3
results_3 <- run_markov(w_03, pop, ens_M, p, T)

w_04 <- c(0.3,0.3,0.4) # initial state 4
results_4 <- run_markov(w_04, pop, ens_M, p, T)

w_05 <- c(0.3,0.3,0.4) # initial state 5
results_5 <- run_markov(w_05, pop, ens_M, p, T)


################################################################################
# Plots 5 curves

a_1 <- cumsum(results_1$logp - results_1$logp[1]) / c(1, seq_len(T))
a_2 <- cumsum(results_2$logp - results_2$logp[1]) / c(1, seq_len(T))
a_3 <- cumsum(results_3$logp - results_3$logp[1]) / c(1, seq_len(T))
a_4 <- cumsum(results_4$logp - results_4$logp[1]) / c(1, seq_len(T))
a_5 <- cumsum(results_5$logp - results_5$logp[1]) / c(1, seq_len(T))


mean_a <- mean(c(a_1[100:(T+1)], a_2[100:(T+1)], a_3[100:(T+1)], a_4[100:(T+1)],
                 a_5[100:(T+1)]))
plot(seq_len(T), a_1[-1], type='l', main='Value of the average growth rate along the
     times steps t', xlab='time step t', ylab='average growth rate a^(t)')
abline(h=mean_a, col="red", lty=2, lwd=2)
lines(seq_len(T), a_2[-1], type='l',col="green")
lines(seq_len(T), a_3[-1], type='l',col="blue")
lines(seq_len(T), a_4[-1], type='l',col="pink")
lines(seq_len(T), a_5[-1], type='l',col="purple")

################################################################################
# Simulation: J chains

J <- 1000

ratio <- numeric(J) # stores the values of the estimated ratio
log_pop <- matrix(0, nrow = T+1, ncol = J) # stores the values of the population
w_T <- matrix(0, nrow = N, ncol = J) # stores the values of w_T

for (j in seq_len(J)){
  
  rdm <- runif(N-1)
  w_0 <- c(rdm, 1-sum(rdm))
  
  results <- run_markov(w_0, pop, ens_M, p, T)
  
  log_pop[,j] <- results$logp
  
  a <- cumsum(results$logp - results$logp[1]) / c(1, seq_len(T))
  ratio[j] <- a[T+1]
  w_T[,j] <- results$w[T+1,]
  
}

################################################################################
# Plots: convergence of mean and variance of aj (CLT)  

a_j <- cumsum(ratio) / seq_len(J)
var_j <- cumsum((ratio - cumsum(ratio)/seq_len(J))^2) / c(1,seq_len(J-1))

plot(seq_len(J), a_j, type='l', main='Mean value of the average growth rate vs
     number of chains', xlab='J', ylab='average growth rate a_J^(t)')
abline(h=mean(ratio), col="red", lty=2, lwd=2)
plot(seq_len(J), var_j, type='l', main='Estimated variance of the average growth
     rate vs number of chains', xlab='J', ylab='variance of a_J^(T)')
abline(h=var(ratio), col="red", lty=2, lwd=2)

################################################################################
# Plots: interval of confidence (95 %)

mean_aj <- mean(ratio)
sup_ic <- mean_aj + 2.365 * sd(ratio)
inf_ic <- mean_aj - 2.365 * sd(ratio)

# plot(seq_len(J), ratio)
# abline(h=mean_aj, col="red", lty=2, lwd=2)
# abline(h=sup_ic, col="red", lty=2, lwd=2)
# abline(h=inf_ic, col="red", lty=2, lwd=2)


a_1 <- cumsum(log_pop[,1] - log_pop[1,1]) / c(1, seq_len(T))
a_2 <- cumsum(log_pop[,2] - log_pop[1,1]) / c(1, seq_len(T))
a_3 <- cumsum(log_pop[,3] - log_pop[1,1]) / c(1, seq_len(T))
a_4 <- cumsum(log_pop[,4] - log_pop[1,1]) / c(1, seq_len(T))
a_5 <- cumsum(log_pop[,5] - log_pop[1,1]) / c(1, seq_len(T))

plot(seq_len(T), a_1[-1], type='l', main='Value of the average growth rate along the
     times steps t', xlab='time step t', ylab='average growth rate a^(t)', 
     ylim=c(-4.4, -4.2))
lines(seq_len(T), a_2[-1], type='l',col="green")
lines(seq_len(T), a_3[-1], type='l',col="blue")
lines(seq_len(T), a_4[-1], type='l',col="pink")
lines(seq_len(T), a_5[-1], type='l',col="purple")
abline(h=mean_aj, col="red", lty=2, lwd=2)
abline(h=sup_ic, col="red", lty=2, lwd=2)
abline(h=inf_ic, col="red", lty=2, lwd=2)

# histogram of the average growth rate
hist(ratio, breaks=15,freq=F, main="Histogram and experimental fitted density of
     the average growth rate", xlab="a_t")
lines(density(ratio), col="red", lwd=2) 
curve(dnorm(x, mean=mean(ratio), sd=sd(ratio)), add=TRUE, lty="dotted",
      col="darkgreen", lwd=3)


################################################################################
E_L <- create_M(b_1, b_2, lambda_1, lambda_2, eps, 1)[[1]]
eig <- eigen(E_L)
if (any(sapply(eig$values, is.complex))){
  a_true <- log(Re(eig$values[1]))
  top_eig <- Re(eig$vectors[,1])
} else {
  a_true <- log(max(eig$values))
  k <- which.max(eig$values)
  top_eig <- eig$vectors[,k]
}

eig_vectors <- top_eigen(ens_M)$vectors

################################################################################
# Plots

#plot(seq_len(T+1), log_pop[,1], type='l')
#print(mean(ratio))
#print(a_true)

plot(seq_len(J), moy, type='l')
#abline(h=a_true, col="blue")
#plot(seq_len(J), var, type='l')

plot(w_T[2,], w_T[3,])
#Arrows(0, 0, top_eig[2], top_eig[3], col="blue")
Arrows(0.37, 0.135, 0, eig_vectors[3,1]/eig_vectors[1,1], col="blue")
#Arrows(c(0,0,0), c(0,0,0), abs(eig_vectors[2,]), abs(eig_vectors[3,]), col="blue")


#abline(a=0, b=eig_vectors[3,1]/eig_vectors[2,1], col="blue")

plot(seq_len(5),seq_len(5))
abline(a=5, b=-1, col="blue")
