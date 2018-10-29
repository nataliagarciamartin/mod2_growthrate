library(rlist)
library(shape)
source("aux_functions.R")

################################################################################
# Parameters

N <- 3 # number og age classes
M <- 2 # number of Leslie matrices

b_1 <- 1.5 # birth rate of class 1
b_2 <- 0.8 # birth rate of class 2
lambda_1 <- 0.9 # survival rate of class 0 (juvenile class)
lambda_2 <- 0.7 # survival rate of class 1

eps <- 0.2 # semi-range of parameters

T <- 1000 # number of time steps
pop <- 100 # initial number of individuals  
K <- 10^6 # maximal number of individuals  
J <- 500 # number of Monte Carlo samples

################################################################################
# Simulations

w_0 <- c(0.2,0.4,0.4) # initial state

ens_M <- create_M(b_1, b_2, lambda_1, lambda_2, eps, M)
p <- c(0.2, 0.7, 0.1) # probability distribution on M

ratio <- numeric(J) # stores the values of the estimated ratio
log_pop <- matrix(0, nrow = T+1, ncol = J) # stores the values of the population
log_pop[1,] <- log(pop) 
w_T <- matrix(0, nrow = N, ncol = J) # stores the values of w_T

for (j in seq_len(J)){
  
  # rdm <- runif(2)
  # w_0 <- c(rdm[1],rdm[1]+rdm[2],1-(rdm[1]+rdm[2])) 
  w <- w_0
  for (t in seq_len(T)){
    
    #L <- sample_continuous_L(b_1, b_2, lambda_1, lambda_2, eps)
    L <- sample_discrete_L(ens_M, p)
    
    # print(v)
    w_new <- L %*% w
    log_pop[t+1,j] <- log(sum(w_new)) + log_pop[t,j]
    
    #print(w_new)
    
    if (t == T){
      ratio[j] <- log(sum(w_new))
      w_T[,j] <- w_new / sum(w_new)
    } else {
      w <- w_new / sum(w_new)
    }
  }
}

################################################################################
# Results

moy <- cumsum(ratio) / seq_len(J)
var <- cumsum((ratio - cumsum(ratio)/seq_len(J))^2) / c(1,seq_len(J-1))

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
