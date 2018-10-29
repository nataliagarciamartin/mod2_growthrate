library(rlist)
source("aux_functions.R")

################################################################################
# Parameters

N <- 3 # number og age classes
M <- 1 # number of Leslie matrices

b_1 <- 1.1 # birth rate of class 1
b_2 <- 0.8 # birth rate of class 2
lambda_1 <- 0.9 # survival rate of class 0 (juvenile class)
lambda_2 <- 0.9 # survival rate of class 1

eps <- 0.1 # semi-range of parameters

T <- 1000 # number of time steps
pop <- 1000 # initial number of individuals  
K <- 10^6 # maximal number of individuals  
J <- 500 # number of Monte Carlo samples

################################################################################
# Simulations

v_0 <- pop * c(0.2,0.4,0.4) # initial state
p <- c(0.2, 0.7, 0.1) # probability distribution on M
ens_M <- create_M(b_1, b_2, lambda_1, lambda_2, eps, M)

ratio <- numeric(J) # stores the values of the estimated ratio
pop_evo <- matrix(0, nrow = T+1, ncol = J) # stores the values of the population
pop_evo[1,] <- pop 
w_T <- matrix(0, nrow = J, ncol = N) # stores the values of w_T

for (j in seq_len(J)){
  
  v <- v_0
  for (t in seq_len(T)){
    
    # L <- sample_continuous_L(b_1, b_2, lambda_1, lambda_2, eps)
    L <- sample_discrete_L(ens_M, p)
    
    # print(v)
    v_new <- L %*% v
    
    # print(v_new)
    
    if (t == T){
      ratio[j] <- log(sum(v_new)) - log(sum(v))
      w_T[j,] <- v_new / sum(v_new)
      pop_evo[t+1,j] <- sum(v_new)
    } else {
      
      if (sum(v_new) > K){
        v_new <- K * v_new / sum(v_new)
      }
      pop_evo[t+1,j] <- sum(v_new)
      v <- v_new
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
} else a_true <- log(max(eig$values))

eig_vectors <- top_eigen(ens_M)$vectors

################################################################################
# Plots

plot(seq_len(T+1), pop_evo[,1], type='l')
#print(mean(ratio))
#print(a_true)

plot(seq_len(J), moy, type='l')
#abline(h=a_true, col="blue")
plot(seq_len(J), var, type='l')

plot(w_T[,2], w_T[,3])
abline(a=eig_vectors[3,1]/eig_vectors[2,1], b=0.240, col="blue")
