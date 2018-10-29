library(rlist)

################################################################################
# Parameters

N <- 3 # number og age classes
M <- 3 # number of Leslie matrices

b_1 <- 1 # birth rate of class 1
b_2 <- 0.3 # birth rate of class 2
lambda_1 <- 0.8 # survival rate of class 0 (juvenile class)
lambda_2 <- 0.7 # survival rate of class 1

eps <- 0.2 # semi-range of parameters

T <- 5000 # number of time steps
pop <- 100 # initial number of individuals  
K <- 10^6 # maximal number of individuals  
J <- 1000 # number of Monte Carlo samples

################################################################################
# Useful functions

create_M  <- function(b, c, d, e, eps, M){
  
  # create set of Leslie matrices
  
  ens_M <- list()
  
  if (M == 1){
    
    L <- matrix(0,3,3)
    L[1,2] = b
    L[1,3] = c
    L[2,1] = min(1, d)
    L[3,2] = min(1, e)
    
    return(list.append(ens_M, L))
    
  }
  
  for (m in seq_len(M)){
    
      L <- matrix(0,3,3)
      L[1,2] = b - eps + (m-1) * 2 * eps / (M - 1)
      L[1,3] = c - eps + (m-1) * 2 * eps / (M - 1)
      L[2,1] = d - eps + min(1, (m-1) * 2 * eps / (M - 1))
      L[3,2] = e - eps + min(1, (m-1) * 2 * eps / (M - 1))
      
      ens_M <- list.append(ens_M,L)
  } 
  
  return(ens_M)

}


top_eigenvectors <- function(ens_M){
  
  M <- length(ens_M)
  top_values <- numeric(M)
  top_vectors <- matrix(0,N,M)
  
  for (i in seq_len(M)){
    
    decomp <- eigen(M[[i]])
    top_values[i] <- max(decomp$values)
    k <- which.max(decomp$values)
    top_vectors[i] <- decomp$vectors[,k]
    
  }
  
  return(list("values"=top_values, "vectors"=top_vectors))
  
}



sample_discrete_L <- function(ens_M, p=0){
  
  # sample from set M with probabilities p (uniformly if p=0)
  
  M <- length(ens_M)
  
  if (M == 1){
    
    return (ens_M[[1]])
    
  }
  
  if (p == 0 || length(p) < M){
    
    rdm <- runif(1)
    k <- as.integer(rdm * M) + 1
    return(ens_M[[k]])
    
  }

  rdm <- sample(x=seq_len(M), size=1, replace=T, prob=p)
  
  return(ens_M[[rdm]])
  
}



sample_continuous_L <- function(b, c, d, e, eps){
  
  # sample a Leslie matrix in a continuous set
  
  L <- matrix(0,3,3)
  rdm <- runif(4, -eps, eps)
  L[1,2] = b + rdm[1]
  L[1,3] = c + rdm[4]
  L[2,1] = min(1, d + rdm[2])
  L[3,2] = min(1, e + rdm[3])
  
  return(L)
    
}


################################################################################
# Simulations

v_0 <- pop * c(0.2,0.4,0.4) # initial state
p <- c(0.2, 0.7, 0.1) # probability distribution on M
ens_M <- create_M(b_1, b_2, lambda_1, lambda_2, eps)

ratio <- numeric(J) # stores the values of the estimated ratio
w_T <- matrix(0, nrow = J, ncol = N) # stores the values of w_T

for (j in seq_len(J)){
  
  v <- v_0
  for (t in seq_len(T)){
    
    # L <- sample_continuous_L(b_1, b_2, lambda_1, lambda_2, eps)
    L <- sample_discrete_L(ens_M, p)
    
    # print(v)
    
    v_new <- L %*% v
    if (sum(v_new) > K){
      
      v_new <- K * v_new / sum(v_new)
      
    }
    
    # print(v_new)
    
    if (t == T){
      
      ratio[j] = log(sum(v_new)) - log(sum(v))
      w_T[j,] = v_new / sum(v_new)
    }
    
    else v <- v_new
    
  }
  
}

moy = cumsum(ratio) / seq_len(J)
var = cumsum((ratio - cumsum(ratio)/seq_len(J))^2) / c(1,seq_len(J-1))
print(mean(ratio))
plot(seq_len(J), moy, type='l')
plot(seq_len(J), var, type='l')

plot(w_T[,2], w_T[,3])
