# File containing useful functions

library(rlist)


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


top_eigen <- function(ens_M){
  
  M <- length(ens_M)
  top_values <- numeric(M)
  top_vectors <- matrix(0,N,M)
  
  for (i in seq_len(M)){
    
    decomp <- eigen(ens_M[[i]])
    
    if (any(sapply(decomp$values, is.complex))){
      top_values[i] <- Re(decomp$values[1])
      top_vectors[,i] <- Re(decomp$vectors[,1])
    } else {
      top_values[i] <- max(decomp$values)
      k <- which.max(decomp$values)
      top_vectors[,i] <- decomp$vectors[,k]
    }
  }
  
  return(list("values"=top_values, "vectors"=top_vectors))
  
}



sample_discrete_L <- function(ens_M, p=0){
  
  # sample from set M with probabilities p (uniformly if p=0)
  
  M <- length(ens_M)
  
  if (M == 1){
    
    return (ens_M[[1]])
    
  }
  
  if (p == 0 || length(p) != M){
    
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