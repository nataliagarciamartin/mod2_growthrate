b_1 <- 1.5 # birth rate of class 1
b_2 <- 0.8 # birth rate of class 2
lambda_1 <- 0.9 # survival rate of class 0 (juvenile class)
lambda_2 <- 0.6 # survival rate of class 1

eps <- 0.3 # semi-range of parameters

################################################################################

run_markov <- function(w0, pop, ens_M, p, T){
  
  # Function running the Markov chain wt from w0 during T times steps with ens_M
  
  N <- length(w0)
  M <- length(ens_M)
  
  log_pop <- numeric(T+1) # stores the increment values of the log of the population
  
  log_pop[1] <- log(pop)
  
  w_t <- matrix(0, nrow = T+1, ncol = N) # stores the values of w_T
  
  w_t[1,] <- w_0
  w <- w_0
  for (t in seq_len(T)){
    
    #L <- sample_continuous_L(b_1, b_2, lambda_1, lambda_2, eps)
    L <- sample_discrete_L(ens_M, p)

    w_new <- L %*% w
    log_pop[t+1] <- log(sum(w_new))
    
    #print(w_new)
    
    w <- w_new / sum(w_new)
    w_t[t+1,] <- w
    
  }
  
  return(list("w"=w_t, "logp"=log_pop))
  
}