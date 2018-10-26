# test3
N = 2
T = 2000

b = 1
c = 1
eps = 0.05


pop = 100

J = 1000

sample_L <- function(b, c, eps){
  
  L <- matrix(0,2,2)
  rdm <- runif(2, -eps, eps)
  L[1,2] = b + rdm[1]
  L[2,1] = c + rdm[2]
  
  return(L)
    
}


v_0 <- c(0.5,0.5)
ratio = numeric(J)
w_T = matrix(0, nrow = J, ncol = 2)

for (j in seq_len(J)){
  
  v <- v_0
  for (t in seq_len(T)){
    
    L <- sample_L(b,c,eps)
    
    if (t == T){
      
      ratio[j] = log(sum(L %*% v)) - log(sum(v))
      w_T[j,] = L %*% v / sum(L %*% v)
    }
    
    else v <- L %*% v

  }
  
}

m_J = mean(ratio)
var_J = 1/(J-1)*sum((ratio-mJ)^2)

moy = cumsum(ratio) / seq_len(J)
var = cumsum((ratio - cumsum(ratio)/seq_len(J))^2) / c(1,seq_len(J-1))
print(mean(ratio))
plot(seq_len(J), moy, type='l')
plot(seq_len(J), var, type='l')

hist(w_T[,1], breaks=100)
