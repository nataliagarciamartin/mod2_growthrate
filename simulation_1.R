# test2
N = 3
T = 2000

b = 1.5
c = 0.9
d = 0.8
e = 0.5
eps = 0.3


pop = 100

J = 1000

sample_L <- function(b, c, d, e, eps){
  
  L <- matrix(0,3,3)
  rdm <- runif(4, -eps, eps)
  L[1,2] = b + rdm[1]
  L[1,3] = e + rdm[4]
  L[2,1] = c + rdm[2]
  L[3,2] = d + rdm[3]
  
  return(L)
    
}


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

m_J = mean(ratio)
var_J = 1/(J-1)*sum((ratio-mJ)^2)

moy = cumsum(ratio) / seq_len(J)
var = cumsum((ratio - cumsum(ratio)/seq_len(J))^2) / c(1,seq_len(J-1))
print(mean(ratio))
plot(seq_len(J), moy, type='l')
plot(seq_len(J), var, type='l')
