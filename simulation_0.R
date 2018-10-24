# test2
M = 2
N = 5
T = 200
b = 1
bprime = 0.7
d = 0.5
rho = 0.7

pop = 100

J = 1000

A_1 <- rbind(c(0,b,b,b,0), cbind(diag(1,N-1), rep(0,N-1)))
A_1[N,N-1] <- d

A_2 <- rbind(c(0,bprime,bprime,bprime,0), cbind(rho * diag(1,N-1), rep(0,N-1)))

v_0 <- pop * c(0.3, 0.2, 0.2, 0.2, 0.1)
ratio = numeric(J)

for (j in seq_len(J)){
  
  v <- v_0
  for (t in seq_len(T)){
    
    rdm = runif(1)
    
    if (rdm < 0.8){
      A = A_1
    }
    
    else A = A_2
    
    if (t == T){
      
      ratio[j] = log(sum(A %*% v)) - log(sum(v))
    }
    
    else v <- A %*% v
  }
  
}

print(mean(ratio))
plot(seq_len(J), cumsum(ratio) / seq_len(J), type='l')