# initial population vector
v0 = c(0.2,0.8)

# projection matrix and its characteristics
A = matrix(c(0.4,0.7,0.8571429,0),2,2)
ev1 = eigen(A)
ev1$values
v = ev1$vectors

# initialising the plot
plot(NULL,main=TeX('Stationary distribution of $\\omega_t$'), xlab=TeX('$\\omega_1$'),ylab=TeX('$\\omega_2$'), xlim=c(0,1), ylim=c(0,1))
arrows(0,0,v[1,1],v[2,1],col='red')

vec = v0
for(i in seq(n-1)){
  
  vec = A%*%vec
  vec = vec/sum(vec)
  plot(NULL,main=TeX('Stationary distribution of $\\omega_t$'),xlab=TeX('$\\omega_1$'),ylab=TeX('$\\omega_2$'), xlim=c(0,1), ylim=c(0,1))
  arrows(0,0,v[1,1],v[2,1],col='red')
  arrows(0,0,vec[1,1],vec[2,1],col='black')
  text(0.5, 0.9, paste("Iteration number:",  i))
  pause(1)
  
}