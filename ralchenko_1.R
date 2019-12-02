library(MASS)

Wien = function(num, delta){
  rn<-rnorm(n=num, m=0, sd=1)
  c(0, cumsum(rn))*sqrt(delta)}

Euler_est = function(delta, t_max, t_0=0){
  t = seq(t_0, t_max, by = delta)                    
  n = length(t)
  x = c(x0)
  w = Wien(n, delta)
  for(k in 1:(n-1))
    x<-c(x, x[k] + a(t[k], x[k])*delta + b(t[k],x[k])*(w[k+1]-w[k]))
  return(x)
}

x0 = -2.5
delta = 0.01 #0.001
a = function(t, x) atan(2 + x*t)
b = function(t, x) log(1 + t + x*x)
t = seq(0, 1, by = delta)                    
n = length(t)

x =  matrix(nrow=10, ncol=n)
for(i in 1:10)
  x[i,] = Euler_est(delta, 1)

matplot(t, t(x[1:10,]), col=c(2:11), type = "l", main = "Euler estimation", xlab = "delta = 0.01", ylab = "")
