library(distr)

N = 1E5
lambda = sample(c(10,1), size = N, replace = TRUE, prob = c(.9,.1))
theta = sapply(lambda, function(l){
  D <- DExp(rate = l) 
  r(D)(1)}) #laplace
Y = sapply(theta, function(x) rnorm(1, x, 1))
p = 2*(1-pnorm(abs(Y)))
p1 = p.adjust(p, method = "fdr")
significant = which(p1<.2)
#things seem to match up for FDR p adjustement

plot(Y[significant], theta[significant])

#model 1- true prior


#model 2- uninformative prior
lines(Y[significant],Y[significant]+1.96)
lines(Y[significant],Y[significant]-1.96)


###how to find posterior? if it is joint?
D1 <- DExp(rate = 1) 
D10 <- DExp(rate = 10) 
likelihood<- function(Y, theta){
  (((.1)*d(D1)(theta)+(.9)*d(D10)(theta))*dnorm(Y-theta))
  
}
#metropolis hasting ish

ci <- NULL
for (y in Y[significant]){
  c = 1
  mu <- 0
  MU <- NULL
  
  for(s in 1:1000){
    mu.star <- rnorm(1, mu, c)
    r = likelihood(y,mu.star)/likelihood(y,mu)
    if(runif(1)<r){
      mu <- mu.star
    }
    MU <- c(MU, mu)
  }
  ci<- rbind(ci, HPD(MU))
  
}
hist(MU)
density(MU)

plot(Y[significant], theta[significant])

#model 1- true prior
for(i in 1:length(significant)){
  lines(c(Y[significant[i]],Y[significant[i]]), model3.ci[i,], pch = 3)
}
lines(Y[significant],ci[,1], pch = 3)


#model 2- uninformative prior
lines(Y[significant],Y[significant]+1.96)
lines(Y[significant],Y[significant]-1.96)



#model 3, spike and slab
getposterior <- function(Y,n.samp, pi = .5){
  n = length(Y)
  odds = (1-pi)/pi
  #bf <- -exp(1)*p*log(p)
  #could also do numerical integration to get BF for precision (but could be unstable)
  #look at BAS zellner
  #mathematica?? *integrate*
  ##plot priors
  bf <-  (n+1)^(-.5)*exp(n^2*mean(Y)^2/(2*(n+1)))
  alt.prob <- odds*bf/(1+odds*bf)
  
  #draws from posterior-flip a coin (ber w prob P(H given Y) and then use that to get draw
  draws = sapply(runif(n.samp), function(x)  {
    ifelse(x<(1-alt.prob), rnorm(1, 0, 0), rcauchy(1,mean(Y)*n/(n+1), sqrt(1/(n+1))))
  })
  #ggplot(data = data.frame(draws))+geom_density(aes(x=draws))
  (list(alt.prob=alt.prob, draws=draws)) 
  
}

model3.ci <- NULL
for (y in Y[significant]){
  posterior <- getposterior(y, 1000, pi = .1)
  model3.ci<- rbind(model3.ci, HPD(posterior$draws))
  
}
model3.ci


m1 = m2 = m3 = 0
for (i in 1:length(significant)){
  m1 = m1+ (ci[i,1]<=theta[significant[i]]& ci[i,2]>=theta[significant[i]])
  m2 = m2+(Y[significant[i]]-1.96<=theta[significant[i]]& Y[significant[i]]+1.96>=theta[significant[i]])
  m3 = m3+ (model3.ci[i,1]<=theta[significant[i]]& model3.ci[i,2]>=theta[significant[i]])
  }
m1/length(significant)
m2/length(significant)
m3/length(significant)
#m3 kind of terrible because of very intense shrinkage (n=1 so ybar/n+1 = ybar/2)
