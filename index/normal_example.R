##
### unconditional w selection (bma)
#p<- 0.05
#pi <- .5 #of H0
#N <-10
#n.samp <- 1000
#prob <- .95

#n.sim <- 100


##calculate posterior
getposterior <- function(Y,n.samp, pi = .5){
  n = length(Y)
  odds = (1-pi)/pi
  #bf <- -exp(1)*p*log(p)
  bf <-  (n+1)^(-.5)*exp(n^2*mean(Y)^2/(2*(n+1)))
  alt.prob <- odds*bf/(1+odds*bf)
  
  #draws from posterior-flip a coin (ber w prob P(H given Y) and then use that to get draw
  draws = sapply(runif(n.samp), function(x)  {
    ifelse(x<(1-alt.prob), rnorm(1, 0, 0), rnorm(1,mean(Y)*n/(n+1), sqrt(1/(n+1))))
  })
  #ggplot(data = data.frame(draws))+geom_density(aes(x=draws))
  (list(alt.prob=alt.prob, draws=draws)) 
  
}


HPD <-function(post, prob  = .95){
  
  #HPD interval- copied from BAS/coda
  obj <- as.matrix(post)
  vals <- apply(obj, 2, sort)
  if (!is.matrix(vals))
    stop("obj must have nsamp > 1")
  nsamp <- nrow(vals)
  npar <- ncol(vals)
  gap <- max(1, min(nsamp - 1, round(nsamp * prob)))
  init <- 1:(nsamp - gap)
  inds <- apply(vals[init + gap, , drop = FALSE] - vals[init,
                                                        , drop = FALSE], 2, which.min)
  (cbind(vals[cbind(inds, 1:npar)], vals[cbind(inds +
                                                 gap, 1:npar)]))
  #look into this and check about continuity of cdf, etc
}


##confidence from cond likelihood

cond.posterior<- function(Y, n.samp){
  cond.likelihood<- function(Y, mu){
    ybar = mean(Y)
    N = length(Y)
    #(abs(ybar-mu)/(sqrt(1/N))>1.96)*
    dnorm(ybar, mu, sqrt(1/N))/
      (pnorm((-1.96*sqrt(1/N)), mu, sqrt(1/N))
       +1-pnorm((1.96*sqrt(1/N)), mu, sqrt(1/N)))
    #not symmetric other than mu = 0
  }
  #x<- seq(-2,2,.01)
  #plot(x, cond.likelihood(Y, x), type="l")
  #metropolis hasting with flat prior
  c=1
  mu <- 0
  MU <- NULL
  
  for(s in 1:n.samp){
    mu.star <- rnorm(1, mu, c)
    r = cond.likelihood(Y,mu.star)/cond.likelihood(Y,mu)
    if(runif(1)<r){
      mu <- mu.star
    }
    MU <- c(MU, mu)
  }
  #plot(MU, type="l")
  #ggplot(data = data.frame(MU))+geom_density(aes(x=MU))
  (MU)
}




##simulations
all <- function(H, N, n.samp, interval,alpha){
mu <- ifelse(H==0, 0, rnorm(1, 0, 1)) ### normally distributed mu
Y <- rnorm(N, mu, 1)
count=0
while(abs(abs(mean(Y))/sqrt(1/N)-qnorm(1-alpha/2))>interval){ #if Z is not in (1.94, 1.98)
  if(count>1000){ #had to add this bc it wouldn't run
    return (c(H, mu , mean(Y), NA  , 
              NA  , NA  , 
              NA , NA , 
              NA , NA , NA ))
  }
  
  Y <- rnorm(N, mu, 1)
count<- count+1
}
post <- getposterior(Y,  n.samp)
alt.prob = post$alt.prob
cred <- HPD(post$draws)
cred.lower = cred[1]
cred.upper = cred[2]
bayes.cov =  (cred.upper>=mu&&cred.lower<=mu)

cond <-cond.posterior(Y, n.samp)
conf <- HPD(cond)
conf.lower = conf[1]
conf.upper = conf[2]

freq.cov =  (conf.upper>=mu&&conf.lower<=mu)
expected.cov <- .95*alt.prob+(conf.upper>=0&&conf.lower<=0)*(1-alt.prob)


(c(H, mu , mean(Y), alt.prob  , 
   cred.lower  , cred.upper  , 
   conf.lower , conf.upper , 
   bayes.cov , freq.cov , expected.cov ))

}


N=100; n.samp = 10000;n.sim=1000; pi = .5

results =data.frame(t(apply(matrix(as.numeric(runif(n.sim)<pi)),1, function(x) all(x, N, n.samp, .2, .05))))
colnames(results)<- c("H", "mu" , "Ybar", "alt.prob" , 
                     "cred.lower"  , "cred.upper"  , 
                     "conf.lower" , "conf.upper" , 
                     "bayes.cov" , "freq.cov" , "expected.cov" )


#for(i in as.numeric(runif(10)<pi)){print(i);print(all(i, N, n.samp))}

dim(results)

#conf int w 0
summary(results)
a <- which(results$conf.upper>=0&results$conf.lower<=0)
summary(results[a,])
summary(results[-a,])

# true nulls
b <- which(results$H==0)
summary(results[b,])
summary(results[-b,])
ggplot(data = na.omit(results))+geom_jitter(aes(y = mu, x= Ybar))+geom_violin(aes(y = mu, x = Ybar, alpha = .2))
#everything contains 0!!!!!!!- bayes factor/alt probability are too low in this range
#calculated coverage and true coverage are not close bc of selection?

table(results$conf.upper>=0&results$conf.lower<=0,results$H==0)
table(results$cred.upper>=0&results$cred.lower<=0,results$H==0)

#hard to get the stuff in the interval
#really bad for H1 for both
#Bayesian shrinkage actually makes things worse

ggplot(data = results)+ geom_point(aes(x = conf.upper-conf.lower, y = cred.upper-cred.lower, colour = H))+geom_abline(intercept = 0, slope = 1)
#almost all credible intervals are smaller

results2 =data.frame(t(apply(matrix(as.numeric(runif(n.sim)<pi)),1, function(x) all(x, N, n.samp, .2, .01))))
colnames(results2)<- c("H", "mu" , "Ybar", "alt.prob" , 
                      "cred.lower"  , "cred.upper"  , 
                      "conf.lower" , "conf.upper" , 
                      "bayes.cov" , "freq.cov" , "expected.cov" )




results1[8,]
results2[33,]
