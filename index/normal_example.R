##
### unconditional w selection (bma)
p<- 0.05
pi <- .5 #of H0
N <-100
n.samp <- 1000
prob <- .95

n.sim <- 100


##calculate posterior
getposterior <- function(Y,n.samp){
  n = length(Y)
  odds = (1-pi)/pi
  #bf <- -exp(1)*p*log(p)
  bf <-  (n+1)^(-.5)*exp(n^2*mean(Y)^2/(2*(n+1)))
  alt.prob <- odds*bf/(1+odds*bf)
  
  #draws from posterior-flip a coin (ber w prob P(H given Y) and then use that to get draw
  draws = sapply(runif(n.samp), function(x)  {
    ifelse(x<(1-alt.prob), rnorm(1, 0, 0), rnorm(1,mean(Y)*n/(n+1), sqrt(1/(n+1))))
  })
  #densityplot(draws)
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
  #stieljes integral????
  #look into this and check about continuity of cdf, etc
}


##confidence from cond likelihood

conf.int<- function(Y, p, prob, n.samp){
  cond.likelihood<- function(mu){
    dnorm(mean(Y), mu, sqrt(1/N))/(pnorm(-1.96+mu)+pnorm(-1.96-mu)) 
  }
  #metropolis hasting with flat prior
  S=10000; c=1
  mu <- 0
  MU <- NULL
  
  for(s in 1:S){
    mu.star <- rnorm(1, mu, c)
    r = cond.likelihood(mu.star)/cond.likelihood(mu)
    if(runif(1)<r){
      mu <- mu.star
    }
    MU <- c(MU, mu)
  }
  #densityplot(MU)
  (HPD(MU))
}




##simulations

  all<- function(H){
  
  mu = ifelse(H==0, 0, rnorm(1, 1, 1)) ### normally distributed mu
  Y <- rnorm(N, mu, 1) 
  post <- getposterior(Y,  n.samp)
  alt.prob = post$alt.prob
  cred <- HPD(post$draws)
  cred.lower = cred[1]
  cred.upper = cred[2]
  bayes.cov =  (cred.upper>mu&&cred.lower<mu)
  
  accepted <-0
  conf.lower<-conf.upper<-freq.cov<-expected.cov<-NA
  if(abs((mean(Y))/sqrt(1/N)-1.96)<.2){ #if Z is in (1.94, 1.98)
    accepted = 1
    conf <- conf.int(Y, p, prob, n.samp)
    conf.lower = conf[1]
    conf.upper = conf[2]
    
    freq.cov =  (conf.upper>mu&&conf.lower<mu)
    expected.cov <- .95*alt.prob+(conf.upper>0&&conf.lower<0)*(1-alt.prob)
    
  }
  
  (c(H, mu , mean(Y), alt.prob , accepted , 
              cred.lower  , cred.upper  , 
              conf.lower , conf.upper , 
              bayes.cov , freq.cov , expected.cov ))
  }

results =data.frame(t(apply(matrix(as.numeric(runif(n.sim)<pi)),1,all)))
colnames(results)<- c("H", "mu" , "mean(Y)", "alt.prob" , "accepted" , 
                      "cred.lower"  , "cred.upper"  , 
                      "conf.lower" , "conf.upper" , 
                      "bayes.cov" , "freq.cov" , "expected.cov" )
