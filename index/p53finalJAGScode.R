p53.model = function() {
  for (j in 1:J) {
    CaseCon[j] ~ dbern(theta[j])
    logit(theta[j]) <- beta.site[site[j]] + beta.p53[site[j]]*p53[j] +
      beta.Age[site[j]]*Age[j] + beta.BC*BC[j]
  }
  for (k in 1:n.sites) {
    beta.site[k] ~ dnorm(mu.site, phi.site)
    beta.p53[k] ~ dnorm(mu.p53, phi.p53)
    beta.Age[k] ~ dnorm(mu.Age, phi.Age)
  }
  
  beta.BC ~ dnorm(0, 3)
  
  mu.site ~ dnorm(0, .1)
  phi.site <- pow(sigma.site, -2)
  sigma.site ~ dunif(0,5)
  
  mu.p53 ~ dnorm(0, .1)
  phi.p53 <- pow(sigma.p53, -2)
  sigma.p53 ~ dunif(0, 5)
  
  mu.Age ~ dnorm(0, .1)
  phi.Age <- pow(sigma.Age, -2)
  sigma.Age  ~ dunif(0, 5)
  
}

p53.newmodel = function() {
  for (j in 1:J) {
    CaseCon[j] ~ dbern(theta[j])
    logit(theta[j]) <- beta.site[site[j]] + beta.p53[site[j]]*p53[j] +
      beta.Age[site[j]]*Age[j] + beta.BC*BC[j]    }
  
  for (l in 1:n.sites) {
    beta.site[l] ~ dnorm(mu.site, phi.site)
    beta.p53.1[l] ~ dnorm(mu.p53, phi.p53)
    beta.p53[l] <- beta.p53.1[l]*assoc
    beta.Age[l] ~ dnorm(mu.Age, phi.Age)
  }
  beta.BC ~ dnorm(0, .1)
  
  mu.site ~ dnorm(0, .1)
  phi.site ~ dgamma(1,.05)
  sigma.site <- pow(phi.site, -.5)
  
  #E[prec] 20
  #(based on range .5 to 2 for OR => range = 1.4 = 6 sigma  sigma = 1.4/6 ~= .2 
  mu.p53<- mu1.p53*assoc
  mu1.p53 ~ dt(0,.1, 1)
  phi.p53 <- pow(sigma.p53, -2)
  sigma.p53 ~ dt(0,1,1)%_%T(0,)
  
  mu.Age ~ dnorm(0, .1)
  phi.Age ~ dgamma(1, .05)
  sigma.Age  <- pow(phi.Age, -.5)
  
  assoc ~ dbern(pind) 
  pind ~ dbeta(.5,.5)
}


#mixed effect conditional likelihood
p53.normal = function() {
  for (j in 1:J) {
    CaseCon[j] ~ dbern(theta[j])
    logit(theta[j]) <- beta.site[site[j]] + beta.p53[site[j]]*p53[j] +
      beta.Age[site[j]]*Age[j] + beta.BC*BC[j]    }
  
  for (l in 1:n.sites) {
    beta.site[l] ~ dnorm(mu.site, phi.site)
    beta.p53[l] ~ dnorm(mu.p53, phi.p53)
    beta.Age[l] ~ dnorm(mu.Age, phi.Age)
  }
  C<-1000
  
  for(m in (n.sites+1):(n.sites+n.discovery)){
    beta.p53[m] ~ dnorm(mu.p53, phi.p53)
  }
  
  for (k in 1:n.discovery){
    #zeroes trick for MLE~cond prob
    tau[k]<- pow(SE[k], -2)
    L[k]<- dnorm(MLE[k],beta.p53[k+n.sites], tau[k])/
      (pnorm(-q*SE[k], beta.p53[k+n.sites], tau[k]) +
         1-pnorm(q*SE[k], beta.p53[k+n.sites], tau[k]))
    phi[k]<- -log(L[k])+C
    zeroes[k]~dpois(phi[k])
  }
  
  beta.BC ~ dnorm(0, .1)
  
  mu.site ~ dnorm(0, .1)
  phi.site ~ dgamma(1,.05)
  sigma.site <- pow(phi.site, -.5)
  
  #E[prec] 20
  #(based on range .5 to 2 for OR => range = 1.4 = 6 sigma  sigma = 1.4/6 ~= .2 
  mu.p53 ~ dt(0,.1,1)
  phi.p53 <- pow(sigma.p53, -2)
  sigma.p53 ~ dt(0,1,1)%_%T(0,)
  
  mu.Age ~ dnorm(0, .1)
  phi.Age ~ dgamma(1, .05)
  sigma.Age  <- pow(phi.Age, -.5)
  #assoc~dbern(pind) *assoc
  pind ~ dbeta(.5,.5)
}


#eplogp approximation of BF for posterior prob of association
p53.bf.approx = function() {
  for (j in 1:J) {
    CaseCon[j] ~ dbern(theta[j])
    logit(theta[j]) <- beta.site[site[j]] + beta.p53[site[j]]*p53[j] +
      beta.Age[site[j]]*Age[j] + beta.BC*BC[j]    }
  
  for (l in 1:n.sites) {
    beta.site[l] ~ dnorm(mu.site, phi.site)
    beta.p53.1[l] ~ dnorm(mu.p53, phi.p53)
    beta.p53[l] <- beta.p53.1[l]*(assoc)
    beta.Age[l] ~ dnorm(mu.Age, phi.Age)
  }
  beta.BC ~ dnorm(0, .1)
  mu.site ~ dnorm(0, .1)
  phi.site ~ dgamma(1,.05)
  sigma.site <- pow(phi.site, -.5)
  
  #E[prec] 20
  #(based on range .5 to 2 for OR => range = 1.4 = 6 sigma  sigma = 1.4/6 ~= .2 
  mu.p53.1 ~ dt(0,.1,1)
  mu.p53<-mu.p53.1*assoc
  phi.p53 <- pow(sigma.p53, -2)
  sigma.p53 ~ dt(0,1,1)%_%T(0,)
  
  mu.Age ~ dnorm(0, .1)
  phi.Age ~ dgamma(1, .05)
  sigma.Age  <- pow(phi.Age, -.5)
  
  assoc~dbern(post.ind)
  prior.odds<- (pind)/(1-pind) #a/0
  BF<- 1/(-exp(1)*p*log(p)) #eplogp is BF h0/h1
  post.ind<-prior.odds*BF/(1+prior.odds*BF)
  pind ~ dbeta(.9,.9)
}
