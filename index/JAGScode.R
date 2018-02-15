#original model
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

#mixture with association variable
p53.newmodel3 = function() {
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
  mu.p53.1 ~ dnorm(0,.1)
  mu.p53<-mu.p53.1*assoc
  phi.p53 ~ dgamma(1, .05)
  #    phi.p53 ~ dgamma(2, .02)
  sigma.p53 <- pow(phi.p53, -.5)
  
  mu.Age ~ dnorm(0, .1)
  phi.Age ~ dgamma(1, .05)
  sigma.Age  <- pow(phi.Age, -.5)
  
  assoc ~ dbern(pind) 
  pind ~ dbeta(2,6)
}

#conditional likelihood prior (normal approx, using zeroes trick)
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
  tau<- pow(SE, -2)
  
  for (k in 1:n.discovery){
    #zeroes trick for MLE~cond prob
    L[k]<- dnorm(MLE[k],mu.p53, tau)/(pnorm(-q*SE, mu.p53, tau) +
                                        1-pnorm(q*SE, mu.p53, tau))
    phi[k]<- -log(L[k])+C
    zeroes[k]~dpois(phi[k])
  }
  
  beta.BC ~ dnorm(0, .1)
  
  mu.site ~ dnorm(0, .1)
  phi.site ~ dgamma(1,.05)
  sigma.site <- pow(phi.site, -.5)
  
  #E[prec] 20
  #(based on range .5 to 2 for OR => range = 1.4 = 6 sigma  sigma = 1.4/6 ~= .2 
  mu.p53 ~ dnorm(0,.1)
  phi.p53 ~ dgamma(1, .05)
  #    phi.p53 ~ dgamma(2, .02)
  sigma.p53 <- pow(phi.p53, -.5)
  
  mu.Age ~ dnorm(0, .1)
  phi.Age ~ dgamma(1, .05)
  sigma.Age  <- pow(phi.Age, -.5)
  #assoc~dbern(pind) *assoc
  pind ~ dbeta(2,6)
}

#mixture using zeroes trick (without assoc latent var)
p53.newmodel4 = function() {
  for (j in 1:J) {
    CaseCon[j] ~ dbern(theta[j])
    logit(theta[j]) <- beta.site[site[j]] + beta.p53[site[j]]*p53[j] +
      beta.Age[site[j]]*Age[j] + beta.BC*BC[j]    }
  
  for (l in 1:n.sites) {
    beta.site[l] ~ dnorm(mu.site, phi.site)
    beta.Age[l] ~ dnorm(mu.Age, phi.Age)
    beta.p53.1[l] ~ dnorm(mu.p53, phi.p53)
    beta.p53[l] <- beta.p53.1[l]*(1-mu.p53.iszero)
  }
  beta.BC ~ dnorm(0, .1)
  
  mu.site ~ dnorm(0, .1)
  phi.site ~ dgamma(1,.05)
  sigma.site <- pow(phi.site, -.5)
  
  #zeroes trick
  C<-1000
  epsilon<-0.001
  tau<- pow(epsilon,-2)
  #if mu geq 0 and mu leq 0, prior is from point mass
  #mu in (-epsilon, epsilon)
  mu.p53.iszero<- step(epsilon-abs(mu.p53))
  L<- (mu.p53.iszero*(1-pind))+ #(dnorm(mu.p53,0,tau)*(1-pind))+
    (dnorm(mu.p53,0,1)*pind)
  #need pind for mixing
  
  phi<- -log(L)+C
  zero~dpois(phi)
  mu.p53 ~ dunif(-2,2) #not sure how big this interval should be, just picked one large enough for .1 sd, started at 10, but even 1 has too large jumps
  
  phi.p53 ~ dgamma(1, .05)
  #    phi.p53 ~ dgamma(2, .02)
  sigma.p53 <- pow(phi.p53, -.5)
  
  mu.Age ~ dnorm(0, .1)
  phi.Age ~ dgamma(1, .05)
  sigma.Age  <- pow(phi.Age, -.5)
  #assoc~dbern(pind)
  pind ~ dbeta(2,10)
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
  C<-1000
  tau<- pow(SE, -2)
  beta.BC ~ dnorm(0, .1)
 
  mu.site ~ dnorm(0, .1)
  phi.site ~ dgamma(1,.05)
  sigma.site <- pow(phi.site, -.5)
  
  #E[prec] 20
  #(based on range .5 to 2 for OR => range = 1.4 = 6 sigma  sigma = 1.4/6 ~= .2 
  mu.p53.1 ~ dnorm(0,.1)
  mu.p53<-mu.p53.1*assoc
  phi.p53 ~ dgamma(1, .05)
  #    phi.p53 ~ dgamma(2, .02)
  sigma.p53 <- pow(phi.p53, -.5)
  
  mu.Age ~ dnorm(0, .1)
  phi.Age ~ dgamma(1, .05)
  sigma.Age  <- pow(phi.Age, -.5)
  
  assoc~dbern(post.ind)
  prior.odds<- (1-pind)/pind
  BF<- -exp(1)*p*log(p)
  post.ind<-prior.odds*BF/(1+prior.odds*BF)
  pind ~ dbeta(2,6)
}