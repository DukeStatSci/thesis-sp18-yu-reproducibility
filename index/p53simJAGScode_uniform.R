cond.likelihood.model = function() {
  for (j in 1:J) {
    CaseCon[j] ~ dbern(theta[j])
    logit(theta[j]) <-  mu.p53  }
  
  C<-1000
  for (k in 1:n.discovery){
    #zeroes trick for MLE~cond prob
    tau[k]<- pow(SE[k], -2)
    L[k]<- dnorm(MLE[k],mu.p53, tau[k])/(pnorm(-q*SE, mu.p53, tau[k]) +
                                           1-pnorm(q*SE, mu.p53, tau[k]))
    phi[k]<- -log(L[k])+C
    zeroes[k]~dpois(phi[k])
  }
  
  mu.p53 ~ dnorm(0,.1)
  
}

cond.likelihood.re.model = function() {
  for (j in 1:J) {
    CaseCon[j] ~ dbern(theta[j])
    logit(theta[j]) <-  beta.p53[site[j]]  }
  
  C<-1000
  for (k in 1:n.discovery){
    #zeroes trick for MLE~cond prob
    tau[k]<- pow(SE[k], -2)
    L[k]<- dnorm(MLE[k],beta.p53[discovery.sites[k]], tau[k])/
      (pnorm(-q*SE, beta.p53[discovery.sites[k]], tau[k]) +
         1-pnorm(q*SE, beta.p53[discovery.sites[k]], tau[k]))
    phi[k]<- -log(L[k])+C
    zeroes[k]~dpois(phi[k])
  }
  for (i in 1:n.sites){ #total sites
    beta.p53[i] ~dnorm(mu.p53,phi.p53)
  }
  
  mu.p53 ~ dnorm(0,.1)
  # phi.p53 ~ dgamma(1, .05)
  # sigma.p53 <- pow(phi.p53, -.5)
  phi.p53 <- pow(sigma.p53, -2)
  sigma.p53 ~ dunif(0,1) #dt(0,1,1)%_%T(0,)
}

ones.cauchy.model= function() {
  for (j in 1:J) {
    CaseCon[j] ~ dbern(theta[j])
    logit(theta[j]) <-  beta.p53[site[j]]  }
  
  for (l in 1:n.sites) {
    beta.p53.1[l] ~ dnorm(mu.p53, phi.p53)
    beta.p53[l] <- beta.p53.1[l]*(mu.p53.notzero)
  }
  #ones trick
  C<-1e6
  epsilon<-0.01
  tau<- pow(epsilon,-2)
  mu.p53.notzero<- step(temp)
  temp<-dt(mu.p53,0,1, 1)-dnorm(mu.p53,0,tau) 
  #cauchy is same as t with df=1
  L<-(dnorm(mu.p53,0,tau)*(1-pind))+(dt(mu.p53,0,1,1)*pind)
  p <- L/ C
  one~ dbern(p)
  mu.p53 ~ dunif(-10,10) 
  # phi.p53 ~ dgamma(1, .05)
  # sigma.p53 <- pow(phi.p53, -.5)
  phi.p53 <- pow(sigma.p53, -2)
  sigma.p53 ~ dunif(0,1) #dt(0,1,1)%_%T(0,)
  pind ~ dbeta(.5,.5)
}

ones.normal.model= function() {
  for (j in 1:J) {
    CaseCon[j] ~ dbern(theta[j])
    logit(theta[j]) <-  beta.p53[site[j]]  }
  
  for (l in 1:n.sites) {
    beta.p53.1[l] ~ dnorm(mu.p53, phi.p53)
    beta.p53[l] <- beta.p53.1[l]*(mu.p53.notzero)
  }
  #ones trick
  C<-1e6
  epsilon<-0.01
  tau<- pow(epsilon,-2)
  mu.p53.notzero<- step(temp)
  temp<-dnorm(mu.p53,0,1)-dnorm(mu.p53,0,tau) 
  L<-(dnorm(mu.p53,0,tau)*(1-pind))+(dnorm(mu.p53,0,1)*pind)
  p <- L/ C
  one ~ dbern(p)
  mu.p53 ~ dunif(-10,10) 
  # phi.p53 ~ dgamma(1, .05)
  # sigma.p53 <- pow(phi.p53, -.5)
  phi.p53 <- pow(sigma.p53, -2)
  sigma.p53 ~ dunif(0,1) #dt(0,1,1)%_%T(0,)
  pind ~ dbeta(.5,.5)
}

latent.normal.model= function() {
  for (j in 1:J) {
    CaseCon[j] ~ dbern(theta[j])
    logit(theta[j]) <-  beta.p53[site[j]]  }
  
  for (l in 1:n.sites) {
    beta.p53.1[l] ~ dnorm(mu.p53, phi.p53)
    beta.p53[l] <- beta.p53.1[l]*(mu.p53.notzero)
  }
  
  mu.p53<- mu1.p53*mu.p53.notzero
  mu1.p53 ~ dnorm(0,1)
  # phi.p53 ~ dgamma(1, .05)
  # sigma.p53 <- pow(phi.p53, -.5)
  phi.p53 <- pow(sigma.p53, -2)
  sigma.p53 ~ dunif(0,1) #dt(0,1,1)%_%T(0,)
  mu.p53.notzero~dbern(pind)
  pind ~ dbeta(.5,.5)
}

latent.cauchy.model= function() {
  for (j in 1:J) {
    CaseCon[j] ~ dbern(theta[j])
    logit(theta[j]) <-  beta.p53[site[j]]  }
  
  for (l in 1:n.sites) {
    beta.p53.1[l] ~ dnorm(mu.p53, phi.p53)
    beta.p53[l] <- beta.p53.1[l]*(mu.p53.notzero)
  }
  
  mu.p53<- mu1.p53*mu.p53.notzero
  mu1.p53 ~ dt(0,1, 1)
  # phi.p53 ~ dgamma(1, .05)
  # sigma.p53 <- pow(phi.p53, -.5) #half cauchy, uniform to 1 or to 5, take a guess on sd
  phi.p53 <- pow(sigma.p53, -2)
  sigma.p53 ~ dunif(0,1) #dt(0,1,1)%_%T(0,)
  mu.p53.notzero~dbern(pind)
  pind ~ dbeta(.5,.5)
}

fixed.model = function() {
  for (j in 1:J) {
    CaseCon[j] ~ dbern(theta[j])
    logit(theta[j]) <-  mu.p53  
  }
  mu.p53<- mu1.p53*mu.p53.notzero
  mu1.p53 ~ dnorm(0,.1)
  mu.p53.notzero~dbern(pind)
  pind ~ dbeta(.5,.5)
}

original.model = function() {
  for (j in 1:J) {
    CaseCon[j] ~ dbern(theta[j])
    logit(theta[j]) <-  beta.p53[site[j]]  }
  
  for (l in 1:n.sites) {
    beta.p53[l] ~ dnorm(mu.p53, phi.p53)
  }
  
  mu.p53 ~ dnorm(0,1)
  # phi.p53 ~ dgamma(1, .05)
  # sigma.p53 <- pow(phi.p53, -.5)
  phi.p53 <- pow(sigma.p53, -2)
  sigma.p53 ~ dunif(0,1) #dt(0,1,1)%_%T(0,)
}