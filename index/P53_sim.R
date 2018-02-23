set.seed(1)
#fix everything (even betas)
mu.p53=.5
phi.p53= 10
observations=100
n.sites=7
beta.p53.fixed = rnorm(n.sites,mu.p53,phi.p53^(-.5)) #"fixed" because seed
Y <-site <- rep(NA, observations*n.sites)
for(i in 1:n.sites){
    Y[((i-1)*observations+1):(i*observations)]<-rbinom(observations, 1, exp(beta.p53[i])/(1+exp(beta.p53[i])))
    site[((i-1)*observations+1):(i*observations)]<- rep(i, observations)
}
sim.data1.1<-list(CaseCon=Y, site=site,  J=n.sites*observations ,n.sites=n.sites)

beta.p53 = rep(0,n.sites) #"fixed" because seed
Y <-site <- rep(NA, observations*n.sites)
for(i in 1:n.sites){
  Y[((i-1)*observations+1):(i*observations)]<-rbinom(observations, 1, exp(beta.p53[i])/(1+exp(beta.p53[i])))
  site[((i-1)*observations+1):(i*observations)]<- rep(i, observations)
}
sim.data1.0<-list(CaseCon=Y, site=site,  J=n.sites*observations ,n.sites=n.sites)
#this might just be 0?

#fix mu, assoc
S=10
Y <-site <- rep(NA, observations*n.sites*S)
for(j in 1:S){
  beta.p53 = rnorm(n.sites,mu.p53,phi.p53^(-.5)) #"fixed" because seed
  for(i in 1:n.sites){
    Y[((i-1)*observations+1):(i*observations)]<-rbinom(observations, 1, exp(beta.p53[i])/(1+exp(beta.p53[i])))
    site[((i-1)*observations+1):(i*observations)]<- rep(i, observations)
  }
}
sim.data2.1<-list(CaseCon=Y, site=site,  J=n.sites*observations*S ,n.sites=n.sites)

#fix assoc

#more than one dataset

#model
p53.test.full = function() {
  for (j in 1:J) {
    CaseCon[j] ~ dbern(theta[j])
    logit(theta[j]) <-  beta.p53[site[j]]  }
  
  for (l in 1:n.sites) {
    beta.p53.1[l] ~ dnorm(mu.p53, phi.p53)
    beta.p53[l] <- beta.p53.1[l]*(1-mu.p53.iszero)
  }
  #zeroes trick
  C<-1000
  epsilon<-0.01
  tau<- pow(epsilon,-2)
  #if mu geq 0 and mu leq 0, prior is from point mass
  #mu in (-epsilon, epsilon)
  mu.p53.iszero<- step(epsilon-abs(mu.p53))
  L<- (dnorm(mu.p53,0,tau)*(1-pind))+ #(mu.p53.iszero*(1-pind))+ #
    (dnorm(mu.p53,0,1)*pind)
  #need pind for mixing
  
  phi<- -log(L)+C
  zero~dpois(phi)
  mu.p53 ~ dunif(-10,10) #not sure how big this interval should be, just picked one large enough for .1 sd, started at 10, but even 1 has too large jumps
  
  phi.p53 ~ dgamma(1, .05)
  #    phi.p53 ~ dgamma(2, .02)
  sigma.p53 <- pow(phi.p53, -.5)
  
  #assoc~dbern(pind)
  pind ~ dbeta(2,10)
}

#hierarchical model
p53.fulljags1.1 = jags(data=sim.data1.1, ####this
                    inits=NULL, parameters.to.save =c("pind", "mu.p53", "phi.p53","mu.p53.iszero"),
                    model = p53.test.full)
par(mar=rep(1,4))
plot(as.mcmc(p53.fulljags1.1))
HPDM(p53.fulljags1.1$BUGSoutput$sims.list$mu.p53)
mean(p53.fulljags1.1$BUGSoutput$sims.list$mu.p53)
median(p53.fulljags1.1$BUGSoutput$sims.list$mu.p53)
#glm
p53.fullglm1.1 = glm(CaseCon ~ factor(site), data=sim.data1.1,  ####also this
               family=binomial, x=T)

coef(p53.fullglm1.1)
#some observations: does terribly for mu=0, p good for true mu
#pind stays the same no matter what

#glm also doesn't get the intercepts???







