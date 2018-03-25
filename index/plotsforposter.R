cond.likelihood<- function(ybar, mu, SE, q=1.96){
  dnorm(ybar, mu, SE)/
    (pnorm(-q*SE, mu, SE) +
       1-pnorm(q*SE, mu, SE))
}

p = 0.00325
LOR = log(1.65)
q = -qnorm(p/2)
# 1.21-2.25
SE =  (log(2.25 ) - LOR)/1.96
#exp(LOR + c(-1, 1)*1.96*SE)
BF = 1/(-exp(1)*p*log(p))
mu <- seq(-.2,1.25,.01)

#CL density

CLdata<-data.frame(mu, cond.likelihood(LOR, mu, SE, q = -qnorm(0.005/2)), 
                   cond.likelihood(LOR, mu, SE, q = -qnorm(0.05/2)),
                   cond.likelihood(LOR, mu, SE, q = -qnorm(0.01/2)),
                   dnorm(LOR, mu, SE))
colnames(CLdata)<- c("mu", "CL.005","CL.05","CL.01","uncond")
melted <- melt(CLdata, id.vars=c("mu"))
clplot<- ggplot(melted, aes(mu, value, color = variable)) + geom_line(aes(group=variable)) +
  labs(title = "Conditional vs. Unconditional Likelihood", 
       x = expression(theta), y = "Likelihood")+  theme_minimal()+ 
  theme(legend.justification=c(1,1), legend.position=c(1,1)) +
  scale_colour_few(name="", labels=c(expression(paste(alpha," = .005")), 
                                     expression(paste(alpha," = .05")), 
                                     expression(paste(alpha," = .01")),
                                     "unconditional"))


#####bayes approx
postind.dens<-function(pind, p=.0035){
  prior.odds<- (pind)/(1-pind) #a/0
  BF<- 1/(-exp(1)*p*log(p)) #eplogp is BF h0/h1
  post.ind<-prior.odds*BF/(1+prior.odds*BF)
  return(post.ind)
}
pind<- seq(.001, .999,by=.001)

BFapproxdata<-data.frame(pind, postind.dens(pind, .05),postind.dens(pind, .01),
                         postind.dens(pind, .005),postind.dens(pind, .0005),
                   dbeta(pind,.5,.5))
colnames(BFapproxdata)<-c("prior", "p = .05", "p = .01", "p = .005", "p = .0005","theta")
ggplot(data=melt(BFapproxdata,id.vars = "theta"), aes(y=theta))+
  geom_line(aes(x=value, group=variable,color=variable))+
  labs(title = "Posterior Distribution of Probability of Effect\n with Bayes Factor Approximation", 
       x = expression(theta), y = "Density")+  theme_minimal()+ 
  theme(legend.justification=c(1,1), legend.position=c(1,1)) +
  scale_colour_few(name="p-value")

#ggsave(filename="clplot.png", plot=clplot)

# Bayes density? a generic example
# set.seed(1100)
# samp<- c(rnorm(2000, .2, .1),rep(0,1000))
# p<-plotvar(samp, gg=TRUE)


bplot<-ggplot(data = data.frame(maxy = dnorm(.4, .4,.2)/.6, x = mu, y = dnorm(mu, .4,.2), xl=0, xu=0,yl=0,yu=.4))+
  geom_line(aes(x,y/maxy))+geom_segment(aes(x=xl,y=yl,xend=xu,yend=yu))+
  labs(title= "Mixture Model with Point Mass at 0",x=expression(theta),y="density")+  theme_minimal() + scale_colour_few()


ggsave(filename="twoplots.png",plot=grid.arrange(bplot,clplot))



# posterior data plots
 #mu confint
ciplot<-ggp+labs(title=expression(paste("Quantile-Based Confidence Intervals of ",mu[p53])),
       x="Model",y= "Odds Ratio")
post1<-data.frame(plotvar(exp(p53.simnew$BUGSoutput$sims.list$mu.p53),pointmass=1,gg=TRUE), model="Bayesian")
post2<-data.frame(plotvar(exp(p53.simnormal$BUGSoutput$sims.list$mu.p53),pointmass=1,gg=TRUE), model="Conditional Likelihood")
post3<-data.frame(plotvar(exp(p53.sim$BUGSoutput$sims.list$mu.p53),pointmass=1,gg=TRUE), model="Original")
postdf<- rbind(post1,post2,post3)
postplot<-ggplot(data.frame(postdf), aes(group=model))+
  geom_line(aes(x,y/maxy, color=model))+geom_segment(aes(x=xl,y=yl,xend=xu,yend=yu, color=model))+
  labs(x="Odds Ratio",y="Density",title=expression(paste("Posterior Density of ", mu[p53])))+  theme_minimal()+ 
  theme(legend.justification=c(1,1), legend.position=c(1,1)) + scale_colour_few()

ggsave(filename="dataplots.png",plot=grid.arrange(ciplot,postplot,))

#simulation plots
museplot<-plots1.c$plots[[1]]+  theme_minimal() + scale_colour_few()
muintplot<-plots1.c$plots[[2]]+  theme_minimal() + scale_colour_few()

rmse<-cbind(plots4.c$table[[1]][,2],plots1.c$table[[1]][,2])
rownames(rmse)<-c("Bayesian", "Original", "Conditional Likelihood")
colnames(rmse)<-c("diffuse true distribution", "dense true distribution")
mutable<-grid.arrange(tableGrob(rmse))

#mu4table<-grid.arrange(tableGrob(plots4.c$table[[1]]))
#plots0$plots[[1]] mayb enot


clseplot<-clplots1$plots[[1]]+  theme_minimal() +
  scale_colour_few() + 
clintplot<-clplots1$plots[[2]]+  theme_minimal() + scale_colour_few()
