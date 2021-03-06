
####not sure if this is necessary
Modes <- function(x, min.size) {
  ### Initial Checks
  if(missing(x)) stop("The x argument is required.")
  x <- as.vector(as.numeric(as.character(x)))
  x <- x[is.finite(x)]
  ### Amodal
  if(sd(x)==0)
    return(list(modes=NA, mode.dens=NA, size=1))
  ### Differentiate kernel density by x
  length(density(x)$y)
  dens.y.diff <- density(x)$y[-1] - density(x)$y[-length(density(x)$y)]
  incr <- dens.y.diff
  incr[which(dens.y.diff > 0)] <- 1
  incr[which(dens.y.diff <= 0)] <- 0
  ### Kernel density by increasing/decreasing density regions
  begin <- 1; count <- 1
  for (i in 2:length(incr)) {
    if(incr[i] != incr[i-1]) {
      count <- count + 1
      begin <- c(begin, i)}
  }
  begin <- c(begin, length(incr))
  size <- modes <- mode.dens <- rep(0, count/2)
  init <- 1
  dens <- density(x); sumdens <- sum(dens$y)
  if(incr[1] == 0) {
    size[1] <- sum(dens$y[1:begin[2]]) / sumdens
    init <- 2}
  j <- init
  for (i in init:length(size)) {
    size[i] <- sum(dens$y[begin[j]:begin[j+2]]) / sumdens
    kde <- dens
    kde$x <- kde$x[begin[j]:begin[j+2]]
    kde$y <- kde$y[begin[j]:begin[j+2]]
    modes[i] <- kde$x[kde$y == max(kde$y)]
    mode.dens[i] <- kde$y[kde$y == max(kde$y)]
    j <- j + 2
  }
  ### Order everything by density
  size <- size[order(mode.dens, decreasing=TRUE)]
  modes <- modes[order(mode.dens, decreasing=TRUE)]
  mode.dens <- mode.dens[order(mode.dens, decreasing=TRUE)]
  ### Remove modes with size < 10%
  if(any(size < min.size)) {
    modes <- modes[-which(size < min.size)]
    mode.dens <- mode.dens[-which(size < min.size)]
    size <- size[-which(size < min.size)]
  }
  if(sum(size) > 1) size <- size / sum(size)
  #Output
  return(list(modes=modes, mode.dens=mode.dens, size=size))
}

is.multimodal <- function(x, min.size=0.01)
{
  if(length(Modes(x, min.size)[[1]]) > 1) return(TRUE)
  else return(FALSE)
}



HPDM <- function(obj, e = 0, prob=0.95, min.size=.01, pointmass = 0, plot=FALSE, print=FALSE){
  if(is.vector(obj)) obj<-as.matrix(obj)
  vals <- apply(obj, 2, sort)
  if(!is.matrix(vals)) stop("obj must have nsamp > 1.")
  nsamp <- nrow(vals)
  npar <- ncol(vals)
  gap <- max(1, min(nsamp - 1, round(nsamp * prob)))
  init <- 1:(nsamp - gap)
  inds <- apply(vals[init + gap, , drop=FALSE] -
                  vals[init, , drop=FALSE], 2, which.min)
  ansmm <- cbind(vals[cbind(inds, 1:npar)],
                 vals[cbind(inds + gap, 1:npar)])
  dimnames(ansmm) <- list(colnames(obj), c("Lower", "Upper"))
  
  mm <- apply(obj, 2, is.multimodal, min.size)
  if(any(mm)) {
    if(print) {cat("\n\nPotentially multimodal column vectors:\n",
        which(mm),"\n")}
    vals <- apply(obj, 2, sort)
    if(!is.matrix(vals)) stop("obj must have nsamp > 1.")
    for (m in which(mm)) {
      X<- vals[,m]
      n<- length(X)
      zeroes<- which(abs(X-pointmass)<=e)
      if(length(zeroes)>=length(X)-1){ #because need 2 points for density
        ansmm[m,]<-c(pointmass,pointmass)
      }
      else{
        if(length(zeroes)==0){ 
          pi<- length(X)/n
          kde<-density(X)
          dens <- rbind(data.frame(approx(kde$x, kde$y, X)),data.frame(x=pointmass,y=0))
          dens$mix <- dens$y
        }
        else{ 
          d<-X[-zeroes]
          epsilons<-c(abs(X[min(max(zeroes)+1,length(X))]),abs(X[max(min(zeroes)-1,1)]))
          epsilon<-ifelse(min(epsilons)==0, max(epsilons),min(epsilons))/20
          pi<- length(d)/n
          kde <- density(d)
          dens <- rbind(data.frame(approx(kde$x, kde$y, d)), 
                        data.frame(x = rep(pointmass,(n-length(d))),y= rep(0,(n-length(d)))))
          #mix of normals
          dens$mix<- pi*dens$y+(1-pi)*dnorm(dens$x, pointmass, epsilon)/dnorm(pointmass, pointmass, epsilon)
        }
        
        dens<-dens[order(dens$x),]
        dens.ind <- dens$mix >= as.vector(quantile(dens$mix,
                                                   probs=1-prob)) * 1
        
        ints <- ""
        count <- 1
        for (i in 1:nrow(vals)) {
          if((i == 1) & (dens.ind[i] == 1)) {
            ints <- paste("(",round(vals[i,m],3),",",sep="")
            if(count > ncol(ansmm)) ansmm <- cbind(ansmm,NA)
            ansmm[m,count] <- vals[i,m]
            count <- count + 1
          }
          if(i > 1) {
            if((dens.ind[i] == 0) & (dens.ind[i-1] == 1)) {
              ints <- paste(ints,round(vals[i-1,m],3),")",sep="")
              if(count > ncol(ansmm)) ansmm <- cbind(ansmm,NA)
              ansmm[m,count] <- vals[i-1,m]
              count <- count + 1
            }
            if((dens.ind[i] == 1) & (dens.ind[i-1] == 0))  {
              ints <- paste(ints," (",round(vals[i,m],3),",",sep="")
              if(count > ncol(ansmm)) ansmm <- cbind(ansmm,NA)
              ansmm[m,count] <- vals[i,m]
              count <- count + 1
            }
          }
        }
        if((dens.ind[i] == 1) & (dens.ind[i-1] == 1)) {
          ints <- paste(ints,round(vals[i,m],3),")",sep="")
          if(count > ncol(ansmm)) ansmm <- cbind(ansmm,NA)
          ansmm[m,count] <- vals[i,m]
          count <- count + 1
        }
        if(print){cat("\nColumn", m, "multimodal intervals:", ints, "\n")}
        if(plot){
          plot(dens$x, dens$mix, type = "l")
          #plotvar(X,e)
          points(ansmm, dens$mix[sapply(ansmm, function(a) which(dens$x==a)[1])], col="red")
        }
      }
      
    }
    return(ansmm)
  }
  else{
    return(ansmm)
  }
}

plotvar = function(x, e = 1e-04, nsteps = 500, pointmass = 0, newplot=TRUE,gg=FALSE) {
  zeroes = which(abs(x-pointmass)<e)
  prob0=length(zeroes)/length(x)
  xne0= x
  if(prob0>0){
    xne0=x[-zeroes]
  }
  if(prob0==1){
    xlower = pointmass-e
    xupper = pointmass+e
    xmax = 1
  }
  m=mean(xne0)
  s= sd(xne0)
  #qmin = min(qnorm(e/2, m, s ))
  #qmax = max(qnorm(1 - e/2, m, s))
  #xlower = min(qmin, 0)
  #xupper = max(0, qmax)
  xlower=min(max(qnorm(e/2, m, s ),min(x)),pointmass)
  xupper=max(min(max(x),qnorm(1 - e/2, m, s)),pointmass)
  
  xx = seq(xlower, xupper, length.out = nsteps)
  yy = rep(0, times = length(xx))
  maxyy = 1
  if (prob0 < 1 ) {
    # kdeneg<- density(xne0[xne0<0])
    # kdepos<-density(xne0[xne0>0])
    # yyneg<- approx(kdeneg$x, kdeneg$y, xx)$y*length(xne0[xne0<0])/length(xne0)
    # yyneg[is.na(yyneg)]<-0
    # 
    # yypos<- approx(kdepos$x, kdepos$y, xx)$y*length(xne0[xne0>0])/length(xne0)
    # yypos[is.na(yypos)]<-0
    # yy = yyneg+yypos
    kde<- density(xne0)
    yy= approx(kde$x, kde$y, xx,yleft=0, yright=0)$y
    
    #yy = dt(x=(x-m)/s, df=)/s
    maxyy = max(yy)/(1-prob0)
  }
  
  ymax = max(prob0, 1 - prob0)
  if(newplot){
    plot(c(xlower, xupper), c(0, ymax), type = "n",
         xlab = "", ylab = "")
  }
  lines(c(pointmass, pointmass), c(0, prob0), lty = 1, lwd = 3,col=as.numeric(newplot)+1)
  lines(xx, (1 - prob0) * yy/maxyy, lty = 1, lwd = 1,col=as.numeric(newplot)+1)
  
  ##ggplot
  if(gg==TRUE){
  return(
                  data.frame(maxy = maxyy, x = xx, y = yy, xl=pointmass, xu=pointmass,yl=0,yu=prob0))
  }
  
  #invisible()
}
