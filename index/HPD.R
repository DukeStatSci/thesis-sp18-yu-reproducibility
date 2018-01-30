X<- sapply(runif(1000), function(x) ifelse(x<.3, rnorm(1,0,0), rnorm(1,3,1)))
hist(X)
epsilon <- .001
HPD <- function(X, prob=.95){
  ##remove 0- can probably find/approx with some check
  n<- length(X)
  d<-X[X!=0]
  pi<- length(d)/n
  
  kde <- density(d)
  dens <- rbind(data.frame(approx(kde$x, kde$y, d)), 
                data.frame(x = rep(0,(n-length(d))),y= rep(0,(n-length(d)))))
  #mix of normals
  dens$mix<- pi*dens$y+(1-pi)*dnorm(dens$y, 0, epsilon)
  min <- quantile(dens$mix, c(1-prob))
  dens<-dens[order(dens$x),]
  
  interval <- which(dens$mix>=min)
  plot(dens$x[-which(dens$x==0)], dens$mix[-which(dens$x==0)], type = "l")
  
  plot(dens$x, dens$mix, type = "l")
  points(dens$x[interval], dens$mix[interval], col="red")
  
  setdiff(1:n, interval)
  c(dens$x[3],dens$x[527],dens$x[553],dens$x[977])
  
  
  #https://cran.r-project.org/web/packages/LaplacesDemon/LaplacesDemon.pdf
}


HPDM <- function(obj, prob=0.95){          
  mm <- apply(obj, 2, is.multimodal)
  if(any(mm)) {
    cat("\n\nPotentially multimodal column vectors:\n",
        which(mm),"\n")
    vals <- apply(obj, 2, sort)
    if(!is.matrix(vals)) stop("obj must have nsamp > 1.")
    for (m in which(mm)) {
      kde <- density(vals[,m])
      dens <- approx(kde$x, kde$y, vals[,m])$y
      dens.ind <- dens >= as.vector(quantile(dens,
                                             probs=1-prob)) * 1
      ints <- ""
      count <- 1
      for (i in 1:nrow(vals)) {
        if((i == 1) & (dens.ind[i] == 1)) {
          ints <- paste("(",round(vals[i,m],3),",",sep="")
          #if(count > ncol(ansmm)) ansmm <- cbind(ansmm,NA)
          #ansmm[m,count] <- vals[i,m]
          count <- count + 1
        }
        if(i > 1) {
          if((dens.ind[i] == 0) & (dens.ind[i-1] == 1)) {
            ints <- paste(ints,round(vals[i-1,m],3),")",sep="")
            #if(count > ncol(ansmm)) ansmm <- cbind(ansmm,NA)
            #ansmm[m,count] <- vals[i-1,m]
            count <- count + 1
          }
          if((dens.ind[i] == 1) & (dens.ind[i-1] == 0))  {
            ints <- paste(ints," (",round(vals[i,m],3),",",sep="")
            #if(count > ncol(ansmm)) ansmm <- cbind(ansmm,NA)
            #ansmm[m,count] <- vals[i,m]
            count <- count + 1
          }
        }
      }
      if((dens.ind[i] == 1) & (dens.ind[i-1] == 1)) {
        ints <- paste(ints,round(vals[i,m],3),")",sep="")
        #if(count > ncol(ansmm)) ansmm <- cbind(ansmm,NA)
        #∫ansmm[m,count] <- vals[i,m]
        count <- count + 1
      }
      cat("\nColumn", m, "multimodal intervals:", ints, "\n")
    }
  }
}
####not sure if this is necessary
Modes <- function(x, min.size=0.1) {
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



HPDM2 <- function(obj, prob=0.95){          
  mm <- apply(obj, 2, is.multimodal)
  if(any(mm)) {
    cat("\n\nPotentially multimodal column vectors:\n",
        which(mm),"\n")
    vals <- apply(obj, 2, sort)
    if(!is.matrix(vals)) stop("obj must have nsamp > 1.")
    for (m in which(mm)) {
      X<- vals[,m]
      n<- length(X)
      zeroes<- which(X==0)
      d<-X[-zeroes]
      pi<- length(d)/n
      kde <- density(d)
      dens <- rbind(data.frame(approx(kde$x, kde$y, d)), 
                    data.frame(x = rep(0,(n-length(d))),y= rep(0,(n-length(d)))))
      dens<-dens[order(dens$x),]
      #mix of normals
      epsilon<- min(abs(X[max(zeroes)+1]),abs(X[min(zeroes)-1]))/20
      dens$mix<- pi*dens$y+(1-pi)*dnorm(dens$x, 0, epsilon)
      dens.ind <- dens$mix >= as.vector(quantile(dens$mix,
                                             probs=1-prob)) * 1
      
      ints <- ""
      count <- 1
      for (i in 1:nrow(vals)) {
        if((i == 1) & (dens.ind[i] == 1)) {
          ints <- paste("(",round(vals[i,m],3),",",sep="")
          #if(count > ncol(ansmm)) ansmm <- cbind(ansmm,NA)
          #ansmm[m,count] <- vals[i,m]
          count <- count + 1
        }
        if(i > 1) {
          if((dens.ind[i] == 0) & (dens.ind[i-1] == 1)) {
            ints <- paste(ints,round(vals[i-1,m],3),")",sep="")
            #if(count > ncol(ansmm)) ansmm <- cbind(ansmm,NA)
            #ansmm[m,count] <- vals[i-1,m]
            count <- count + 1
          }
          if((dens.ind[i] == 1) & (dens.ind[i-1] == 0))  {
            ints <- paste(ints," (",round(vals[i,m],3),",",sep="")
            #if(count > ncol(ansmm)) ansmm <- cbind(ansmm,NA)
            #ansmm[m,count] <- vals[i,m]
            count <- count + 1
          }
        }
      }
      if((dens.ind[i] == 1) & (dens.ind[i-1] == 1)) {
        ints <- paste(ints,round(vals[i,m],3),")",sep="")
        #if(count > ncol(ansmm)) ansmm <- cbind(ansmm,NA)
        #∫ansmm[m,count] <- vals[i,m]
        count <- count + 1
      }
      cat("\nColumn", m, "multimodal intervals:", ints, "\n")
      plot(dens$x, dens$mix, type = "l")
      points(dens$x[dens.ind==1], dens$mix[interval], col="red")

    }
  }
}

##testing
X<- sapply(runif(1000), function(x) ifelse(x<.3, rnorm(1,0,0), rnorm(1,3,1)))
HPDM2(matrix(X))

plot(density(X))

X<- sapply(runif(1000), function(x) ifelse(x<.5, rnorm(1,0,0), rnorm(1,3,1.5)))
HPDM2(matrix(X))


X<- sapply(runif(1000), function(x) ifelse(x<.5, rnorm(1,0,0), rnorm(1,2,1.5)))
HPDM2(matrix(X))

X<- sapply(runif(1000), function(x) ifelse(x<.8, rnorm(1,0,0), rnorm(1,2,1.5)))
HPDM2(matrix(X))


X<- sapply(runif(1000), function(x) ifelse(x<.05, rnorm(1,0,0), rnorm(1,2,.5)))
HPDM2(matrix(X))

#what to do about this? would expect 0 not to be included bc <0.01 
X<- sapply(runif(1000), function(x) ifelse(x<.01, rnorm(1,0,0), rnorm(1,2,.5)))
HPDM2(matrix(X))

