score <- array(0,dim=c(100,8))
l <- 0
rej <- 0
for (c in 1:100)
{
 l <- l+1
 for (d in 1:8)
 {
  if (d == 1)
  {
   path <- "C:\\Users\\Abhi\\Documents\\R\\minpointsdba\\"
   concat <- paste(c(path,c,',',d,'.txt'),collapse='')
   min1 <- read.table(concat,  sep="\t", col.names=c("birow","bicol","type"), fill=FALSE, strip.white=TRUE)
   minutiae1x <- min1[,1]
   minutiae1y <- min1[,2]
   minutiae1t <- min1[,3]
  }
  if (d > 1)
  {
    path <- "C:\\Users\\Abhi\\Documents\\R\\minpointsdba\\"
    concat <- paste(c(path,c,',',d,'.txt'),collapse='')
    min2 <- read.table(concat,  sep="\t", col.names=c("birow","bicol","type"), fill=FALSE, strip.white=TRUE)
    minutiae2x <- min2[,1]
    minutiae2y <- min2[,2]
    minutiae2t <- min2[,3]
    minnew1 <- rbind(minutiae1x,minutiae1y,minutiae1t) 
    minnew2 <- rbind(minutiae2x,minutiae2y,minutiae2t) 
    trigmat <- array(dim=c(3,1))
    xnew <- c()
    ynew <- c()
    thetanew <- c()
    nmin1 <- nrow(min1)
    nmin2 <- nrow(min2)
    if (nmin1 < nmin2)
    {
      cnt <-nmin1
    } else {
      cnt <- nmin2
    }
    if (cnt >20)
    {
      for (i in 1:cnt)
      {
        trigmat <- matrix(c(cos(min1[1,3]),-sin(min1[1,3]),0,sin(min1[1,3]),cos(min1[1,3]),0,0,0,1),nrow=3,ncol=3) %*% matrix(c((min1[1,1]-min1[i,1]),(min1[1,2]-min1[i,2]),(min1[1,3]-min1[i,3])))
        xnew[i] <- trigmat[1,1]
        ynew[i] <- trigmat[2,1]
        thetanew[i] <- trigmat[3,1]
      }
      minutiae1new <- cbind(xnew,ynew,thetanew)
      trigmat <- array(dim=c(3,1))
      xnew <- c()
      ynew <- c()
      thetanew <- c()
      for (i in 1:cnt)
      {
       trigmat <- matrix(c(cos(min2[1,3]),-sin(min2[1,3]),0,sin(min2[1,3]),cos(min2[1,3]),0,0,0,1),nrow=3,ncol=3) %*% matrix(c((min2[1,1]-min2[i,1]),(min2[1,2]-min2[i,2]),(min2[1,3]-min2[i,3])))
       xnew[i] <- trigmat[1,1]
       ynew[i] <- trigmat[2,1]
       thetanew[i] <- trigmat[3,1]
      }
      minutiae2new <- cbind(xnew,ynew,thetanew)
      sd <- c()
      dd <- c()
      for (i in 1:cnt)
      {
        sd[i] <- sqrt((minutiae1new[i,1]-minutiae2new[i,1])^2+(minutiae1new[i,2]-minutiae2new[i,2])^2)
        dd[i] <- min(abs(minutiae2new[i,3]-minutiae1new[i,3]),360-abs(minutiae2new[i,3]-minutiae1new[i,3]))
      } 
      mm <- 0
      for (i in 1:cnt)
      {
       if (sd[i] < 15)
       {
        mm <- mm+1
       }
       if (dd[i] < 14)
       {
        mm <- mm+1
       }
      }
      score[l,d] <- (mm/cnt)-0.5
      
     }
    }else {
      rej <- rej+1
      score[l,d] <- 0
    }
   }
  }
 path <- "C:\\Users\\Abhi\\Documents\\R\\score\\scoredb_a.txt"
 write.table(score, path, sep="\t")
 scorenew <- array(0,dim=c(100,10000))
 l <- 1
 d <- 1
 nira <- 0
 for (a in 1:100)
 {
  for (b in 1:100)
  {
    path <- "C:\\Users\\Abhi\\Documents\\R\\minpointsdba\\"
    concat <- paste(c(path,a,',',1,'.txt'),collapse='')
    min1 <- read.table(concat,  sep="\t", col.names=c("birow","bicol","type"), fill=FALSE, strip.white=TRUE)
    minutiae1x <- min1[,1]
    minutiae1y <- min1[,2]
    minutiae1t <- min1[,3]
    if (a!=b)
    {
      path <- "C:\\Users\\Abhi\\Documents\\R\\minpointsdba\\"
      c <- a+1
      concat <- paste(c(path,b,',',1,'.txt'),collapse='')
      min2 <- read.table(concat,  sep="\t", col.names=c("birow","bicol","type"), fill=FALSE, strip.white=TRUE)
      minutiae2x <- min2[,1]
      minutiae2y <- min2[,2]
      minutiae2t <- min2[,3]
      minnew1 <- rbind(minutiae1x,minutiae1y,minutiae1t) 
      minnew2 <- rbind(minutiae2x,minutiae2y,minutiae2t) 
      trigmat <- array(dim=c(3,1))
      xnew <- c()
      ynew <- c()
      thetanew <- c()
      nmin1 <- nrow(min1)
      nmin2 <- nrow(min2)
      if (nmin1 < nmin2)
      {
        cnt <-nmin1
      } else {
        cnt <- nmin2
      }
      if (cnt > 25)
      {
        for (i in 1:cnt)
        {
          trigmat <- matrix(c(cos(min1[1,3]),-sin(min1[1,3]),0,sin(min1[1,3]),cos(min1[1,3]),0,0,0,1),nrow=3,ncol=3) %*% matrix(c((min1[1,1]-min1[i,1]),(min1[1,2]-min1[i,2]),(min1[1,3]-min1[i,3])))
          xnew[i] <- trigmat[1,1]
          ynew[i] <- trigmat[2,1]
          thetanew[i] <- trigmat[3,1]
        }
        minutiae1new <- cbind(xnew,ynew,thetanew)
        trigmat <- array(dim=c(3,1))
        xnew <- c()
        ynew <- c()
        thetanew <- c()
        for (i in 1:cnt)
        {
          trigmat <- matrix(c(cos(min2[1,3]),-sin(min2[1,3]),0,sin(min2[1,3]),cos(min2[1,3]),0,0,0,1),nrow=3,ncol=3) %*% matrix(c((min2[1,1]-min2[i,1]),(min2[1,2]-min2[i,2]),(min2[1,3]-min2[i,3])))
          xnew[i] <- trigmat[1,1]
          ynew[i] <- trigmat[2,1]
          thetanew[i] <- trigmat[3,1]
        }
        minutiae2new <- cbind(xnew,ynew,thetanew)
        sd <- c()
        dd <- c()
        for (i in 1:cnt)
        {
          sd[i] <- sqrt((minutiae1new[i,1]-minutiae2new[i,1])^2+(minutiae1new[i,2]-minutiae2new[i,2])^2)
          dd[i] <- min(abs(minutiae2new[i,3]-minutiae1new[i,3]),360-abs(minutiae2new[i,3]-minutiae1new[i,3]))
        } 
        mm <- 0
        for (i in 1:cnt)
        {
          if (sd[i] < 15)
          {
            mm <- mm+1
          }
          if (dd[i] < 14)
          {
            mm <- mm+1
          }
        }
        scorenew[a,b] <- (mm/cnt)+15                                                                                                     
      } else {
        scorenew[l,d] <- 0
        nira <- nira+1
        
      }
    }
  }
}
