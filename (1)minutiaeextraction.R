db <- "C:\\Users\\Abhi\\Desktop\\research\\fingerprint recognition\\DB 2002\\DB1_B\\"
for (a in 1:100)
{
 for (b in 1:8)
 {
  concat <- paste(c(db,101,'_',1,'.tif'),collapse='')
  image<-readTIFF(concat)
  imagefile <- image*255
  pimage(imagefile)
    
              ######STEP1 : HISTOGRAM EQUALIZATION########
    
  GIm <- imagefile
  numofpixels <- nrow(GIm)*ncol(GIm)
  row <- nrow(GIm)
  col <- ncol(GIm)
  HIm <- array(0,dim=c(row,col))
  freq <- array(0,dim=c(256,1))
  probf <- array(0,dim=c(256,1))
  probc <- array(0,dim=c(256,1))
  cum <- array(0,dim=c(256,1))
  output <- array(0,dim=c(256,1))
  for (i in 1:row)
  {
   for (j in 1:col)
   {
    value <- GIm[i,j]
    freq[value+1] <- freq[value+1]+1
    probf[value+1] <- freq[value+1]/numofpixels
   }
  }
  sum <- 0
  no_bins <-255
  for (i in 1:length(probf))
  {
   sum <- sum+freq[i]
   cum[i] <- sum
   probc[i] <- cum[i]/numofpixels
   output[i] <- round(probc[i]*no_bins)
  }
  for (i in 1:row)
  {
   for (j in 1:col)
   {
    HIm[i,j] <- output[GIm[i,j]+1]
   }
  }
  pimage(HIm)
    
  ####### STEP 2: FUZZIFICATION ########
    
  y <- HIm
  min <- min(y)
  max <- max(y)
  row <- nrow(y)
  col <- ncol(y)
  mug <- array(0,dim=c(row,col))
  for(i in 1:row)
  {
   for(j in 1:col)
   {
    if((y[i,j] >= 0)&(y[i,j] < min))
    {
     mug[i,j] <- 0
    }
    if ((y[i,j] >= min)&(y[i,j] <= max))
    {
     mug[i,j] <- (y[i,j]-min)/(max-min)  
    }
    if((y[i,j] > max)&(y[i,j] <= 255))
    {
     mug[i,j] <- 1
    }
   }
  }
  mu1 <- array(0,dim=c(row,col))
  for (i in 1:row)
  {
   for(j in 1:col)
   {
    if((mug[i,j] >= 0)&(mug[i,j]<=0.5))
    {
      mu1[i,j] <- 2*(mug[i,j]^2)
    }
    if((mug[i,j] > 0.5)&(mug[i,j]<=1.0))
    {
      mu1[i,j] <- 1-(2*((1-mug[i,j])^2))
    }
   }
  }
  mu2 <- array(0,dim=c(row,col))
  for (i in 1:row)
  {
   for(j in 1:col)
   {
    if((mu1[i,j] >= 0)&(mu1[i,j]<=0.5))
    {
     mu2[i,j] <- 2*(mu1[i,j]^2)
    }
    if((mu1[i,j] > 0.5)&(mu1[i,j]<=1.0))
    {
      mu2[i,j] <- 1-(2*((1-mu1[i,j])^2))
    }
   }
  }
  mu3 <- array(0,dim=c(row,col))
  for (i in 1:row)
  {
   for(j in 1:col)
   {
    if((mu2[i,j] >= 0)&(mu2[i,j]<=0.5))
    {
      mu3[i,j] <- 2*(mu2[i,j]^2)
    }
    if((mu2[i,j] > 0.5)&(mu2[i,j]<=1.0))
    {
      mu3[i,j] <- 1-(2*((1-mu2[i,j])^2))
    }
   }
  }
  defuzzy1 <- array(0,dim=c(row,col))
  for(i in 1:row)
  {
   for(j in 1:col)
   {
    if((y[i,j]>=0)&(y[i,j]<min))
    {
     defuzzy1[i,j] <- y[i,j]
    }
    if((y[i,j]>=min)&(y[i,j]<=max))
    {
     defuzzy1[i,j] <- (max-min)*mu3[i,j]
    }
    if((y[i,j]>max)&(y[i,j]<=255))
    {
      defuzzy1[i,j] <- y[i,j]
    }
   }
  }
  pimage(defuzzy1)
    
    ######## STEP3 : Binarization ########
    
 mean <- mean(defuzzy1)
 for (i in 1:row)
 {
  for (j in 1:col)
  {
   if (defuzzy1[i,j]>mean)
   {
    defuzzy1[i,j] <- 1
   }
   else
   {
    defuzzy1[i,j] <- 0
   }
  }
 }
 pimage(defuzzy1)
    
    ####### STEP4 : THINNING #################
    
 thin <- thinImage(defuzzy1)
 absDiff <- function(matrix1,matrix2)
 {
  r <- nrow(matrix1)
  c <- ncol(matrix1)
  destMatrix <- matrix1
  for(r in 0:r-1)
  {
    for(c in 0:c-1)
    {
     destMatrix[r,c] <- abs(matrix1[r,c]-matrix1[r,c])
    }
  }
  return(destMatrix)
 }
 countNonZero <- function(inputMatrix)
 {
  return(length(inputMatrix[inputMatrix > 0]))
 }
 thinningIteration <- function(imageMatrix, iter)
 {
  imageInput <- imageMatrix
  r <- nrow(imageInput) - 1
  c <- ncol(imageInput) - 1
  for(i in 2:r)
  {
    for(j in 2:c)
    {
      p2 <- imageInput[i-1, j]
      p3 <- imageInput[i-1, j+1]
      p4 <- imageInput[i, j+1]
      p5 <- imageInput[i+1, j+1]
      p6 <- imageInput[i+1, j]
      p7 <- imageInput[i+1, j-1]
      p8 <- imageInput[i, j-1]
      p9 <- imageInput[i-1, j-1]
      A  <- (p2 == 0 && p3 == 1) + (p3 == 0 && p4 == 1) + 
            (p4 == 0 && p5 == 1) + (p5 == 0 && p6 == 1) + 
            (p6 == 0 && p7 == 1) + (p7 == 0 && p8 == 1) +
            (p8 == 0 && p9 == 1) + (p9 == 0 && p2 == 1)
      B  <- p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9
      if(iter == 0)
      {
        m1 <- (p2 * p4 * p6)
        m2 <- (p4 * p6 * p8)
      }
      else {
         m1 <- (p2 * p4 * p8)
         m2 <- (p2 * p6 * p8)
       }
      if (A == 1 && (B >= 2 && B <= 6) && m1 == 0 && m2 == 0)
      {
        imageInput[i,j] <- 0
      }
     }
    }
    return(imageInput)
   }
   thinImage <- function(imageMatrix)
   {
    im <- imageMatrix
    prev <- im
    repeat
    {
     im <- thinningIteration(im, 0)
     im <- thinningIteration(im, 1)
     diff <- absDiff(im, prev)
     prev <- im
     if(countNonZero(diff) <= 0)
     {
      break
     }
    } 
    pimage(im)
    return(im)
   }
    
    #########STEP5 : ENHANCED THINNING  ###############
    
    t <- c()
    r <- row-1
    c <- col-1
    for (i in 2:r)
    {
      for (j in 2:c)
      {
        sum <- 0
        if (thin[i,j]==1)
        {
          t[1] <- thin[i-1,j]
          t[2] <- thin[i,j-1]
          t[3] <- thin[i+1,j]
          t[4] <- thin[i,j+1]
        }
        for (k in 1:4)
        {
          sum <- sum+t[k]
        }
        if (sum > 2)
        {
          thin[i,j] <- 0
        }
      }
    }
    pimage(thin)
    
    ####### STEP6 : THOROUGH THINNING #########
    
    q <- c()
    p <- c()
    r <- row-2
    c <- col-2
    for (i in 3:r)
    {
      for (j in 3:c)
      {
        q[1] <- thin[i-2,j]
        q[2] <- thin[i-2,j+1]
        q[3] <- thin[i-2,j+2]
        q[4] <- thin[i-1,j+2]
        q[5] <- thin[i,j+2]
        q[6] <- thin[i+1,j+2]
        q[7] <- thin[i+2,j+2]
        q[8] <- thin[i+2,j+1]
        q[9] <- thin[i+2,j]
        q[10] <- thin[i+2,j-1]
        q[11] <- thin[i+2,j-2]
        q[12] <- thin[i+1,j-2]
        q[13] <- thin[i,j-2]
        q[14] <- thin[i-1,j-2]
        q[15] <- thin[i-2,j-2]
        q[16] <- thin[i-2,j-1]
        p <- thin[i,j]
        p[1] <- thin[i-1,j]
        p[2] <- thin[i-1,j+1]
        p[3] <- thin[i, j+1]
        p[4] <- thin[i+1, j+1]
        p[5] <- thin[i+1, j]
        p[6] <- thin[i+1, j-1]
        p[7] <- thin[i, j-1]
        p[8] <- thin[i-1, j-1]
        if ((q[16]==0)&(q[1]==0)&(q[2]==0))
        {
          thin[i-1,j] <- 0
        }
        if ((q[3]==0)&(q[4]==0)&(q[5]==0))
        {
          thin[i-1,j+1] <- 0
        }
        if ((q[4]==0)&(q[5]==0)&(q[6]==0))
        {
          thin[i,j+1] <- 0
        }
        if ((q[7]==0)&(q[8]==0)&(q[9]==0))
        {
          thin[i+1,j+1] <- 0
        }
        if ((q[8]==0)&(q[9]==0)&(q[10]==0))
        {
          thin[i+1,j] <- 0
        }
        if ((q[11]==0)&(q[12]==0)&(q[13]==0))
        {
          thin[i+1,j-1]<- 0
        }
        if ((q[12]==0)&(q[13]==0)&(q[14]==0))
        {
          thin[i,j-1] <- 0
        }
        if ((q[15]==0)&(q[16]==0)&(q[1]==0))
        {
          thin[i-1,j-1] <- 0
        }
      }
    }
    
    ######## STEP7 : minutiae extration  ################
    
    s <- dim(thin)
    cnmat <- array(5,dim=c(row,col))
    p <- c()
    r <- row-1
    co <- col-1
    for (i in 2:r)
    {
      for (j in 2:co)
      {
        sum <- 0
        p[1] <- thin[i,j+1]
        p[2] <- thin[i-1,j+1]
        p[3] <- thin[i-1,j]
        p[4] <- thin[i-1,j-1]
        p[5] <- thin[i,j-1]
        p[6] <- thin[i+1,j-1]
        p[7] <- thin[i+1,j]
        p[8] <- thin[i+1,j+1]
        p[9] <- p[1]
        for (k in 1:8)
        {
          sig <- abs(p[k]-p[k+1] )
          sum <- sum+sig
        }
        cnmat[i,j] <- 0.5*sum
      }
    }
    newcnmat <- array(0,dim=c(row,col))
    p1 <- c()
    r <- row-4
    co <- col-4
    for (i in 4:r)
    {
      for (j in 4:co)
      {
        p1[1] <- thin[i-3,j-3]
        p1[2] <- thin[i-3,j-2]
        p1[3] <- thin[i-3,j-1]
        p1[4] <- thin[i-3,j]
        p1[5] <- thin[i-3,j+1]
        p1[6] <- thin[i-3,j+2]
        p1[7] <- thin[i-3,j+3]
        p1[8] <- thin[i-2,j+3]
        p1[9] <- thin[i-1,j+3]
        p1[10] <- thin[i,j+3]
        p1[11] <- thin[i+1,j+3]
        p1[12] <- thin[i+2,j+3]
        p1[13] <- thin[i+3,j+3]
        p1[14] <- thin[i+3,j+2]
        p1[15] <- thin[i+3,j+1]
        p1[16] <- thin[i+3,j]
        p1[17] <- thin[i+3,j-1]
        p1[18] <- thin[i+3,j-2]
        p1[19] <- thin[i+3,j-3]
        p1[20] <- thin[i+2,j-3]
        p1[21] <- thin[i+1,j-3]
        p1[22] <- thin[i,j-3]
        p1[23] <- thin[i-1,j-3]
        p1[24] <- thin[i-2,j-3]
        count <- 0
        for (k in 2:5)
        {
          if (((p1[k]==0)&(p1[k+1]==1))|((p1[k]==1)&(p1[k+1]==0)))
          {
            count <- count+1
          }
        }
        for (k in 8:11)
        {
          if (((p1[k]==0)&(p1[k+1]==1))|((p1[k]==1)&(p1[k+1]==0)))
          {
            count <- count+1
          }
        }
        for (k in 14:17)
        {
          if (((p1[k]==0)&(p1[k+1]==1))|((p1[k]==1)&(p1[k+1]==0)))
          {
            count <- count+1
          }
        }
        for (k in 20:23)
        {
          if (((p1[k]==0)&(p1[k+1]==1))|((p1[k]==1)&(p1[k+1]==0)))
          {
            count <- count+1
          }
        }
        sum <- 0
        if (count == 2)
        {
          for (k in 1:24)
          {
            sum <- sum+p1[k]
          }
          avg <- sum/24
          newcnmat[i,j] <- avg 
        }
      }
    }
    newmean <- mean(newcnmat)
    newmeanbi <- 1-newmean
    for (i in 1:row)
    {
      for (j in 1:col)
      {
        if (cnmat[i,j]==1)
        {
          if (newcnmat[i,j]<newmean)
          {
            cnmat[i,j] <- 0
          }
        }
        if (cnmat[i,j]==3)
        {
          if (newcnmat[i,j]>newmeanbi)
          {
            cnmat[i,j] <- 0
          }
        }
       }
      }
    ridgerow <- c()
    ridgecol <- c()
    birow <- c()
    bicol <- c()
    x<-1
    y<-1
    k<-1
    l<-1
    bicount <- 0
    termcnt <- 0
    for (i in 1:row)
    {
      for (j in 1:col)
      {
        if (cnmat[i,j]==1)
        {
          ridgerow[x] <- i
          x <- x+1
          ridgecol[y] <- j
          y <- y+1
          termcnt <- termcnt+1
        }
        if (cnmat[i,j]==3)
        {
          birow[k] <- i
          k <- k+1
          bicol[l] <- j
          l <- l+1
          bicount <- bicount+1
        }
      }
    }
    
    bifurge <- cbind(birow,bicol)
    term <- cbind(ridgerow,ridgecol)
    
    ########STEP8 : FALSE MINUTIAE REMOVAL #################
    
    rowval <- c()
    sumrow <- 0
    for (i in 1:row)
    {
      for (j in 1:col)
      {
        if (cnmat[i,j]==1)
        {
          sumrow <- sumrow+1
        }
      }
      val <- sumrow/225
      rowval <- c(rowval,val)
    }
    ird <- sum(rowval)/225
    if (bicount!=0)
    {
      for (x in 1:bicount)
      {
        for  (i in 1:bicount)
        {
          if ((bifurge[i,1]!=bifurge[x,1]) & (bifurge[i,2]!=bifurge[x,2]))
          {
            dist <- sqrt((bifurge[i,1]-bifurge[x,1])^2+(bifurge[i,2]-bifurge[x,2])^2)
            if (dist < ird)
            {
              bifurge[i,1] <- 0
              bifurge[i,2] <- 0
              bifurge[x,1] <- 0
              bifurge[x,2] <- 0
            }
          }
        }
      }
    }
    if (termcnt!=0)
    {
      for (x in 1:termcnt )
      {
        for  (i in 1:termcnt )
        {
          if ((term[i,1]!=term[x,1]) & (term[i,2]!=term[x,2]))
          {
            dist <- sqrt((term[i,1]-term[x,1])^2+(term[i,2]-term[x,2])^2)
            if (dist < ird)
            {
              term[i,1] <- 0
              term[i,2] <- 0
              term[x,1] <- 0
              term[x,2] <- 0
            }
          }
        }
      }
    }
    if (bicount!=0&termcnt!=0)
    {
      for (x in 1:bicount)
      {
        for (i in 1:termcnt)
        {
          if (((bifurge[x,1]!=term[i,1]) & (bifurge[x,2]!=term[i,2])))
          {
            dist <- sqrt((bifurge[x,1]-term[i,1])^2+(bifurge[x,2]-term[i,2])^2)
            if (dist < ird)
            {
              bifurge[x,1] <- 0
              bifurge[x,2] <- 0
              term[i,1] <- 0
              term[i,2] <- 0
            }  
          }
        }
      }
    }
    
    bicnt <- 0
    if (bicount!=0)
    {
      for (i in 1:bicount)
      {
        if ((bifurge[i,1]!=0)&(bifurge[i,2]!=0))
        {
          bicnt <- bicnt+1
        }
      }
    }
    termcount <- 0
    if (termcnt!=0)
    {
      for (i in 1:termcnt)
      {
        if ((term[i,1]!=0)&(term[i,2]!=0))
        {
          termcount <- termcount+1
        }
      }
    }
    minutiae <- rbind(bifurge,term)
    type1 <- c()
    
    if (bicnt!=0)
    {
      for (i in 1:bicnt)
      {
        type1[i] <- 3
      }
    }
    type2 <- c()
    if (termcount!=0)
    {
      for (i in 1:termcount)
      {
        type2[i] <- 1
      }
    }
    type <- c(type1,type2)
    minutiae <- cbind(minutiae,type)
    path <- "C:\\Users\\Abhi\\Documents\\R\\minutiaepoints\\"
    concat <- paste(c(path,a,',',b,'.txt'),collapse='')
    write.table(minutiae, concat, sep="\t")
  }
}
