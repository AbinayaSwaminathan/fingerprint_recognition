fnmr <- c()
k <- 1
ngra <- 800-rej
for (i in 1:100)
{
  for (j in 1:8)
  {
    if (score[i,j]<threshold)
    {
      fnmr[k] <- ((score[i,j]+rej)/ngra)
      k <- k+1
    }
  }
}
fnmrfrr <-mean(fnmr)
score1 <- scorenew
for (i in 1:100)
{
  for (j in 1:100)
  {
    if (scorenew[i,j]!=0)
    {
      scorenew[i,j] <- score1[i,j]+15
    }
  }
}
k <- 1
fmr <- c()
for (i in 1:100)
{
  for (j in 1:100)
  {
    if (scorenew[i,j]>=threshold)
    {
      fmr[k] <- scorenew[i,j]/rej
      k <- k+1
    }
  }
}
fmrfar <- mean(fmr)
