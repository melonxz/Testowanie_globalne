set.seed(123)
library("car")
library(goftest)
library(knitr)
library(ddpca)
 
# Sprawdzenie rozkładu p-wartości przy prawdziwym H_0
  options(scipen = 0.1)
  par(mfrow = c(2,3))
  n <- c(100,1000)
  l <- 500
  lambda0 <- c(0.5,5,50)
  alpha  <- 0.05
  for(i in n)
  {    
  for(l0 in lambda0)
  {

        x <- c()
        X <- rpois(i,l*lambda0)
        x <- 1- ppois(X,l*lambda0)
        qqPlot(x,distribution = "unif", col.lines = "red", main=sprintf("n= %d, lambda = %g", i,l0), pch = 16, ylab = "kwantyle p-wartości",xlab="kwantyle teoretyczne",id = FALSE)
    }
  }

# Błąd I rodzaju
  
  l <- 500 #ustalona długość wektorów sumujących się do statystyk testowych
  lambda0 <- seq(0.1,20,0.1)
  alpha  <- 0.05 #poziom istotnosci
  n <- 1000 #liczba hipotez
  sbonf <- c()
  sf <- c()
  sks <- c()
  sad <- c()
  scvm <-c()
  shc <- c()
  
f <- function(n) 
{
  for(l1 in lambda0)
  {
    bonf <- c()
    KS <- c()
    F <- c()
    AD <- c()
    cvm<-c()
    hc <- c()
    for(razy in 1:1000)
    {
      x <- c()
      X <- rpois(n,l*l1) #wektor statystyk testowych
      x <- 1- ppois(X,l*l1) #p-wartosci
      c <- cvm.test(x, 'punif')
      cvm <- c(cvm, c$p.value < alpha )
      hc <- c(hc, HCdetection(x,alpha)$H)
      bonf <- c(bonf,min(x) <= (alpha/n))
      F <- c(F,-sum(2*log(x))>qchisq(1-alpha,2*n))
      a = ks.test(x,'punif')
      KS <- c(KS,a$p.value < alpha)
      b = ad.test(x, null = 'punif') 
      AD <- c(AD,b$p.value < alpha)
    }
    sbonf <- c(sbonf,sum(bonf))
    sf <- c(sf,sum(F))
    sks <- c(sks,sum(KS))
    sad <- c(sad,sum(AD))
    shc <- c(shc, sum(hc))
    scvm <- c(scvm, sum(cvm))
  }
  x = data.frame( "Lambda" = lambda0,"Bonferroniego" = sbonf/1000, "Fishera" =  sf/1000, "KS" = sks/1000, "CvM" = scvm/1000, "AD" = sad/1000, "HC" = shc/1000)
  return(x)
}

x1 <- f(1000)

c1 <- c(1,5,10,15,20,30,40,50,75)
kable(x1[c1,],"latex") # tabela wybranych wartości dla n= 1000 w zal. od lambda_0
par(mfrow = c(1,1))
# błąd I rodzaju dla n = 1000 w zal. od lambda_0
plot(x1[,2]~lambda0, type = "l",  col = "red", ylim= c(0,1), ylab = "Frakcja fałszywych odkryć", xlab = "Lambda")
abline(h=0.05, lty =3)
lines(x1[,3]~lambda0, type = "l", col = "orange", ylab = "Frakcja fałszywych odkryć", xlab = "Lambda", main = "F")
lines(x1[,4]~lambda0, type = "l", col = "green2", ylab = "Frakcja fałszywych odkryć", xlab = "Lambda", main = "KS")
lines(x1[,5]~lambda0, type = "l", col = "purple", ylab = "Frakcja fałszywych odkryć", xlab = "Lambda", main = "CvM" )
lines(x1[,6]~lambda0, type = "l", col = "cyan2", ylab = "Frakcja fałszywych odkryć", xlab = "Lambda", main = "AD")
lines(x1[,7]~lambda0, type = "l", col = "slategray3", ylab = "Frakcja fałszywych odkryć", xlab = "Lambda", main = "HC" )
legend("topright", c("Bonf.","F","KS","CvM","AD","HC"), col = c("red","orange","green2","purple","cyan2","slategray3"), lwd = 2)

#błąd I rodzaju w zależności od n

lambda0 <-  5
n <- c(5,10,20,50,100,250,500,1000,1500,2500,4000,6000,8000,10000)
l <- 500
f2 <- function() 
{
  for(n1 in n)
  {
    bonf <- c()
    KS <- c()
    F <- c()
    AD <- c()
    cvm<-c()
    hc <- c()
    for(razy in 1:1000)
    {
      x <- c()
      X <- rpois(n1,l*lambda0)
      x <- 1- ppois(X,l*lambda0)
      c <- cvm.test(x, 'punif')
      cvm <- c(cvm, c$p.value < alpha )
      hc <- c(hc, HCdetection(x,alpha)$H)
      bonf <- c(bonf,min(x) <= (alpha/n1))
      F <- c(F,-sum(2*log(x))>qchisq(1-alpha,2*n1))
      a = ks.test(x,'punif')
      KS <- c(KS,a$p.value < alpha)
      b = ad.test(x, null = 'punif') 
      AD <- c(AD,b$p.value < alpha)
    }
    sbonf <- c(sbonf,sum(bonf))
    sf <- c(sf,sum(F))
    sks <- c(sks,sum(KS))
    sad <- c(sad,sum(AD))
    shc <- c(shc, sum(hc))
    scvm <- c(scvm, sum(cvm))
  }
  x = data.frame( "Liczba hipotez" = n,"Bonferroniego" = sbonf/1000, "Fishera" =  sf/1000, "KS" = sks/1000, "CvM" = scvm/1000, "AD" = sad/1000, "HC" = shc/1000)
  return(x)
}

x3 <- f2()


c1 <- c(1,3,5,7,8,11,14)
kable(x3[c1,],"latex") # tabela wybranych wartości dla lambda_0 = 5


par(mfrow = c(1,1))

# błąd I rodzaju dla lambda_0 = 5 w zal. od n
plot(x3[,2]~n, type = "l", col = "red", ylim = c(0,1), ylab = "Frakcja fałszywych odkryć", xlab = "n")
abline(h=0.05, lty =3)
lines(x3[,3]~n, type = "l", col = "orange", ylab = "Frakcja fałszywych odkryć")
lines(x3[,4]~n, type = "l", col = "green2", ylab = "Frakcja fałszywych odkryć")
lines(x3[,5]~n, type = "l", col = "purple", ylab = "Frakcja fałszywych odkryć")
lines(x3[,6]~n, type = "l", col = "cyan2", ylab = "Frakcja fałszywych odkryć")
lines(x3[,7]~n, type = "l", col = "slategray3", ylab = "Frakcja fałszywych odkryć")
legend(x = 7850, y = 0.6, c("Bonf.","F","KS","CvM","AD","HC"), col = c("red","orange","green2","purple","cyan2","slategray3"), lwd = 2)

# Igła w stogu siana


l <- 500
lambda0 <- seq(0.1,4,0.1)
l0 <- 5
alpha  <- 0.05
n<- 1000
sbonf <- c()
sf <- c()
sks <- c()
sad <- c()
scvm <-c()
shc <- c()

f3 <- function() 
{
  for(l1 in lambda0)
  {
    bonf <- c()
    KS <- c()
    F <- c()
    AD <- c()
    cvm<-c()
    hc <- c()
    for(razy in 1:1000)
    {
      k <- sample(seq(1,n),1)
      X <- rpois(n,l*l0)
      X[k] <- rpois(1,l*(l0+l1))
      x <- 1- ppois(X,l*l0)
      c <- cvm.test(x, 'punif')
      cvm <- c(cvm, c$p.value < alpha )
      hc <- c(hc, HCdetection(x,alpha)$H)
      bonf <- c(bonf,min(x) <= (alpha/n))
      F <- c(F,-sum(2*log(x))>qchisq(1-alpha,2*n))
      a = ks.test(x,'punif')
      KS <- c(KS,a$p.value < alpha)
      b = ad.test(x, null = 'punif') 
      AD <- c(AD,b$p.value < alpha)
    }
    sbonf <- c(sbonf,sum(bonf))
    sf <- c(sf,sum(F))
    sks <- c(sks,sum(KS))
    sad <- c(sad,sum(AD))
    shc <- c(shc, sum(hc))
    scvm <- c(scvm, sum(cvm))
  }
  x = data.frame( "Lambda" = lambda0,"Bonferroniego" = sbonf/1000, "Fishera" =  sf/1000, "KS" = sks/1000, "CvM" = scvm/1000, "AD" = sad/1000, "HC" = shc/1000)
  return(x)
}

x4 <- f3()

c1 <- c(1,5,10,20,30,40)
kable(x4[c1,],"latex") # tabela wybranych wartości dla lambda_0 = 5 n= 1000


par(mfrow = c(1,1))

# moc dla lambda_0 = 5 w zal. od lambda_1 - lambda_0
plot(x4[,2]~lambda0, type = "l", col = "red", ylim = c(0,1), ylab = "Frakcja odkryć", xlab = "Wartość sygnału")
abline(h=0.05, lty =3)
lines(x4[,3]~lambda0, type = "l", col = "orange")
lines(x4[,4]~lambda0, type = "l", col = "green2")
lines(x4[,5]~lambda0, type = "l", col = "purple")
lines(x4[,6]~lambda0, type = "l", col = "cyan2")
lines(x4[,7]~lambda0, type = "l", col = "slategray3")
legend("topleft", c("Bonf.","F","KS","CvM","AD","HC"), col = c("red","orange","green2","purple","cyan2","slategray3"), lwd = 2, cex = 0.9)


# igła w stogu siana ze zmiennym lambda_0 

l <- 500
lambda0 <-  seq(0.2,10,0.2)
alpha  <- 0.05
n<- 1000
sbonf <- c()
sf <- c()
sks <- c()
sad <- c()
scvm <-c()
shc <- c()

f5 <- function() 
{
  for(l0 in lambda0)
  {
    bonf <- c()
    KS <- c()
    F <- c()
    AD <- c()
    cvm<-c()
    hc <- c()
    for(razy in 1:1000)
    {
      k <- sample(seq(1,n),1)
      X <- rpois(n,l*l0)
      X[k] <- rpois(1,l*(l0+0.5))
      x <- 1- ppois(X,l*l0)
      c <- cvm.test(x, 'punif')
      cvm <- c(cvm, c$p.value < alpha )
      hc <- c(hc, HCdetection(x,alpha)$H)
      bonf <- c(bonf,min(x) <= (alpha/n))
      F <- c(F,-sum(2*log(x))>qchisq(1-alpha,2*n))
      a = ks.test(x,'punif')
      KS <- c(KS,a$p.value < alpha)
      b = ad.test(x, null = 'punif') 
      AD <- c(AD,b$p.value < alpha)
    }
    sbonf <- c(sbonf,sum(bonf))
    sf <- c(sf,sum(F))
    sks <- c(sks,sum(KS))
    sad <- c(sad,sum(AD))
    shc <- c(shc, sum(hc))
    scvm <- c(scvm, sum(cvm))
  }
  x = data.frame( "Lambda" = lambda0,"Bonferroniego" = sbonf/1000, "Fishera" =  sf/1000, "KS" = sks/1000, "CvM" = scvm/1000, "AD" = sad/1000, "HC" = shc/1000)
  return(x)
}

x5 <- f5()

c1 <- c(1,4,7,11,14,20,27,34,50)
kable(x5[c1,],"latex") # tabela wybranych wartości dla  n= 1000 ze zmiennym lambda_0


par(mfrow = c(1,1))

# moc w zal. od lambda_0 dla lambda_1 - lambda_0 = 2
plot(x5[,2]~lambda0, type = "l", col = "red", ylim = c(0,1), ylab = "Frakcja odkryć", xlab = "Lambda_0")
abline(h=0.05, lty =3)
lines(x5[,3]~lambda0, type = "l", col = "orange")
lines(x5[,4]~lambda0, type = "l", col = "green2")
lines(x5[,5]~lambda0, type = "l", col = "purple",)
lines(x5[,6]~lambda0, type = "l", col = "cyan2")
lines(x5[,7]~lambda0, type = "l", col = "slategray3")
legend("topright", c("Bonf.","F","KS","CvM","AD","HC"), col = c("red","orange","green2","purple","cyan2","slategray3"), lwd = 2)

# igła w stogu siana ze zmeinnym n

l <- 500
l0 <- 5
alpha  <- 0.05
n <- c(5,10,20,50,100,250,500,1000,1500,2500,4000,7000,10000)
sbonf <- c()
sf <- c()
sks <- c()
sad <- c()
scvm <-c()
shc <- c()

f6 <- function() 
{
  for(n1 in n)
  {
    bonf <- c()
    KS <- c()
    F <- c()
    AD <- c()
    cvm<-c()
    hc <- c()
    for(razy in 1:1000)
    {
      k <- sample(seq(1,n1),1)
      X <- rpois(n1,l*l0)
      X[k] <- rpois(1,l*(l0+0.5))
      x <- 1- ppois(X,l*l0)
      c <- cvm.test(x, 'punif')
      cvm <- c(cvm, c$p.value < alpha )
      hc <- c(hc, HCdetection(x,alpha)$H)
      bonf <- c(bonf,min(x) <= (alpha/n1))
      F <- c(F,-sum(2*log(x))>qchisq(1-alpha,2*n1))
      a = ks.test(x,'punif')
      KS <- c(KS,a$p.value < alpha)
      b = ad.test(x, null = 'punif') 
      AD <- c(AD,b$p.value < alpha)
    }
    sbonf <- c(sbonf,sum(bonf))
    sf <- c(sf,sum(F))
    sks <- c(sks,sum(KS))
    sad <- c(sad,sum(AD))
    shc <- c(shc, sum(hc))
    scvm <- c(scvm, sum(cvm))
  }
  x = data.frame( "n" = n,"Bonferroniego" = sbonf/1000, "Fishera" =  sf/1000, "KS" = sks/1000, "CvM" = scvm/1000, "AD" = sad/1000, "HC" = shc/1000)
  return(x)
}

x6 <- f6()

c1 <- c(1,3,5,7,8,11,13)
kable(x6[c1,],"latex") # tabela wybranych wartości dla problemu igły w stogu siana ze zmiennym n


par(mfrow = c(1,1))

# moc w zal. n
plot(x6[,2]~n, type = "l", col = "red", ylim = c(0,1), ylab = "Frakcja odkryć", xlab = "n")
abline(h=0.05, lty =3)
lines(x6[,3]~n, type = "l", col = "orange")
lines(x6[,4]~n, type = "l", col = "green2")
lines(x6[,5]~n, type = "l", col = "purple")
lines(x6[,6]~n, type = "l", col = "cyan2")
lines(x6[,7]~n, type = "l", col = "slategray3")
legend(x = 2000, y= 0.65, c("Bonf.","F","KS","CvM","AD","HC"), col = c("red","orange","green2","purple","cyan2","slategray3"), lwd = 2)


#Many small effects ze zmiennym epsilonem

l <- 500
epsilon <- seq(1,30,1)
l0 <- 5
alpha  <- 0.05
n<- 1000
sbonf <- c()
sf <- c()
sks <- c()
sad <- c()
scvm <-c()
shc <- c()

f8 <- function(s) 
{
  for(e in epsilon)
  {
    n1 <- floor(s*n)
    e1 <- e/(s*n)
    bonf <- c()
    KS <- c()
    F <- c()
    AD <- c()
    cvm<-c()
    hc <- c()
    for(razy in 1:1000)
    {
      k <- sample(seq(1,n),n1)
      X <- rpois(n,l*l0)
      X[k] <- rpois(n1,l*(l0+e1))
      x <- 1- ppois(X,l*l0)
      c <- cvm.test(x, 'punif')
      cvm <- c(cvm, c$p.value < alpha )
      hc <- c(hc, HCdetection(x,alpha)$H)
      bonf <- c(bonf,min(x) <= (alpha/n))
      F <- c(F,-sum(2*log(x))>qchisq(1-alpha,2*n))
      a = ks.test(x,'punif')
      KS <- c(KS,a$p.value < alpha)
      b = ad.test(x, null = 'punif') 
      AD <- c(AD,b$p.value < alpha)
    }
    sbonf <- c(sbonf,sum(bonf))
    sf <- c(sf,sum(F))
    sks <- c(sks,sum(KS))
    sad <- c(sad,sum(AD))
    shc <- c(shc, sum(hc))
    scvm <- c(scvm, sum(cvm))
  }
  x = data.frame( "epsilon" = epsilon,"Bonferroniego" = sbonf/1000, "Fishera" =  sf/1000, "KS" = sks/1000, "CvM" = scvm/1000, "AD" = sad/1000, "HC" = shc/1000)
  return(x)
}


x81 <- f8(0.5)


k <- c(1,5,10,15,20,25,30)

kable(x81[k,],"latex") # tabela wybranych wartości dla lambda_0 = 5 n= 1000


par(mfrow = c(1,1))

# moc dla lambda_0 = 5 w zal. od lambda_1 - lambda_0
plot(x81[,2]~epsilon, type = "l", col = "red", ylim = c(0,1), ylab = "Frakcja odkryć")
lines(x81[,3]~epsilon, type = "l", col = "orange")
lines(x81[,4]~epsilon, type = "l", col = "green2")
lines(x81[,5]~epsilon, type = "l", col = "purple")
lines(x81[,6]~epsilon, type = "l", col = "cyan2")
lines(x81[,7]~epsilon, type = "l", col = "slategray3")
abline(h=0.05, lty =3)
legend(x = 23.4, y= 0.8, c("Bonf.","F","KS","CvM","AD","HC"), col = c("red","orange","green2","purple","cyan2","slategray3"), lwd = 2)



#Many small effects ze zmiennym lambda_0

l <- 500
e <- 15
lambda0 <-  seq(0.2,15,0.2)
alpha  <- 0.05
n<- 1000
sbonf <- c()
sf <- c()
sks <- c()
sad <- c()
scvm <-c()
shc <- c()

f9 <- function(s) 
{
  for(l0 in lambda0)
  {
    n1 <- floor(s*n)
    e1 <- e/(s*n)
    bonf <- c()
    KS <- c()
    F <- c()
    AD <- c()
    cvm<-c()
    hc <- c()
    for(razy in 1:1000)
    {
      k <- sample(seq(1,n),n1)
      X <- rpois(n,l*l0)
      X[k] <- rpois(n1,l*(l0+e1))
      x <- 1- ppois(X,l*l0)
      c <- cvm.test(x, 'punif')
      cvm <- c(cvm, c$p.value < alpha )
      hc <- c(hc, HCdetection(x,alpha)$H)
      bonf <- c(bonf,min(x) <= (alpha/n))
      F <- c(F,-sum(2*log(x))>qchisq(1-alpha,2*n))
      a = ks.test(x,'punif')
      KS <- c(KS,a$p.value < alpha)
      b = ad.test(x, null = 'punif') 
      AD <- c(AD,b$p.value < alpha)
    }
    sbonf <- c(sbonf,sum(bonf))
    sf <- c(sf,sum(F))
    sks <- c(sks,sum(KS))
    sad <- c(sad,sum(AD))
    shc <- c(shc, sum(hc))
    scvm <- c(scvm, sum(cvm))
  }
  x = data.frame( "lambda_0" = lambda0,"Bonferroniego" = sbonf/1000, "Fishera" =  sf/1000, "KS" = sks/1000, "CvM" = scvm/1000, "AD" = sad/1000, "HC" = shc/1000)
  return(x)
}

x9 <- f9(0.5)

k <- c(1,5,10,20,35,50,60,75)

kable(x9[k,],"latex") # tabela wybranych wartości 

par(mfrow = c(1,1))
plot(x9[,2]~lambda0, type = "l", col = "red", ylim = c(0,1), ylab = "Frakcja odkryć")
lines(x9[,3]~lambda0, type = "l", col = "orange")
lines(x9[,4]~lambda0, type = "l", col = "green2")
lines(x9[,5]~lambda0, type = "l", col = "purple")
lines(x9[,6]~lambda0, type = "l", col = "cyan2")
lines(x9[,7]~lambda0, type = "l", col = "slategray3")
abline(h=0.05, lty =3)
legend(x = 12.08, y= 0.65, c("Bonf.","F","KS","CvM","AD","HC"), col = c("red","orange","green2","purple","cyan2","slategray3"), lwd = 2)

#Many small effects ze zmiennym n

l <- 500
e <- 15
n <- c(10,20,50,100,250,500,1000,1500,2500,4000,6000,10000)
l0 <- 5
alpha  <- 0.05
sbonf <- c()
sf <- c()
sks <- c()
sad <- c()
scvm <-c()
shc <- c()

f10 <- function(s) 
{
  for(n0 in n)
  {
    
    n1 <- floor(s*n0)
    e1 <- e/(s*n0)
    bonf <- c()
    KS <- c()
    F <- c()
    AD <- c()
    cvm<-c()
    hc <- c()
    for(razy in 1:1000)
    {
      k <- sample(seq(1,n0),n1)
      X <- rpois(n0,l*l0)
      X[k] <- rpois(n1,l*(l0+e1))
      x <- 1- ppois(X,l*l0)
      c <- cvm.test(x, 'punif')
      cvm <- c(cvm, c$p.value < alpha )
      hc <- c(hc, HCdetection(x,alpha)$H)
      bonf <- c(bonf,min(x) <= (alpha/n0))
      F <- c(F,-sum(2*log(x))>qchisq(1-alpha,2*n0))
      a = ks.test(x,'punif')
      KS <- c(KS,a$p.value < alpha)
      b = ad.test(x, null = 'punif') 
      AD <- c(AD,b$p.value < alpha)
    }
    sbonf <- c(sbonf,sum(bonf))
    sf <- c(sf,sum(F))
    sks <- c(sks,sum(KS))
    sad <- c(sad,sum(AD))
    shc <- c(shc, sum(hc))
    scvm <- c(scvm, sum(cvm))
  }
  x = data.frame( "n" = n,"Bonferroniego" = sbonf/1000, "Fishera" =  sf/1000, "KS" = sks/1000, "CvM" = scvm/1000, "AD" = sad/1000, "HC" = shc/1000)
  return(x)
}

x10 <- f10(0.5)

k <- c(1,2,4,6,7,10,12)

kable(x10[k,],"latex") # tabela wybranych wartości 
par(mfrow = c(1,1))

plot(x10[,2]~n, type = "l", col = "red", ylim = c(0,1), ylab = "Frakcja odkryć")
lines(x10[,3]~n, ylim = c(0,1), type = "l", col = "orange", main = "F")
lines(x10[,4]~n, ylim = c(0,1), type = "l", col = "green2", main = "KS")
lines(x10[,5]~n,  ylim = c(0,1),type = "l", col = "purple", main = "CvM")
lines(x10[,6]~n, ylim = c(0,1), type = "l", col = "cyan2", main = "AD")
lines(x10[,7]~n, ylim = c(0,1), type = "l", col = "slategray3",main = "HC")
legend(x =16000,y = 0.5, c("Bonf.","F","KS","CvM","AD","HC"), col = c("red","orange","green2","purple","cyan2","slategray3"), lwd = 2)
abline(h=0.05, lty =3)

#Many small effects dla HC:

l <- 500
e <- 15
n <- seq(100,500,20)
l0 <- 5
alpha  <- 0.05
sbonf <- c()
sf <- c()
sks <- c()
sad <- c()
scvm <-c()
shc <- c()

f10 <- function(s) 
{
  for(n0 in n)
  {
    
    n1 <- floor(s*n0)
    e1 <- e/(s*n0)
    hc <- c()
    for(razy in 1:1000)
    {
      k <- sample(seq(1,n0),n1)
      X <- rpois(n0,l*l0)
      X[k] <- rpois(n1,l*(l0+e1))
      x <- 1- ppois(X,l*l0)
      hc <- c(hc, HCdetection(x,alpha)$H)
    }
    shc <- c(shc, sum(hc))
  }
  x = shc/1000
  return(x)
}

x10 <- f10(0.5)
par(mfrow = c(1,1))

plot(x10~n, type = "l", col = "slategray3", ylim = c(0,1), ylab = "Moc", main = "HC")
abline(h=0.05, lty =3)


# test własny
#zal. od lambda0
alpha <- 0.05
lambda0 <- seq(0.1,20,0.1)
l <- 500
funk <- function(n)
{
  wynik <- c()
  for(l0 in lambda0)
  {
    il <- c()
    il2 <- c()
    for(i in 1:1000)
    {
      X <- rpois(n,l *l0)
      p <- qpois(1-alpha,l *n*l0)
      il <- c(il,sum(X) > p) #odrzucamy dla il = 1
    }
    wynik <- c(wynik, sum(il)/1000)
  }
  return(wynik)
}

x <- funk(1000)


# zal. od n
alpha <- 0.05
l0 <- 5
n <- c(10,20,50,100,250,500,1000,1500,2500,4000,6000,10000)
l <- 500
funk2 <- function(l0)
{
  wynik <- c()
  for(n0 in n)
  {
    il <- c()
    il2 <- c()
    for(i in 1:1000)
    {
      X <- rpois(n0,l *l0)
      p <- qpois(1-alpha,l *n0*l0)
      il <- c(il,sum(X) > p) #odrzucamy dla il = 1
    }
    wynik <- c(wynik, sum(il)/1000)
  }
  return(wynik)
}

y <- funk2(5)

par(mfrow = c(2,1))
plot(x~lambda0, type = "l",ylim = c(0,0.1), col = "cornflowerblue", ylab = "Frakcja odkryć")
abline(h=0.05, lty =3)
plot(y~n, type = "l",ylim = c(0,0.1), col = "cornflowerblue", ylab = "Frakcja odkryć")
abline(h=0.05, lty =3)

# test własny, moc testu igła w stogu
par(mfrow = c(3,1))

alpha <- 0.05
e <- seq(0.1,4,0.1)
l0 <- 5
n <- 1000
l <- 500
funk0 <- function()
{
  wynik <- c()
  for(eps in e)
  {
    il <- c()
    il2 <- c()
    for(i in 1:1000)
    {
      k <- sample(seq(1,n),1)
      X <- rpois(n,l *l0)
      X[k] <- rpois(1,l *(l0+eps))
      p <- qpois(1-alpha,l *n*l0)
      il <- c(il,sum(X) > p) #odrzucamy dla il = 1
    }
    wynik <- c(wynik, sum(il)/1000)
  }
  return(wynik)
}

x0 <- funk0()


#zal. od lambda0
e<- 0.5
alpha <- 0.05
lambda0 <- seq(0.1,20,0.1)
l <- 500
funk <- function(n)
{
  wynik <- c()
  for(l0 in lambda0)
  {
    il <- c()
    il2 <- c()
    for(i in 1:1000)
    {
      k <- sample(seq(1,n),1)
      X <- rpois(n,l *l0)
      X[k] <- rpois(1,l *(l0+e))
      p <- qpois(1-alpha,l *n*l0)
      il <- c(il,sum(X) > p) #odrzucamy dla il = 1
    }
    wynik <- c(wynik, sum(il)/1000)
  }
  return(wynik)
}

x <- funk(1000)


#zal. od n
e <- 0.5
alpha <- 0.05
l0 <- 5
n <- c(10,20,50,100,250,500,1000,1500,2500,5000)
l <- 500
funk2 <- function(l0)
{
  wynik <- c()
  for(n0 in n)
  {
    il <- c()
    il2 <- c()
    for(i in 1:1000)
    {
      k <- sample(seq(1,n0),1)
      X <- rpois(n0,l *l0)
      X[k] <- rpois(1,l *(l0+e))
      p <- qpois(1-alpha,l *n0*l0)
      il <- c(il,sum(X) > p) #odrzucamy dla il = 1
    }
    wynik <- c(wynik, sum(il)/1000)
  }
  return(wynik)
}

y <- funk2(5)

par(mfrow = c(3,1))
e <- seq(0.1,4,0.1)
plot(x0~e, type = "l",ylim = c(0,0.5), col = "cornflowerblue", ylab = "Frakcja odkryć")
abline(h=0.05, lty =3)
plot(x~lambda0, type = "l",ylim = c(0,0.5), col = "cornflowerblue", ylab = "Frakcja odkryć")
abline(h=0.05, lty =3)
plot(y~n, type = "l",ylim = c(0,0.5), col = "cornflowerblue", ylab = "Frakcja odkryć")
abline(h=0.05, lty =3)

# test własny, moc testu many small
par(mfrow = c(3,1))

alpha <- 0.05
e <- seq(1,30,1)
l0 <- 5
n <- 1000
l <- 500
funk0 <- function()
{

  wynik <- c()
  for(eps in e)
  {  
    n1 <- floor(0.5*n)
    e1 <- eps/(0.5*n)
    il <- c()
    il2 <- c()
    for(i in 1:1000)
    {
      k <- sample(seq(1,n),n1)
      X <- rpois(n,l *l0)
      X[k] <- rpois(n1,l *(l0+e1))
      p <- qpois(1-alpha,l *n*l0)
      il <- c(il,sum(X) > p) #odrzucamy dla il = 1
    }
    wynik <- c(wynik, sum(il)/1000)
  }
  return(wynik)
}

x0 <- funk0()


#zal. od lambda0
eps<- 15
alpha <- 0.05
lambda0 <- seq(0.2,15,0.2)

l <- 500
funk <- function(n)
{
  wynik <- c()
  for(l0 in lambda0)
  {
    
    n1 <- floor(0.5*n)
    e1 <- e/(0.5*n)
    il <- c()
    il2 <- c()
    for(i in 1:1000)
    {
      k <- sample(seq(1,n),n1)
      X <- rpois(n,l *l0)
      X[k] <- rpois(n1,l *(l0+e1))
      p <- qpois(1-alpha,l *n*l0)
      il <- c(il,sum(X) > p) #odrzucamy dla il = 1
    }
    wynik <- c(wynik, sum(il)/1000)
  }
  return(wynik)
}

x <- funk(1000)


#zal. od n
e <- 15
alpha <- 0.05
l0 <- 5
n <- c(10,20,50,100,250,500,1000,2500,5000)
l <- 500
funk2 <- function(l0)
{
  wynik <- c()
  for(n0 in n)
  {
    e1 <- e/(0.5*n0)
    n1 <- floor(0.5*n0)
    il <- c()
    il2 <- c()
    for(i in 1:1000)
    {
      k <- sample(seq(1,n0),n1)
      X <- rpois(n0,l *l0)
      X[k] <- rpois(n1,l *(l0+e1))
      p <- qpois(1-alpha,l *n0*l0)
      il <- c(il,sum(X) > p) #odrzucamy dla il = 1
    }
    wynik <- c(wynik, sum(il)/1000)
  }
  return(wynik)
}

y <- funk2(5)

par(mfrow = c(3,1))

n <- c(10,20,50,100,250,500,1000,2500,5000)
e <- seq(1,30,1)
plot(x0~e, type = "l",ylim = c(0,1), col = "cornflowerblue", ylab = "Frakcja odkryć")
abline(h=0.05, lty =3)
plot(x~lambda0, type = "l",ylim = c(0,1), col = "cornflowerblue", ylab = "Frakcja odkryć")
abline(h=0.05, lty =3)
plot(y~n, type = "l",ylim = c(0,1), col = "cornflowerblue", ylab = "Frakcja odkryć")
abline(h=0.05, lty =3)
