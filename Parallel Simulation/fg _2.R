library(numDeriv)  
library(cmprsk) 
library(survminer)
library(crskdiag)
library(survival)

gamma_0 <- function(t) 0.001 * exp(-0.001*t/log(1.3))
lambda1_0 <- function(t) 0.001
lambda2_0 <- function(t) gamma_0(t)-lambda1_0(t) + 0.001/ln(1.3) 

gamma_x <- function(t,xi,X) gamma_0(t)*exp(sum(xi*X))
lambda1_x <- function(t,xi,X)
  lambda1_0(t)*exp(sum(xi*X)*exp(-0.0005*t))


xi <- log(c(1.15,0.9,2))
# Alternatively – to save computation time – the cause-specific hazard
# for the competing event can be determined analytically
lambda2_x <- function(t,xi,X)
  0.001*exp(-0.001*t/log(1.3)+sum(xi*X)) -
  0.001*exp(sum(xi*X*exp(-0.0005*t))) + 0.001/log(1.3) -
  sum(xi*X*0.0005*exp(-0.0005*t))

fg_2 <- function(RR,n, xi,lambda1_0,lambda2_0,lambda1_x,lambda2_x){
  
  nx <- length((xi))
  X1 <- runif(n,0,3)
  X2 <- pmax(pmin(rnorm(n,5,1),10),0)
  X3 <- sample(0:1,n,repl=T)
  XX <- cbind(X1,X2,X3)
  
  Evstat <- rep(0,n)
  Evtime <- c()
  
  start_time <- Sys.time() 
  for(j in 1:n)
  {
    Ti <- 0 
    while(Evstat[j]==0) 
    {
      Ti <- Ti+1
      Prob <- lambda1_x(t=Ti,xi=xi,X=XX[j,]) + lambda2_x(t=Ti,xi=xi,X=XX[j,])
      Event <- sample(0:1,1,prob=c(1-Prob,Prob))
      if(Event==1)
        Evstat[j] <- sample(1:2,1,prob=c(
          lambda1_x(t=Ti,xi=xi,X=XX[j,]),
          lambda2_x(t=Ti,xi=xi,X=XX[j,]))/(lambda1_x(t=Ti,xi,X=XX[j,]) + lambda2_x(t=Ti,xi,X=XX[j,])))
    }
    Evtime[j] <- Ti
  }
  end_time <- Sys.time()
  
  dat <- data.frame(Time=Evtime,Stat=Evstat,X1=X1,X2=X2,X3=X3)
  
  test.cox <- summary(coxph(Surv(Time, Stat == 1) ~ X1+X2+X3,data=dat))
  test.crr <- summary(crr(dat$Time,dat$Stat,cov1=cbind(dat$X1,dat$X2,dat$X3)))
  
  TIME <- start_time - end_time
  COX.1 <- test.cox$coefficients[,2]
  CICOX <- test.cox$coefficients[,5]<0.05
  CRR <- test.crr$coef[,2]
  CICRR <- test.crr$coef[,5]<0.05
  
  return(c(TIME,COX.1,CRR,as.numeric(CICOX),as.numeric(CICRR))) 
  
}

# testing:
fg_2(RR=1,n=1000,xi=xi,lambda1_0=lambda1_0,lambda2_0=lambda2_0,lambda1_x=lambda1_x,lambda2_x=lambda2_x)



# let's parallel

library(parallel)
library(doSNOW)
clnum<-detectCores() 
clnum
cl <- makeCluster(8);
registerDoSNOW(cl)

n <- 500
set.seed(123)

n <- 1000
set.seed(456)

t1 <- Sys.time()
result <- foreach(i=1:5000,.combine=rbind,.packages=c('survival','cmprsk')) %dopar% {
  fg_2(RR=i,n=n,xi=xi,
         lambda1_0=lambda1_0,lambda2_0=lambda2_0,
         lambda1_x=lambda1_x,lambda2_x=lambda2_x)
}
t2 <- Sys.time()
t2-t1

# c(TIME,COX.1,CRR,as.numeric(CICOX),as.numeric(CICRR))
TIME <- result[,1] 
COX.1 <- result[,2:4] 
CRR <- result[,5:7] 
CICOX <- result[,8:10] 
CICRR <- result[,11:13] 

stopCluster(cl)

# Finally!! HERE IS THE SUMMARY STATISTICS
exp(xi)

mean(TIME)

(MEAN.COX.1 <- apply(COX.1,2,mean,na.rm=T))
(MEAN.CRR <- apply(CRR,2,mean,na.rm=T))

(MED.COX.1 <- apply(COX.1,2,median,na.rm=T))
(MED.CRR <- apply(CRR,2,median,na.rm=T))

(SD.COX.1 <- apply(COX.1,2,sd,na.rm=T))
(SD.CRR <- apply(CRR,2,sd,na.rm=T))

(Q1.COX.1 <- apply(COX.1,2,quantile,probs=c(0.025,0.975),na.rm=T))
(Q1.CRR <- apply(CRR,2,quantile,probs=c(0.025,0.975),na.rm=T))

(CI.COX.1 <- apply(CICOX,2,mean,na.rm=T))
(CI.CRR <- apply(CICRR,2,mean,na.rm=T))
