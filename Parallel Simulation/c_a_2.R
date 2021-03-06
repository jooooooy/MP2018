library(numDeriv)  
library(cmprsk) 
library(survminer)
library(crskdiag)
library(survival)


cox_a2 <- function(RR,n, xi,lambda1_0,lambda2_0,lambda1_x,lambda2_x){
  
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
    while(Evstat[j]==0) # repeat until an event was observed
    {
      Ti <- Ti+1
      # Probability of an event for individual j at time Ti
      Prob <- lambda1_x(t=Ti,xi=xi,X=XX[j,]) + lambda2_x(t=Ti,xi=xi,X=XX[j,])
      # Determine, if any event happened at Ti
      Event <- sample(0:1,1,prob=c(1-Prob,Prob))
      # if an event happened, determine which type
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
  
  test.cox
  test.crr
  
  TIME <- end_time - start_time
  COX.1 <- test.cox$coefficients[,2]
  CICOX <- test.cox$coefficients[,5]<0.05
  CRR <- test.crr$coef[,2]
  CICRR <- test.crr$coef[,5]<0.05
  
  
  return(c(TIME,COX.1,CRR,as.numeric(CICOX),as.numeric(CICRR))) 
  
}

# testing:
cox_a2(RR=1,n=1000,xi=xi,lambda1_0=lambda1_0,lambda2_0=lambda2_0,lambda1_x=lambda1_x,lambda2_x=lambda2_x)


# let's parallel

library(parallel)
library(doSNOW)

lambda1_0 <- function(t=0) 0.06
lambda2_0 <- function(t=0) 0.01
lambda1_x <- function(t=0,xi=xi,X) lambda1_0(t)*exp(sum(xi*X))
lambda2_x <- function(t=0,xi=xi,X) lambda2_0(t)*exp(sum(xi*X))
xi <- log(c(1.15,0.9,2))

clnum<-detectCores() 
clnum
cl <- makeCluster(6);
registerDoSNOW(cl)

n <- 500
set.seed(123)

n <- 1000
set.seed(456)

t1 <- Sys.time()
result <- foreach(i=1:5000,.combine=rbind,.packages=c('survival','cmprsk')) %dopar% {
  cox_a2(RR=i,n=n,xi=xi,
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


