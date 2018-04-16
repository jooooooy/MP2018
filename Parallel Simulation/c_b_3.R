library(numDeriv)  
library(cmprsk) 
library(survminer)
library(crskdiag)
library(survival)


cox_b3 <- function(RR,n, xi,lambda1_0,lambda2_0,lambda1_x,lambda2_x){
  
  nx <- length((xi))
  X1 <- runif(n,0,3)
  X2 <- pmax(pmin(rnorm(n,5,1),10),0)
  X3 <- sample(0:1,n,repl=T)
  XX <- cbind(X1,X2,X3)
  
  Evstat <- rep(0,n)
  Evtime <- c()
  
  start_time <- Sys.time() 
  
  U <- runif(n,0,1)
  J <- 1:n
  CDF <- function(t,xi,X,uni,j){
    L1 <- 0.09*log(t+1)*exp(sum(xi*X[j,]))
    L2 <- 0.012*t*t*exp(sum(xi*X[j,]))
    return(1-exp(-L1-L2)-uni)
  } 
  
  # sapply(U, function(x) uniroot(CDF,uni=x,xi=xi,X=XX,j=2,c(1e-6,1e6)))
  ANS <- mapply(function(x,counts) uniroot(CDF,uni=x,xi=xi,X=XX,j=counts,c(0,1e6),extendInt = 'yes'),U,J)
  Evtime <- unlist(ANS[1,])
  Evstat <- sapply(J, function(j) sample(1:2,1,
                                         prob = c(lambda1_x(t=Evtime[j],xi,X=XX[j,]),
                                                  lambda2_x(t=Evtime[j],xi,X=XX[j,]))/(lambda1_x(t=Evtime[j],xi,X=XX[j,]) + lambda2_x(t=Evtime[j],xi,X=XX[j,]))))
  
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
cox_b3(RR=1,n=1000,xi=xi,
       lambda1_0=lambda1_0,lambda2_0=lambda2_0,lambda1_x=lambda1_x,lambda2_x=lambda2_x)



lambda1_0 <- function(t=0) 0.09/(t+1)
lambda2_0 <- function(t=0) 0.024*t
lambda1_x <- function(t=0,xi=xi,X) lambda1_0(t)*exp(sum(xi*X))
lambda2_x <- function(t=0,xi=xi,X) lambda2_0(t)*exp(sum(xi*X))


xi <- log(c(1.15,0.9,2))


# let's parallel

library(parallel)
library(doSNOW)
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
  cox_b3(RR=i,n=n,xi=xi,
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
