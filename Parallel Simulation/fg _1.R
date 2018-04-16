library(numDeriv)  
library(cmprsk) 
library(survminer)
library(crskdiag)
library(survival)
library(flexsurv)

#Gompertz Distribution
a1 <- 0.25
a2 <- 0.5
b1 <- 0.15
b2 <- 0.175

gamma1_0 <- function(t,a1,b1) b1*exp(a1*t)
gamma2_0 <- function(t,a2,b2) b2*exp(a2*t)

gamma1_x <- function(t,xi,X,a1,b1) gamma1_0(t,a1,b1)*exp(sum(xi*X))
gamma2_x <- function(t,xi,X,a2,b2) gamma2_0(t,a2,b2)*exp(sum(xi*X))
GAMMA1_X <- function(t,xi,X,a1,b1) (b1/a1)*(exp(a1*t)-1)*exp(sum(xi*X))
GAMMA2_X <- function(t,xi,X,a2,b2) (b2/a2)*(exp(a2*t)-1)*exp(sum(xi*X))

xi <- log(c(1.15,0.9,2))
RR <- 10
n <- 1000

fg_1 <- function(RR,n,xi,gamma1_0,gamma2_0,gamma1_x,gamma2_x,GAMMA1_X,GAMMA2_X,a1,a2,b1,b2){
  
  nx <- length((xi))
  X1 <- runif(n,0,3)
  X2 <- pmax(pmin(rnorm(n,5,1),5),0)
  X3 <- sample(0:1,n,repl=T)
  XX <- cbind(X1,X2,X3)
  
  
  rate1 <- as.vector(b1*exp(XX%*%xi))
  rate2 <- as.vector(b2*exp(XX%*%xi))
  
  start_time <- Sys.time()
  time1 <- sapply(rate1, function(x) rgompertz(1,rate = rate1, shape = a1))
  time2 <- sapply(rate2, function(x) rgompertz(1,rate = rate2, shape = a2))
  Evtime <- pmin(time1, time2)
  end_time <- Sys.time()
  
  Evstat <- ifelse(time1<time2,1,2)
  
  
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
