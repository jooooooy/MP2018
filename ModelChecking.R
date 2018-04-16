library(numDeriv)  
library(cmprsk) 
library(survminer)
library(crskdiag)

# Data dat



# Data Visualization 

fit <- cuminc(dat$Time,dat$Stat)
ggcompetingrisks(fit,
                 palette = "Dark2",
                 legend = "top",
                 ggtheme = theme_bw(), 
                 conf.int = TRUE)+ 
  geom_line(size=1.4)

fit2 <- survfit(coxph(Surv(Time, Stat == 1) ~ X1+X2+X3, data = dat))

ggsurvplot(fit2, 
           data = dat, 
           censor.shape="|", 
           censor.size = 4,
           conf.int = TRUE,             
           risk.table = 'nrisk_cumcensor',
           risk.table.height = 0.25,
           ggtheme = theme_bw(),
           surv.median.line='hv')

# Model Checking for Cox

fit3 <- coxph(Surv(Time, Stat == 1) ~ X1+X2+X3, data = dat)
cox.zph.fit <- cox.zph(fit3) #A p-value less than 0.05 indicates a violation of the proportionality assumption.
cox.zph.fit
# plot all variables
ggcoxzph(cox.zph.fit,
         font.main = 12,
         font.x = 12,
         font.y = 12,
         point.alpha = 0.3,
         caption = "Graphical Test of Proportional Hazards Assumption")+ 
  geom_line(size=2)

# Model Chechking for Fine and Gray

out2 <- diag_crr(Crsk(Time,Stat)~X1+X2+X3,data=dat,test="prop")
print(out2)
plot(out2)

out1 <- diag_crr(Crsk(Time,Stat)~X1+X2+X3,data=dat,test="lin")
print(out1)
plot(out1)

