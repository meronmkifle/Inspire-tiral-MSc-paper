

# Article title 

# The prognostic utility of soluble fms-like tyrosine kinase-1 (sFlt-1) and placental growth factor (PlGF)
# biomarkers for predicting preeclampsia: A secondary analysis of data from the INSPIRE trial data   

# libraries
pacman::p_load(readxl, tidyverse,rms,reshape2, mfp, plyr,Hmisc,MASS,pROC,DescTools,ggrepel,gridExtra,boot)

# Load dataset


Inspire <- read_excel("Inspire.xls")
View(Inspire)

attach(Inspire)

# Question-1. When adjusted for the trial arm, what are the calibration and  discrimination 
# performances of the three biomarkers as continous values?
# 1.1) sFlt-1, 
# 1.2) PIGF 
# 1.3) sflt-1/PIGF ratio 

# Question-2. When adjusted for the trial arm, what are the calibration and  discrimination 
# performances of the ratio of the biomarkers when used a continuous vs when used as a binary in predicting PE?
# 2.1) sflt-1/PIGF ratio 
# 2.2) sflt-1/PIGF cut off of 38  



# Descriptive analyisis
# age at recruitment 

Inspire %>% group_by(PE=pet_after_visit1_within7days) %>% 
  summarise(mean=mean(age_at_recruitment), sd=sd(age_at_recruitment),)


# Gestational age at recruitment

Inspire %>% group_by(PE=pet_after_visit1_within7days) %>% 
  summarise(median=median(gestweeks_recruitment), 
            quart25 = quantile(gestweeks_recruitment, 0.25),
            quart75 = quantile(gestweeks_recruitment, 0.75))

# BMI

Inspire %>% group_by(PE=pet_after_visit1_within7days) %>% 
  summarise(median=median(bmi), 
            quart25 = quantile(bmi, 0.25),
            quart75 = quantile(bmi, 0.75))


# Median and IQR of biomarker values (log transformed) by PE within 7 days
# sflt-1  by PE within 7 days

Inspire %>% group_by(PE= pet_after_visit1_within7days) %>% 
  summarise(median_sflt=median(logsflt1vi), 
            quart25 = quantile(logsflt1vi, 0.25),
            quart75 = quantile(logsflt1vi, 0.75))



# PIGF  by PE within 7 days
Inspire %>%
  group_by(PE=pet_after_visit1_within7days) %>%
  summarise(median_pigf = median(logpigfv1),
            quart25= quantile(logpigfv1, 0.25),
            quart75=quantile(logpigfv1, 0.75) )


# sflt-1/PIGF ratio  by PE within 7 days
Inspire %>%
  group_by(PE=pet_after_visit1_within7days) %>%
  summarise(median  = median(log10_sFlt1_PIGF_atvisit1),
            quart25=quantile(log10_sFlt1_PIGF_atvisit1, 0.25),
            quat75=quantile(log10_sFlt1_PIGF_atvisit1, 0.75))

# PE outcome by biomarker values (log transformed)
par(mfrow=c(1,3))
boxplot(logsflt1vi ~ pet_after_visit1_within7days, xlab = "Preeclampsia", ylab = "sFlt values (pg/mL)")
boxplot(logpigfv1 ~ pet_after_visit1_within7days, xlab = "Preeclampsia", ylab = "PIGF values (pg/mL)")
boxplot(log10_sFlt1_PIGF_atvisit1 ~ pet_after_visit1_within7days, xlab = "Preeclampsia", ylab = "sflt/PIGF ratio values (pg/mL)")

# sflt-1/PIGF cut off proportion  by PE within 7 days

Inspire %>% janitor:: tabyl(ratio_class_atvisit1,  pet_after_visit1_within7days) 

#######################################################
## 1.1) Prognostic model with (sflt-1) + (trial arm) ##
#######################################################
# Fitting regression model on PE with (sflt-1) + (trial arm) 

sflt_mod <- glm(pet_after_visit1_within7days2 ~  factor(trial_arm)+ logsflt1vi ,family="binomial")
summary(sflt_mod)


# Obtain the predicted probabilities for each patient
prob <- predict(sflt_mod,type="response")

# Obtain the linear predictor/PI for each patient
lin_pred <- predict(sflt_mod,type="link")


# Associations between predicted probabilities of PE and log transformed sFlt-1 values
plot(prob~logsflt1vi, xlab="sflt-1 (pg/mL)", ylab = "Predicted probabilities of PE",
     main= "Observed sflt-1 values Vs predicted probabilities", pch= 20)

# Apparent model performance
# Overall Model fit 
# Estimate R squared(PseudoR2)
PseudoR2(sflt_mod,"Nagelkerke")
# And the Brier score
BIC(sflt_mod)

# Model calibration
# CITL 
mod_log_1 <- glm(pet_after_visit1_within7days2~offset(lin_pred),family="binomial")
mod_log_1$coefficients


# Slope
mod_log_2<-glm(pet_after_visit1_within7days2~lin_pred,family="binomial",x=TRUE,y=TRUE)
mod_log_2$coef

# Model discriminaton
# c-index
c1 <- roc(pet_after_visit1_within7days2~prob,plot=TRUE)
c1

# Calibration plot
# create 10 risk groups
groups <- cut(lin_pred,breaks=quantile(lin_pred, prob = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(Inspire,groups,prob)
obs <- (ddply(gpdata,~groups,summarise,mean=mean(as.numeric(pet_after_visit1_within7days2)))[,2])
exp <- ddply(gpdata,~groups,summarise,mean=mean(prob))
obsn <- table(pet_after_visit1_within7days2,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# graphing a simple calibration plot over 10 risk groups
par(pty="s")
plot(obs~exp[,2],xlim=c(0,1),ylim=c(0,1),col="red", ylab="Observed",xlab="Predicted")
lines(c(0,1),c(0,1),lty=2)
for(i in 1:10){lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col="green")}

h <- hist(prob, breaks=50, plot=F)
for(i in 1:length(h$mids)){lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                                          1-((h$counts[i]/max(h$counts))/10)))}

# Add a loess smoother to the plot
obs_all <- predict(loess(pet_after_visit1_within7days2~prob,span=1))
lines_data <- data.frame(prob,obs_all)
lines_data2 <- lines_data[order(prob),] 
lines(lines_data2[,1],lines_data2[,2],col="blue")
legend(0.0,0.99988,c("Risk groups","Reference line","95% CI","Loess"),
       col=c("red","black","green","blue"),lty=c(0,2,1,1),pch=c(1,NA,NA,NA), bty="n")





# Internal validation  
# For C-statistic (Dxy/2 + 0.5) & c-slope  

mod1 <-  lrm(Inspire$pet_after_visit1_within7days  ~ logsflt1vi , x=TRUE , y=TRUE)
set.seed(231398)
boot_1 <- validate(mod1,method="boot",B=1000, pr=F)
boot_1



# Uniform adjusted final model (optimism adjusted, intercept re estimated)

# Shrinkage value = bootstrap adjusted  slope 
# Multiplying the Betas by  the shrinkage value, but the LP of a logistic model 
# includes the intercept term which must be re-estimated post shrinkage of the betas.
# Therefore we must first remove the constant/interceptt from the LP variable and
# then multiply the betas (lp) by the shrinkage factor (0.98)

lin_pred_no_constant = lin_pred - sflt_mod$coef[1]
lin_pred_no_constant_adjusted = 0.98*lin_pred_no_constant

# Now to re-estimate a new intercept value allowing for the shrunken betas
# We fit a logistic model with the shrunken lp set as an offset to fix these shrunken beta coefficients.
# We then add this additional intercept term to the linear predictor to obtain our final shrunken model LP.

mod_log2 <- glm(pet_after_visit1_within7days2~offset(lin_pred_no_constant_adjusted),family="binomial")
lin_pred_adjusted = lin_pred_no_constant_adjusted + mod_log2$coef


# CITL of adjusted model
mod_log3 <- glm(pet_after_visit1_within7days2~offset(lin_pred_adjusted),family="binomial")
mod_log3$coefficients

#BIC
BIC(mod_log3)
# Obtain the predicted probabilities for each patient
probb <- predict(mod_log3,type="response") 


# Associations between predicted probabilities of PE and log transformed sFlt-1 values

plot(probb~logsflt1vi, xlab="sflt-1 (pg/mL)", ylab = "Predicted probabilities of PE",
     main= "Observed sflt-1 values Vs predicted probabilities", pch= 20)

# Model discriminaton
# c-index
c1 <- roc(pet_after_visit1_within7days2~probb,ci=T,plot=TRUE, legacy.axes=TRUE, print.auc=TRUE)
c1



par(mfrow=c(1,3))
# Calibration plot
# create 10 risk groups
groups <- cut(probb,breaks=quantile(probb, 
     prob = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)),labels=c(1:10),include.lowest=TRUE)



h <- hist(probb, breaks=50, plot=F)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i], 
                                 1-((h$counts[i]/max(h$counts))/10)))}


# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(Inspire,groups,probb)
obs <- (ddply(gpdata,~groups,summarise,mean=mean(as.numeric(pet_after_visit1_within7days2)))[,2])
exp <- ddply(gpdata,~groups,summarise,mean=mean(probb))
obsn <- table(pet_after_visit1_within7days2,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# graph a simple calibration plot over 10 risk groups
par(pty="s")
plot(obs~exp[,2],xlim=c(0,1),ylim=c(0,1),col="red", ylab="Observed",xlab="Predicted
Sflt model")
lines(c(0,1),c(0,1),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col="green")
}

h <- hist(probb, breaks=50, plot=F)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i], 
                                 1-((h$counts[i]/max(h$counts))/10)))}

# Add a loess smoother to the plot
obs_all <- predict(loess(pet_after_visit1_within7days2~probb,span=1))
lines_data <- data.frame(probb,obs_all)
lines_data2 <- lines_data[order(probb),] 
lines(lines_data2[,1],lines_data2[,2],col="blue")
legend(0.0,0.99,c("Risk groups","Reference line","95% CI","Loess"),
       col=c("red","black","green","blue"),lty=c(0,2,1,1),pch=c(1,NA,NA,NA), bty="n")






##################################




#######################################################
## 1.2) Prognostic model with (PIGF) + (trial arm) ##
#######################################################

# Fitting regression model on PE with (PIGF) + (trial arm) 

PIGF_mod <- glm(pet_after_visit1_within7days2~  factor(trial_arm)+ logpigfv1 ,family="binomial")
summary(PIGF_mod)

# Obtain the predicted probabilities for each patient
prob <- predict(PIGF_mod,type="response")

# Obtain the linear predictor/PI for each patient
lin_pred <- predict(PIGF_mod,type="link")


# Associations between predicted probabilities of PE and log transformed sFlt-1 values
plot(prob~logpigfv1, xlab="PIGF (pg/mL)", ylab = "Predicted probabilities of PE",
     main= "Observed PIGF values Vs predicted probabilities", pch= 20)


# Apparent model performance
# Overall Model fit 
# Estimate R squared(PseudoR2)
PseudoR2(PIGF_mod,"Nagelkerke")
# And the Brier score
BIC(PIGF_mod)

# Model calibration
# CITL 
mod_log_1 <- glm(pet_after_visit1_within7days2~offset(lin_pred),family="binomial")
mod_log_1$coefficients


# Slope
mod_log_2<-glm(pet_after_visit1_within7days2~lin_pred,family="binomial",x=TRUE,y=TRUE)
mod_log_2$coef

# Model discriminaton
# c-index
c1 <- roc(pet_after_visit1_within7days2~prob,ci=TRUE)
c1


# Calibration plot
# create 10 risk groups
groups <- cut(lin_pred,breaks=quantile(lin_pred, prob = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(Inspire,groups,prob)
obs <- (ddply(gpdata,~groups,summarise,mean=mean(as.numeric(pet_after_visit1_within7days2)))[,2])
exp <- ddply(gpdata,~groups,summarise,mean=mean(prob))
obsn <- table(pet_after_visit1_within7days2,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# graphing a simple calibration plot over 10 risk groups
par(pty="s")
plot(obs~exp[,2],xlim=c(0,1),ylim=c(0,1),col="red", ylab="Observed",xlab="Predicted")
lines(c(0,1),c(0,1),lty=2)
for(i in 1:10){lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col="green")}

h <- hist(prob, breaks=50, plot=F)
for(i in 1:length(h$mids)){lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                                          1-((h$counts[i]/max(h$counts))/10)))}

# Add a loess smoother to the plot
obs_all <- predict(loess(pet_after_visit1_within7days2~prob,span=1))
lines_data <- data.frame(prob,obs_all)
lines_data2 <- lines_data[order(prob),] 
lines(lines_data2[,1],lines_data2[,2],col="blue")
legend(0.0,0.99988,c("Risk groups","Reference line","95% CI","Loess"),
       col=c("red","black","green","blue"),lty=c(0,2,1,1),pch=c(1,NA,NA,NA), bty="n")





# Internal validation  
# For C-statistic (Dxy/2 + 0.5) & c-slope  

mod1 <-  lrm(Inspire$pet_after_visit1_within7days ~ logpigfv1 , x=TRUE , y=TRUE)
set.seed(231398)
boot_1 <- validate(mod1,method="boot",B=1000, pr=F)
boot_1



# Uniform adjusted final model (optimism adjusted, intercept re estimated)

# Shrinkage value = bootstrap adjusted  slope 
# Multiplying the Betas by  the shrinkage value, but the LP of a logistic model 
# includes the intercept term which must be re-estimated post shrinkage of the betas.
# Therefore we must first remove the constant/interceptt from the LP variable and
# then multiply the betas (lp) by the shrinkage factor (0.98)

lin_pred_no_constant = lin_pred - PIGF_mod$coef[1]
lin_pred_no_constant_adjusted = 0.99*lin_pred_no_constant

# Now to re-estimate a new intercept value allowing for the shrunken betas
# We fit a logistic model with the shrunken lp set as an offset to fix these shrunken beta coefficients.
# We then add this additional intercept term to the linear predictor to obtain our final shrunken model LP.

mod_log2 <- glm(pet_after_visit1_within7days2~offset(lin_pred_no_constant_adjusted),family="binomial")
lin_pred_adjusted = lin_pred_no_constant_adjusted + mod_log2$coef


# CITL of adjusted model
mod_log3 <- glm(pet_after_visit1_within7days2~offset(lin_pred_adjusted),family="binomial")
mod_log3$coefficients

# Obtain the predicted probabilities for each patient
probb <- predict(mod_log3,type="response") 
#BIC
BIC(mod_log3)
# Associations between predicted probabilities of PE and log transformed sFlt-1 values

plot(probb~logpigfv1, xlab="PIGF (pg/mL)", ylab = "Predicted probabilities of PE",
     main= "Observed PIGF values Vs predicted probabilities", pch= 20)



# Calibration plot
# create 10 risk groups
groups <- cut(probb,breaks=quantile(probb, 
                                    prob = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)),labels=c(1:10),include.lowest=TRUE)

h <- hist(probb, breaks=50, plot=F)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i], 
                                 1-((h$counts[i]/max(h$counts))/10)))}


# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(Inspire,groups,probb)
obs <- (ddply(gpdata,~groups,summarise,mean=mean(as.numeric(pet_after_visit1_within7days2)))[,2])
exp <- ddply(gpdata,~groups,summarise,mean=mean(probb))
obsn <- table(pet_after_visit1_within7days2,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# graph a simple calibration plot over 10 risk groups
par(pty="s")
plot(obs~exp[,2],xlim=c(0,1),ylim=c(0,1),col="red", ylab="Observed",xlab="Predicted")
lines(c(0,1),c(0,1),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col="green")
}

h <- hist(probb, breaks=50, plot=F)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i], 
                                 1-((h$counts[i]/max(h$counts))/10)))}

# Add a loess smoother to the plot
obs_all <- predict(loess(pet_after_visit1_within7days2~probb,span=1))
lines_data <- data.frame(probb,obs_all)
lines_data2 <- lines_data[order(probb),] 
lines(lines_data2[,1],lines_data2[,2],col="blue")
legend(0.0,0.999,c("Risk groups","Reference line","95% CI","Loess"),
       col=c("red","black","green","blue"),lty=c(0,2,1,1),pch=c(1,NA,NA,NA), bty="n")





##############################



## 1.3) Prognostic model with (sflt-1/PIGF ratio) + (trial arm) ##
#######################################################


# Fitting regression model on PE with (sflt-1) + (trial arm) 

sflt_PIGF_mod <- glm(pet_after_visit1_within7days2~  factor(trial_arm)+ log10_sFlt1_PIGF_atvisit1 ,family="binomial")
summary(sflt_PIGF_mod)
 # BIC
BIC(sflt_PIGF_mod)


# Obtain the linear predictor/PI for each patient
lin_pred <- predict(sflt_PIGF_mod,type="link")

# Obtain the predicted probabilities for each patient
prob <- plogis(lin_pred)


# Associations between predicted probabilities of PE and log transformed sFlt-1 values
plot(prob~log10_sFlt1_PIGF_atvisit1, xlab="sflt-1 (pg/mL)", ylab = "Predicted probabilities of PE",
     main= "Observed sflt-1 values Vs predicted probabilities", pch= 20)


# Apparent model performance
# Overall Model fit 
# Estimate R squared(PseudoR2)
PseudoR2(sflt_PIGF_mod,"Nagelkerke")
# And the Brier score
BIC(sflt_PIGF_mod)

# Model calibration
# CITL 
mod_log_1 <- glm(pet_after_visit1_within7days2~offset(lin_pred),family="binomial")
mod_log_1$coefficients


# Slope
mod_log_2<-glm(pet_after_visit1_within7days2~lin_pred,family="binomial",x=TRUE,y=TRUE)
mod_log_2$coef

# Model discriminaton
# c-index
c1 <- roc(pet_after_visit1_within7days2~prob,ci=TRUE)
c1


# Calibration plot
# create 10 risk groups
groups <- cut(lin_pred,breaks=quantile(lin_pred, prob = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(Inspire,groups,prob)
obs <- (ddply(gpdata,~groups,summarise,mean=mean(as.numeric(pet_after_visit1_within7days2)))[,2])
exp <- ddply(gpdata,~groups,summarise,mean=mean(prob))
obsn <- table(pet_after_visit1_within7days2,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# graphing a simple calibration plot over 10 risk groups
par(pty="s")
plot(obs~exp[,2],xlim=c(0,1),ylim=c(0,1),col="red", ylab="Observed",xlab="Predicted")
lines(c(0,1),c(0,1),lty=2)
for(i in 1:10){lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col="green")}

h <- hist(prob, breaks=50, plot=F)
for(i in 1:length(h$mids)){lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],
                                                          1-((h$counts[i]/max(h$counts))/10)))}

# Add a loess smoother to the plot
obs_all <- predict(loess(pet_after_visit1_within7days2~prob,span=1))
lines_data <- data.frame(prob,obs_all)
lines_data2 <- lines_data[order(prob),] 
lines(lines_data2[,1],lines_data2[,2],col="blue")
legend(0.0,0.999,c("Risk groups","Reference line","95% CI","Loess"),
       col=c("red","black","green","blue"),lty=c(0,2,1,1),pch=c(1,NA,NA,NA), bty="n")





# Internal validation  
# For C-statistic (Dxy/2 + 0.5) & c-slope  

mod1 <-  lrm(pet_after_visit1_within7days2 ~ log10_sFlt1_PIGF_atvisit1 , x=TRUE , y=TRUE)
set.seed(231398)
boot_1 <- validate(mod1,method="boot",B=1000, pr=F)
boot_1



# Uniform adjusted final model (optimism adjusted, intercept re estimated)

# Shrinkage value = bootstrap adjusted  slope 
# Multiplying the Betas by  the shrinkage value, but the LP of a logistic model 
# includes the intercept term which must be re-estimated post shrinkage of the betas.
# Therefore we must first remove the constant/interceptt from the LP variable and
# then multiply the betas (lp) by the shrinkage factor (0.98)

lin_pred_no_constant = lin_pred - sflt_PIGF_mod$coef[1]
lin_pred_no_constant_adjusted = 0.94*lin_pred_no_constant

# Now to re-estimate a new intercept value allowing for the shrunken betas
# We fit a logistic model with the shrunken lp set as an offset to fix these shrunken beta coefficients.
# We then add this additional intercept term to the linear predictor to obtain our final shrunken model LP.

mod_log2 <- glm(pet_after_visit1_within7days2~offset(lin_pred_no_constant_adjusted),family="binomial")
lin_pred_adjusted = lin_pred_no_constant_adjusted + mod_log2$coef


# CITL of adjusted model
mod_log3 <- glm(pet_after_visit1_within7days2~offset(lin_pred_adjusted),family="binomial")
mod_log3$coefficients

# Obtain the predicted probabilities for each patient
probb <- predict(mod_log3,type="response") 
#BIC
BIC(mod_log3)
# Associations between predicted probabilities of PE and log transformed sFlt-1 values

plot(probb~log10_sFlt1_PIGF_atvisit1, xlab="sflt-1 (pg/mL)", ylab = "Predicted probabilities of PE",
     main= "Observed sflt/PIGF values Vs predicted probabilities", pch= 20)



# Calibration plot
# create 10 risk groups
groups <- cut(lin_pred_adjusted,breaks=quantile(lin_pred_adjusted, 
                                                probb = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)),labels=c(1:10),include.lowest=TRUE)
h <- hist(probb, breaks=50, plot=F)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i], 
                                 1-((h$counts[i]/max(h$counts))/10)))}


# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(Inspire,groups,probb)
obs <- (ddply(gpdata,~groups,summarise,mean=mean(as.numeric(pet_after_visit1_within7days2)))[,2])
exp <- ddply(gpdata,~groups,summarise,mean=mean(probb))
obsn <- table(pet_after_visit1_within7days2,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# graph a simple calibration plot over 10 risk groups
par(pty="s")
plot(obs~exp[,2],xlim=c(0,1),ylim=c(0,1),col="red", ylab="Observed",xlab="Predicted")
lines(c(0,1),c(0,1),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col="green")
}

h <- hist(probb, breaks=50, plot=F)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i], 
                                 1-((h$counts[i]/max(h$counts))/10)))}

# Add a loess smoother to the plot
obs_all <- predict(loess(pet_after_visit1_within7days2~probb,span=1))
lines_data <- data.frame(probb,obs_all)
lines_data2 <- lines_data[order(probb),] 
lines(lines_data2[,1],lines_data2[,2],col="blue")
legend(0.0,0.999,c("Risk groups","Reference line","95% CI","Loess"),
       col=c("red","black","green","blue"),lty=c(0,2,1,1),pch=c(1,NA,NA,NA), bty="n")






##################



## 1.4) Prognostic model with (sflt-1/PIGF cuttoff) + (trial arm) ##
#######################################################


# Fitting regression model on PE with (sflt-1) + (trial arm) 

cuttoff_mod <- glm(pet_after_visit1_within7days2~  factor(trial_arm)+ ratio_class_atvisit1 ,family="binomial")
summary(cuttoff_mod)
BIC(cuttoff_mod)


# Obtain the linear predictor/PI for each patient
lin_pred <- predict(cuttoff_mod,type="link")

# Obtain the predicted probabilities for each patient
prob <- plogis(lin_pred)

# Associations between predicted probabilities of PE and log transformed sFlt-1 values
plot(prob~cuttoff, xlab="sflt-1 (pg/mL)", ylab = "Predicted probabilities of PE",
     main= "Observed sflt/PIGF cuttoff values Vs predicted probabilities", pch= 20)


# Apparent model performance
# Overall Model fit 
# Estimate R squared(PseudoR2)
PseudoR2(cuttoff_mod,"Nagelkerke")
# And the Brier score
BIC(cuttoff_mod)

# Model calibration
# CITL 
mod_log_1 <- glm(pet_after_visit1_within7days2~offset(lin_pred),family="binomial")
mod_log_1$coefficients


# Slope
mod_log_2<-glm(pet_after_visit1_within7days2~lin_pred,family="binomial",x=TRUE,y=TRUE)
mod_log_2$coef

# Model discriminaton
# c-index
c1 <- roc(pet_after_visit1_within7days2~prob,ci=TRUE)
c1




# Internal validation  
# For C-statistic (Dxy/2 + 0.5) & c-slope  

mod1 <-  lrm(pet_after_visit1_within7days2 ~ ratio_class_atvisit1 , x=TRUE , y=TRUE)
set.seed(231398)
boot_1 <- validate(mod1,method="boot",B=1000, pr=F)
boot_1



# Uniform adjusted final model (optimism adjusted, intercept re estimated)

# Shrinkage value = bootstrap adjusted  slope 
# Multiplying the Betas by  the shrinkage value, but the LP of a logistic model 
# includes the intercept term which must be re-estimated post shrinkage of the betas.
# Therefore we must first remove the constant/interceptt from the LP variable and
# then multiply the betas (lp) by the shrinkage factor (0.98)

lin_pred_no_constant = lin_pred - cuttoff_mod$coef[1]
lin_pred_no_constant_adjusted = 0.88*lin_pred_no_constant

# Now to re-estimate a new intercept value allowing for the shrunken betas
# We fit a logistic model with the shrunken lp set as an offset to fix these shrunken beta coefficients.
# We then add this additional intercept term to the linear predictor to obtain our final shrunken model LP.

mod_log2 <- glm(pet_after_visit1_within7days2~offset(lin_pred_no_constant_adjusted),family="binomial")
lin_pred_adjusted = lin_pred_no_constant_adjusted + mod_log2$coef


# CITL of adjusted model
mod_log3 <- glm(pet_after_visit1_within7days2~offset(lin_pred_adjusted),family="binomial")
mod_log3$coefficients
# Obtain the predicted probabilities for each patient
probb <- predict(mod_log3,type="response") 
#BIC
BIC(mod_log3)
# Associations between predicted probabilities of PE and log transformed sFlt-1 values

plot(probb~logsflt1vi, xlab="sflt-1 (pg/mL)", ylab = "Predicted probabilities of PE",
     main= "Observed sflt-1 values Vs predicted probabilities", pch= 20)









#########################  Delong test for statistical differnce of models AUC curves


df <- data.frame(disease_status = rbinom(n=100, size=1, prob=0.20),
                 test1 = rnorm(100, mean=15, sd=4),
                 test2 = rnorm(100, mean=30, sd=2),
                 test3 = rnorm(100, mean=50, sd=3))

#create roc object for test1, test2, test3
roc.out_test1<-roc(df$disease_status, df$test1, plot=TRUE, smooth = FALSE)
#> Setting levels: control = 0, case = 1
#> Setting direction: controls < cases

roc.out_test2<-roc(df$disease_status, df$test2, plot=TRUE, smooth = FALSE)
#> Setting levels: control = 0, case = 1
#> Setting direction: controls < cases

roc.out_test3<-roc(df$disease_status, df$test3, plot=TRUE, smooth = FALSE)
#> Setting levels: control = 0, case = 1
#> Setting direction: controls < cases

# compare the AUC of test1 and test 2
roc.test(roc.out_test1, roc.out_test2,  method = "delong", na.rm = TRUE)
#> 
#>  DeLong's test for two correlated ROC curves
#> 
#> data:  roc.out_test1 and roc.out_test2
#> Z = 0.60071, p-value = 0.548
#> alternative hypothesis: true difference in AUC is not equal to 0
#> sample estimates:
#> AUC of roc1 AUC of roc2 
#>   0.5840108   0.5216802
