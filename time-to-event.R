remove(list=ls())
set.seed(42) # the answer to life the universe and everything

# simulate data
library(MASS)
N <- 500
Sigma <- outer(1:10, 1:10, function(x,y) 0.5^abs(x-y)) #variance covariance matrix for covariates

x <- mvrnorm(n = N, rep(0, 10), Sigma)
x[,3] <- ifelse(x[,3] > 0.5, 1, 0) #binary predictor
x[,4] <- ifelse(x[,4] > 0, 1, 0) #binary predictor
x[,5] <- ifelse(x[,5] > 1, 1, 0) #binary predictor
x[,8] <- ifelse(x[,8] > 2, 1, 0) #binary auxiliary variable
x[,9] <- ifelse(x[,9] > 2, 1, 0) #binary auxiliary variable
x[,10] <- cut(x[,10], breaks=c(-Inf, -1, 1, 2, Inf)) #categorical with 4 categories

survdat.compl <- data.frame(x)
colnames(survdat.compl) <- paste0("x", 1:10)
rate <- with(survdat.compl, exp(x1+0.4*x1^2+0.4*x2+0.05*x2^2+(x3==2)+0.5*(x4==2)-
            0.5*(x5==1)+rnorm(N,0,0.1))/10)
survdat.compl[,c(3:5, 8:10)] <- lapply(survdat.compl[,c(3:5, 8:10)], factor)

fulltime <- rexp(N, rate = rate)
censtimes <- 5 + 20*runif(N)
mean(censtimes)
survdat.compl$time <- pmin(fulltime, censtimes)
survdat.compl$status <- as.numeric(censtimes > fulltime)
colnames(survdat.compl)[6:10]=c("z1","z2" ,"z3", "z4", "z5")
table(survdat.compl$status)


# introduce missing data for the covariates
missing.matrix=matrix(0, nrow=nrow(survdat.compl), ncol=10)
missing.matrix=matrix(rbinom(length(missing.matrix),1, p=0.1), nrow=nrow(survdat.compl))
survdat=survdat.compl[,1:10]
survdat[missing.matrix==1]=NA
survdat=cbind(survdat, survdat.compl[11:12])

# create clusters
survdat$clust <- factor(sample(1:5, size = N, replace = TRUE, prob = rep(0.2,5)))
survdat <- survdat[order(survdat$clust),]

head(survdat)

# KM curves
library(ggplot2)
library(ggfortify)
library(survival)
model_fit <- survfit(Surv(time, status) ~ 1, data=survdat) 

autoplot(model_fit) + 
  labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
       title = "Survival Times") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(face="bold", size = 12),
        axis.title.y = element_text(face="bold", size = 12),
        legend.title = element_text(face="bold", size = 10))

### plot survival curves for specific levels of covariates 
model_fit2 <- survfit(Surv(time, status) ~ x5, data=survdat) 
autoplot(model_fit2) + 
  labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
       title = "Survival Times") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(face="bold", size = 12),
        axis.title.y = element_text(face="bold", size = 12),
        legend.title = element_text(face="bold", size = 10))+ xlim(0,30)

# Impute missing data
library(Hmisc)
library(mice)
n.impute <- 10
survdat$nelsonaalen <- nelsonaalen(survdat, time, status) 
a <- aregImpute(data=survdat, ~status+x1+x2+I(x3)+I(x4)+I(x5)+
                  z1+z2+I(z3)+I(z4)+I(z5)+clust+nelsonaalen, n.impute=n.impute, nk=3, match='closest')

# get imputed datasets
imputed <- list()
for (i in 1:n.impute){
  imputed[[i]] <- impute.transcan(a, imputation=i, data=survdat, list.out=TRUE, pr=FALSE, check=FALSE)
  imputed[[i]]$time<-survdat$time
}



# Weibull model with smoothing splines
fit.model <- list()
for(i in 1:n.impute){
  fit.model[[i]] <- survreg(Surv(time, status)~pspline(x1, df=2)+pspline(x2, df=2)+ridge(x3,x4, scale=T), 
                            imputed[[i]], dist="weibul")
}

# predict for new patient
new.patient <- data.frame(x1 = 0.5, x2 = 0.2, x3 = as.factor(1), x4 = as.factor(2))
prediction.weibull.time <- function(new.patient, fit.model = NULL){ 
  pct <- 1:99/100  
  predicted.exp <- list()
  for(i in 1:n.impute){
    ptime <- predict(fit.model[[i]], newdata=new.patient, type='quantile', p=pct, se=TRUE)
    predicted.exp[[i]] <- data.frame(surv=pct, pred=ptime$fit,var=ptime$se.fit^2)
  }
  
  predicted.exp.average <- Reduce("+", predicted.exp) / n.impute
  
  means <- predicted.exp.average[,2]
  U <- predicted.exp.average[,3]
  
  errors <- list()
  for(i in 1:n.impute){
    errors[[i]] <- (predicted.exp[[i]]$pred-means)^2
  }
  B <- 1/(n.impute-1)*(Reduce("+", errors))
  varMI <- U + (1+1/n.impute)*B
  results <- data.frame("percentiles.surv"=pct,"means"=means, "sd" = sqrt(varMI), 
                        "lowerCI"=means-1.96*sqrt(varMI), "upperCI"=means+1.96*sqrt(varMI))
  return(results)
}

predictions.time <- prediction.weibull.time(new.patient, fit.model)

ggplot(predictions.time, aes(x = means, y =1 - percentiles.surv )) + 
  geom_line(col="black") + 
  geom_ribbon(aes(xmin = lowerCI, xmax = upperCI), 
              alpha=0.1, 
              linetype="dashed",
              color="grey")+xlab("days")+ylab("Survival probability")+xlim(0,30)


# calculate the linear predictor
prediction.weibull <- function(new.patient, single.fit = NULL, multiple.fit = NULL){ 
  if(!is.null(multiple.fit)){
    
    ff <- function(i){
      predict(multiple.fit[[i]], newdata = new.patient, type = "lp")
    }
    prediction_matrix <- sapply(1:length(multiple.fit), ff)
    
    if(dim(new.patient)[1] == 1){
      prediction <- mean(prediction_matrix)
    } else{
      prediction <- apply(prediction_matrix, 1, mean)
    }
    
  } else if(!is.null(single.fit)){
    prediction <- predict(single.fit, newdata = new.patient, type = "lp")
  }
  return(prediction)
}

prediction.weibull(new.patient, multiple.fit = fit.model)

complete.data <- survdat[complete.cases(survdat[,1:5]),]
predicted.weibull <- prediction.weibull(complete.data, multiple.fit = fit.model)

# Apparent performance (no optimism correction) -------------------
calculate_performance <- function(time = NULL, status = NULL, lp = NULL){
  
  #discrimination
  harrell_C <- concordance(Surv(time, status) ~ lp)
  harrell_C_est <- harrell_C$concordance
  #harrell_C_var <- harrell_C$var
  
  Uno_C <- concordance(Surv(time, status) ~ lp, timewt = "n/G2")
  Uno_C_est <- Uno_C$concordance
  #Uno_C_var <- Uno_C$var
  
  #calibration
  second_model <- survreg(Surv(time, status) ~ lp, dist="weibul")
  calslope <- second_model$coef[2]

  returnVec <- c(harrell_C_est, Uno_C_est, calslope)
  names(returnVec) <- c("Harrel_C", "Uno_C", "calibration.slope")

  return(returnVec)  
}  
apparent.weibull <- calculate_performance(time = complete.data$time, status = complete.data$status,
                      lp = predicted.weibull)
round(apparent.weibull,2)

# calibration in the large
timepoint <- 5
# Observed
obj <- summary(survfit(
  Surv(time, status) ~ 1, 
  data = complete.data),
  times = timepoint)
obs_t <- 1 - obj$surv

# Predicted risk 
prediction.weibull2 <- function(new.patient, fit.model = NULL, time.point = NULL){ 
  
  if(is.null(time.point)){
    pct <- 10*1:300/100
  } else{
    pct <- time.point
  }
  predicted.exp <- list()
  for(i in 1:n.impute){
    mu_hat <- predict(fit.model[[i]], newdata = new.patient, type = "lp")
    predicted.exp[[i]] <- 1 - pweibull(pct, shape = 1/fit.model[[i]]$scale, 
                                       scale = exp(mu_hat))
  }
  predicted.exp.average <- Reduce("+", predicted.exp) / n.impute
  
  results <- data.frame(time = pct, means = predicted.exp.average)
  return(results)
}

predictions <- list()
for (i in 1:dim(complete.data)[1]){
  predictions[[i]] <- prediction.weibull2(complete.data[i,], fit.model, time.point = timepoint)$means
}
est.surv <- unlist(predictions)
pred <- 1 - est.surv

# Expected
exp_t <- mean(pred)
OE_t <- obs_t / exp_t

alpha <- .05
OE_summary <- c(
  "OE" = OE_t,
  "2.5 %" = OE_t * exp(-qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event)),
  "97.5 %" = OE_t * exp(+qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event))
)
round(OE_summary,2)

# drawing a calibration plot
library(rms)
predictions <- list()
for (i in 1:dim(complete.data)[1]){
  predictions[[i]] <- prediction.weibull2(complete.data[i,], fit.model, time.point = complete.data[i,]$time)$means
}
est.surv <- unlist(predictions)
f <- val.surv(S= Surv(complete.data$time, complete.data$status), 
              est.surv = est.surv)
plot(f, xlab="predicted survival", ylab="observed survival")


# drawing calibration plot 2
# calibration plot overall
prediction.weibull2 <- function(new.patient, fit.model = NULL){ 
  pct <- 10*1:300/100
  predicted.exp <- list()
  for(i in 1:n.impute){
    mu_hat <- predict(fit.model[[i]], newdata = new.patient, type = "lp")
    predicted.exp[[i]] <- 1 - pweibull(pct, shape = 1/fit.model[[i]]$scale, 
                                       scale = exp(mu_hat))
  }
  predicted.exp.average <- Reduce("+", predicted.exp) / n.impute
  
  results <- data.frame(time = pct, means = predicted.exp.average)
  return(results)
}

predictions <- list()
for (i in 1:dim(complete.data)[1]){
  predictions[[i]] <- prediction.weibull2(complete.data[i,],fit.model)$means
}

predictions.mean <- apply(do.call(cbind, predictions), 1, mean)

m2 <- survfit(Surv(time, status) ~ 1, data= complete.data) 
res <- summary(m2, censored = T)
dt1 <- with(res, data.frame(time = time, surv = surv, upper = upper,
                           lower = lower))
dt3 <- data.frame(time = c(10*1:300/100, dt1$time), 
               surv = c(predictions.mean, dt1$surv), 
               group = c(rep("pred",300), rep("obs",length(dt1$surv))))

ggplot(data = dt3) + 
  geom_line(aes(x = time, y = surv, group = group, colour = group, 
                linetype = group, size = group)) +
  scale_size_manual(values = c(0.5,0.1)) + 
  scale_color_manual(values = c("red", "black")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  xlim(0,30)

# calibration plot for imputed by levels of covariates
group1 <- complete.data[complete.data$x3 == 0,]
group2 <- complete.data[complete.data$x3 == 1,]

predictions.group1 <- list()
for (i in 1:dim(group1)[1]){
  predictions.group1[[i]] <- prediction.weibull2(group1[i,],fit.model)$means
}

predictions.group2 <- list()
for (i in 1:dim(group2)[1]){
  predictions.group2[[i]] <- prediction.weibull2(group2[i,],fit.model)$means
}

predictions.mean.group1 <- apply(do.call(cbind, predictions.group1), 1, mean)
predictions.mean.group2 <- apply(do.call(cbind, predictions.group2), 1, mean)

predictions.both.groups <- data.frame(means= c(predictions.mean.group1,predictions.mean.group2), 
                                   group= c(rep("x3=1, fitted", 300), rep("x3=2, fitted", 300)), percentiles.surv = 1:300/100)

m2 <- survfit(Surv(time, status) ~ x3, data= complete.data[complete.cases(complete.data[,1:5]),]) 
res <- summary(m2, censored = T)
dt1 <- with(res, data.frame(time = time, surv = surv, upper = upper,
                           lower = lower, strata=strata))

dt2 <- data.frame(time = c(10*1:300/100, 10*1:300/100, dt1$time), 
               surv = c(predictions.both.groups$means, dt1$surv), 
               group = c(predictions.both.groups$group, dt1$strata))
dt2$group[dt2$group == "1"] <- "x3=1, observed"
dt2$group[dt2$group == "2"] <- "x3=2, observed"

ggplot(data = dt2) + 
  geom_line(aes(x = time, y = surv, group = group, colour=group, linetype=group,
                size = group)) +
  scale_size_manual(values=c(0.8,0.1,0.8,0.1))+ 
  scale_color_manual(values=c("red", "red", "black", "black"))+
  scale_linetype_manual(values=c("solid", "dashed","solid", "dashed"))+
  xlim(0,25)


# Internal CV via bootstrapping -------------------
# bootstrapping in each multiply imputed dataset
n.bootstrap <- 10
n.impute <- 10
optimism.weibull.eachbootstrap <-  matrix(NA, n.bootstrap, 3)
optimism.weibull <- matrix(NA, n.impute, 3)

for (i in 1:n.impute){
  for(j in 1:n.bootstrap){
    
    boot.sample <- sample(length(imputed[[i]]$time), replace = T) 
    
    # create bootstrap sample
    imp.boot <- lapply(imputed[[i]], function(x){x[boot.sample]}) 
    
    # fit spline model in bootstrap sample
    fit.weibull.boot <- survreg(Surv(time, status)~pspline(x1, df=2)+pspline(x2, df=2)+ridge(x3,x4, scale=T), 
                              imp.boot, dist="weibul")
    
    # predict in bootstrap
    f1 <- as.data.frame(do.call(cbind, lapply(imp.boot, function(x) {as.numeric(as.character(x))})))
    boot.prediction.weibull <- prediction.weibull(f1, single.fit = fit.weibull.boot)
    
    weibull.boot <- calculate_performance(time = f1$time, status = f1$status,
                                          lp = boot.prediction.weibull)

    # predict in test data
    f2 <- as.data.frame(do.call(cbind, lapply(imputed[[i]], function(x) {as.numeric(as.character(x))})))
    test.prediction.weibull <- prediction.weibull(f2, single.fit = fit.weibull.boot)
    
    weibull.test <- calculate_performance(time = f2$time, status = f2$status, 
                                           lp = test.prediction.weibull)
    
    optimism.weibull.eachbootstrap[j,] <- weibull.boot - weibull.test
  }
  optimism.weibull[i,] <- apply(optimism.weibull.eachbootstrap, 2, mean)

  print(paste0("imputation done: ", i))
}

mean.optimism.weibull <- apply(optimism.weibull, 2, mean)
optimism.corrected.weibull <- apparent.weibull - mean.optimism.weibull
round(optimism.corrected.weibull,2)


# Internal-external CV -------------------
clusters <- unique(survdat$clust)
N.clust <- length(clusters) # 5 clusters in this example
data.in <- data.leftout <- list()

#create the datasets
for(i in 1:N.clust){
  data.in[[i]]<- survdat[survdat$clust!=clusters[i],]
  data.leftout[[i]]<- survdat[survdat$clust==clusters[i],]
  complete.index <- complete.cases(data.leftout[[i]][,c(paste0("x", 1:5))])
  data.leftout[[i]] <- data.leftout[[i]][complete.index,] 
}

n.impute <- 10
imputed <- fit.weibull.CV <- list()
leftout.prediction.weibull <- leftout.performance.weibull <- list()

for (i in 1:N.clust){

  data.in[[i]]$nelsonaalen <- nelsonaalen(data.in[[i]], time, status)
  meth <- make.method(data.in[[i]])
  pred <- make.predictorMatrix(data.in[[i]])
  pred[,"time"] <- 0 # don't include time variable in the imputation;instead we have baseline hazard
  imp.surv <- mice(data.in[[i]], m = n.impute)

  imputed <- list()
  impc <- complete(imp.surv, action="long")
  for(j in 1:n.impute){
    imputed[[j]] <- impc[impc$.imp==j, c(3:7, 13:15)]
  }

  for (j in 1:n.impute){
    fit.weibull.CV[[j]] <- survreg(Surv(time, status)~pspline(x1, df=2)+pspline(x2, df=2)+ridge(x3,x4, scale=T), 
                                  imputed[[j]], dist="weibul")
  }
  
  leftout.prediction.weibull[[i]] <- prediction.weibull(data.leftout[[i]], multiple.fit = fit.weibull.CV)
  
  leftout.performance.weibull[[i]] <- calculate_performance(time = data.leftout[[i]]$time, status = data.leftout[[i]]$status,
                                      lp = leftout.prediction.weibull[[i]])
}

# performance per cluster
IECV.weibull=data.frame(t(leftout.performance.weibull[[1]]))
for( i in 2:N.clust){IECV.weibull=rbind(IECV.weibull, t(leftout.performance.weibull[[i]]))}
IECV.weibull$cluster=1:N.clust
round(IECV.weibull,2)


### decision curve analysis -----
library(dcurves)
timepoint<-10
ptime <- predict(fit.model[[i]], newdata=complete.data, type='quantile', p=10, se=TRUE)

pred.complete=matrix(0, nrow=dim(complete.data)[1], ncol=n.impute)
for(i in 1:n.impute){
  mu_hat <- predict(fit.model[[i]], newdata = complete.data, type = "lp")
  pred.complete[,i]= pweibull(timepoint, shape = 1/fit.model[[i]]$scale, 
                   scale = exp(mu_hat))}


pred.complete2=rowMeans(pred.complete)

for.dca=data.frame(time=complete.data$time, status=complete.data$status, pred=pred.complete2)
dca1<-dca(Surv(time, status) ~ pred, 
          data = for.dca,thresholds = 1:80 / 100,
          time = timepoint) 
plot(dca1,smooth = TRUE)

setwd("G:/My Drive/PROJECT/practical guide to prediction models/R code/time-to-event outcome")
#save.image("G:/My Drive/PROJECT/practical guide to prediction models/R code/time-to-event outcome/t2e.RData")
load(file="t2e.RData")


