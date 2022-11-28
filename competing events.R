# parts of the code copied from https://github.com/survival-lumc/ValidationCompRisks/blob/main/Prediction_CSC_minimal.R
# Original paper https://www.bmj.com/content/bmj/377/bmj-2021-069249.full.pdf

remove(list=ls())
set.seed(42) # the answer to life the universe and everything

## simulate data ------------
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
rate1 <- with(survdat.compl, exp(-0.5+x1+0.4*x2+(x3==2)+0.5*(x4==2)+0.5*(x5==1)+rnorm(N,0,0.1))/10)
rate2 <- with(survdat.compl, exp(-1+0.1*x1+0.2*x2+(x3==2)+0.1*(x4==2)-0.2*(x5==1)+rnorm(N,0,0.1))/10)

survdat.compl[,c(3:5, 8:10)] <- lapply(survdat.compl[,c(3:5, 8:10)], factor)

fulltime <- rexp(N, rate = rate1)
fulltime2 <- rexp(N, rate = rate2)
censtimes <- 10 + 20*runif(N)

mintime <- pmin(fulltime, fulltime2, censtimes)
survdat.compl$time <- mintime

survdat.compl$status <- 0
survdat.compl$status[mintime == fulltime] <- 1
survdat.compl$status[mintime == fulltime2] <- 2

# introduce missing data for the covariates
missing.matrix=matrix(0, nrow=nrow(survdat.compl), ncol=10)
missing.matrix=matrix(rbinom(length(missing.matrix),1, p=0.02), nrow=nrow(survdat.compl))
survdat=survdat.compl[,1:10]
survdat[missing.matrix==1]=NA
survdat=cbind(survdat, survdat.compl[11:12])

#create clusters
survdat$clust <- factor(sample(1:5, size = N, replace = TRUE, prob = rep(0.2,5)))
survdat <- survdat[order(survdat$clust),]

#inspect data 
head(survdat)
table(survdat$status)

# Impute missing data --------------
library(mice)
n.impute <- 10
survdat$logtime <- log(survdat$time)
meth <- make.method(survdat)
pred <- make.predictorMatrix(survdat)
pred[,c("time")] <- 0
imp.surv <- mice(survdat, pred = pred, meth = meth, m = n.impute)

# get imputed datasets
imputed <- list()
impc <- complete(imp.surv, action="long")
for(i in 1:n.impute){
  imputed[[i]] <- impc[impc$.imp==i, c(3:7, 13:15)]
}


# visualize the data -------
library(cmprsk)
ci_fit <- 
  cuminc(ftime = survdat.compl$time,fstatus=survdat.compl$status,cencode=0)
library(survminer)
ggcompetingrisks(ci_fit, xlab = "Days", conf.int = T, multiple_panels = FALSE)

# Fit cause-specific hazards models  ------ 
library(riskRegression)
fit.model <- list()
for(i in 1:n.impute){
  fit.model[[i]] <- CSC(
    formula = Hist(time, status) ~ x1 + x2 + x3 + x4 + x5,
    data = imputed[[i]]
  )
}

# prediction for new patient
prediction.competing <- function(new.patient, single.fit = NULL, 
                                 multiple.fit = NULL, 
                                 time.horizon = NULL, primary.event = NULL){ 
  
  if(!is.null(multiple.fit)){
    
    ff <- function(i){
      predictRisk(
        object = multiple.fit[[i]], 
        cause = primary.event, 
        newdata = new.patient, 
        times = time.horizon
      )  
    }
    prediction_matrix <- sapply(1:length(multiple.fit), ff)
    
    if(dim(new.patient)[1] == 1){
      prediction <- mean(prediction_matrix)
    } else{
      prediction <- apply(prediction_matrix, 1, mean)
    }
    
  } else if(!is.null(single.fit)){
    prediction <- 
      predictRisk(
        object = single.fit, 
        cause = primary.event, 
        newdata = new.patient, 
        times = time.horizon
      )   
  }
  return(prediction)
}

time.horizon <- 10
new.patient <- data.frame(x1 = 1, x2 = -1.2, x3 = as.factor(1), x4 = as.factor(0), x5=as.factor(0))
prediction.competing(new.patient, multiple.fit = fit.model, time.horizon = time.horizon, primary.event = 1)
prediction.competing(new.patient, multiple.fit = fit.model, time.horizon = time.horizon, primary.event = 2)

# plot survival curves
n.points=30
new.patient.prediction=data.frame(time=seq(0,30,length.out=n.points))

for (i in 1:n.points){
  new.patient.prediction$new.pred1[i]<-prediction.competing(new.patient, multiple.fit = fit.model, time.horizon = new.patient.prediction$time[i], primary.event = 1) 
  new.patient.prediction$new.pred2[i]<-prediction.competing(new.patient, multiple.fit = fit.model, time.horizon = new.patient.prediction$time[i], primary.event = 2) 
  
}
predictions.for.plot=data.frame(time=
  rep(new.patient.prediction$time,2), 
                                pred=1-c(new.patient.prediction$new.pred1, 
                                  new.patient.prediction$new.pred2 ), 
                                event=c(rep("1",n.points), rep("2",n.points)))
                                
ggplot(predictions.for.plot, aes(x = time, y = pred, color = event)) +
geom_line(aes(group=factor(event)),size=1)

#predict for the complete dataset
complete.data <- survdat[complete.cases(survdat[,1:5]),]
predicted.competing.1 <- prediction.competing(complete.data, multiple.fit = fit.model, time.horizon = time.horizon, primary.event = 1)
predicted.competing.2 <- prediction.competing(complete.data, multiple.fit = fit.model, time.horizon = time.horizon, primary.event = 2)


#  Assess apparent performance (no optimism correction) -------------------
library(survival)
time.horizon<-10
# event 1
primary.event=1
complete.data$event <- factor(complete.data$status, 0:2, labels=c("censor", "event", "competing"))
  
obj <- summary(survfit(Surv(time, event) ~ 1, data = complete.data), times = time.horizon)
aj <- list("obs" = obj$pstate[, primary.event + 1], "se" = obj$std.err[, primary.event + 1])

# Calculate O/E
OE <- aj$obs / mean(predicted.competing.1)

# For the confidence interval we use method proposed in Debray et al. (2017) doi:10.1136/bmj.i6460
OE_summary1 <- c(
  "OE" = OE,
  "lower" = exp(log(OE - qnorm(0.975) * aj$se / aj$obs)),
  "upper" = exp(log(OE + qnorm(0.975) * aj$se / aj$obs))
)
round(OE_summary1,2)

# event 2
primary.event=2
complete.data$event <- factor(complete.data$status, 0:2, labels=c("censor", "event", "competing"))

obj <- summary(survfit(Surv(time, event) ~ 1, data = complete.data), times = time.horizon)
aj <- list("obs" = obj$pstate[, primary.event + 1], "se" = obj$std.err[, primary.event + 1])

# Calculate O/E
OE <- aj$obs / mean(predicted.competing.2)

# For the confidence interval we use method proposed in Debray et al. (2017) doi:10.1136/bmj.i6460
OE_summary2 <- c(
  "OE" = OE,
  "lower" = exp(log(OE - qnorm(0.975) * aj$se / aj$obs)),
  "upper" = exp(log(OE + qnorm(0.975) * aj$se / aj$obs))
)
round(OE_summary2,2)


library(geepack) 
library(pec) 

calculate_performance <- function(fit.model = NULL, time = NULL, 
                                  status = NULL, data.used = NULL,
                                  time.horizon = NULL, primary.event = NULL){
  
  score_vdata <- Score(
    list("csh_validation" = fit.model),
    formula = Hist(time, status) ~ 1, 
    cens.model = "km", 
    data = data.used, 
    conf.int = TRUE, 
    times = time.horizon,
    metrics = c("auc", "brier"),
    summary = c("ipa"), 
    cause = primary.event,
    plots = "calibration" )

  pseudos <- data.frame(score_vdata$Calibration$plotframe)
  pseudos <- pseudos[order(pseudos$risk), ]
  pseudos$cll_pred <- log(-log(1 - pseudos$risk))  
  
  # Fit model for calibration intercept
  fit_cal_int <- geese(
    pseudovalue ~ offset(cll_pred), 
    data = pseudos,
    id = ID, 
    scale.fix = TRUE, 
    family = gaussian,
    mean.link = "cloglog",
    corstr = "independence", 
    jack = TRUE  )
  
  # Fit model for calibration slope
  fit_cal_slope <- geese(
    pseudovalue ~ offset(cll_pred) + cll_pred,
    data = pseudos,
    id = ID, 
    scale.fix = TRUE, 
    family = gaussian,
    mean.link = "cloglog",
    corstr = "independence", 
    jack = TRUE  )

  AUC <- score_vdata$AUC$score$AUC
  
  cindex_csh <- cindex(
    object = fit.model, 
    formula = Hist(time, status) ~ 1, 
    cause = primary.event, 
    eval.times = time.horizon, 
    data = data.used
  )$AppCindex$CauseSpecificCox
  
  returnVec <- c("calibration intercept" = summary(fit_cal_int)$mean$estimate,
                 "calibration slope" = 1 + summary(fit_cal_slope)$mean["cll_pred",]$estimate,
                 "AUC" = AUC,         
                 "c-index" = cindex_csh)
  
  return(returnVec)
}  

apparent.competing.list <- matrix(NA, nrow = n.impute, ncol = 4)
for(i in 1:n.impute){
  apparent.competing.list[i,]<- calculate_performance(fit.model = fit.model[[i]],
                                                      time = complete.data$time,
                                                      status = complete.data$status,
                                                      data.used = complete.data,
                                                      time.horizon = time.horizon,
                                                      primary.event = 1)}
apparent.competing1 <- apply(apparent.competing.list, 2, mean)
names(apparent.competing1)<-c("calibration intercept", "calibration slope", "AUC", "c-index")
round(apparent.competing1,2)


#  calibration plot
primary.event<-1
pseudos=list()
smooth_pseudos=list()
for(i in 1:n.impute){
  score_vdata <- Score(
    list("csh_validation" = fit.model[[i]]),
    formula = Hist(time, status) ~ 1, 
    cens.model = "km", 
    data = complete.data, 
    conf.int = TRUE, 
    times = time.horizon,
    metrics = c("auc", "brier"),
    summary = c("ipa"), 
    cause = primary.event,
    plots = "calibration"
  ) 
  pseudos[[i]] <- data.frame( score_vdata$Calibration$plotframe)
  pseudos[[i]] <- pseudos[[i]][order(pseudos[[i]]$risk), ]
  smooth_pseudos[[i]] <- predict(
    stats::loess(pseudovalue ~ risk, data = pseudos[[i]], degree = 1, span = 0.33), 
    se = TRUE)
  
}

ps.risk=pseudos[[1]]$risk/n.impute; for(i in 2:n.impute){ps.risk=ps.risk+pseudos[[i]]$risk/n.impute}
ps.val=pseudos[[1]]$pseudovalue/n.impute; for(i in 2:n.impute){ps.val=ps.val+pseudos[[i]]$pseudovalue/n.impute}
ps.fit=smooth_pseudos[[1]]$fit/n.impute; for(i in 2:n.impute){ps.fit=ps.fit+smooth_pseudos[[i]]$fit/n.impute}
ps.df=smooth_pseudos[[1]]$df/n.impute; for(i in 2:n.impute){ps.df=ps.df+smooth_pseudos[[i]]$df/n.impute}
ps.se=smooth_pseudos[[1]]$se/n.impute;for(i in 2:n.impute){ps.se=ps.se+smooth_pseudos[[i]]$se/n.impute}

spike_bounds <- c(-0.075, 0)
bin_breaks <- seq(0, 0.6, length.out = 100 + 1)
freqs <- table(cut(predicted.competing.1, breaks = bin_breaks))
bins <- bin_breaks[-1]
freqs_valid <- freqs[freqs > 0]
freqs_rescaled <- spike_bounds[1] + (spike_bounds[2] - spike_bounds[1]) * 
  (freqs_valid - min(freqs_valid)) / (max(freqs_valid) - min(freqs_valid))


# produce plot
plot(
  x = ps.risk, 
  y = ps.val,
  xlim = c(0, 0.6), 
  ylim = c(spike_bounds[1], 0.6),
  yaxt = "n",
  frame.plot = FALSE,
  xlab = "Estimated risks",
  ylab = "Observed outcome proportions", 
  type = "n"
)
axis(2, seq(0, 0.6, by = 0.1), labels = seq(0, 0.6, by = 0.1))
polygon(
  x = c(ps.risk, rev(ps.risk)),
  y = c(
    pmax(ps.fit - qt(0.975, ps.df) * ps.se, 0),
    rev(ps.fit + qt(0.975, ps.df) * ps.se)
  ),
  border = FALSE,
  col = "lightgray"
)
abline(a = 0, b = 1, col = "gray")
lines(x = ps.risk, y = ps.fit, lwd = 2)
segments(
  x0 = bins[freqs > 0], 
  y0 = spike_bounds[1], 
  x1 = bins[freqs > 0], 
  y1 = freqs_rescaled
)



# Internal CV via bootstrapping -------------------
# bootstrap in each multiply imputed dataset
n.bootstrap <- 10
n.impute <- 10
primary.event<-1
time.horizon<-10
optimism.competing.eachbootstrap <-  matrix(NA, n.bootstrap, 4)
optimism.competing <- matrix(NA, n.impute, 4)

for (i in 1:n.impute){
  for(j in 1:n.bootstrap){
    
    boot.sample <- sample(length(imputed[[i]]$time), replace = T) 
    
    # create bootstrap sample
    imp.boot <- as.data.frame(lapply(imputed[[i]], function(x){x[boot.sample]})) 
    
    # predict in bootstrap
    f1 <- as.data.frame(do.call(cbind, lapply(imp.boot, function(x) {as.numeric(as.character(x))})))
    
    # fit cause specific hazard model
    fit.competing.boot <- CSC(formula = Hist(time, status) ~ x1 + x2 + x3 + x4,
                              data = f1)
    
    competing.boot <- calculate_performance(fit.model = fit.competing.boot,
                                            time = f1$time, 
                                            status = f1$status,
                                            data.used = f1,
                                            time.horizon = time.horizon,
                                            primary.event = primary.event)
    
    # predict in test data
    f2 <- as.data.frame(do.call(cbind, lapply(imputed[[i]], function(x) {as.numeric(as.character(x))})))
 
    competing.test <- calculate_performance(fit.model = fit.competing.boot,
                                            time = f2$time, 
                                            status = f2$status,
                                            data.used = f2,
                                            time.horizon = time.horizon,
                                            primary.event = primary.event)
    
    optimism.competing.eachbootstrap[j,] <- competing.boot - competing.test
  }
  optimism.competing[i,] <- apply(optimism.competing.eachbootstrap, 2, mean)
  
  print(paste0("imputation done: ", i))
}

mean.optimism.competing <- apply(optimism.competing, 2, mean)
optimism.corrected.competing <- apparent.competing1 - mean.optimism.competing
round(apparent.competing1,3)
round(optimism.corrected.competing,3)



# Internal-external CV -------------------
clusters <- unique(survdat$clust)
N.clust <- length(clusters) # 10 clusters in this example
data.in <- data.leftout <- list()

#create the datasets
for(i in 1:N.clust){
  data.in[[i]]<- survdat[survdat$clust!=clusters[i],]
  data.leftout[[i]]<- survdat[survdat$clust==clusters[i],]
  complete.index <- complete.cases(data.leftout[[i]][,c(paste0("x", 1:5))])
  data.leftout[[i]] <- data.leftout[[i]][complete.index,] 
}

n.impute <- 10
imputed <- fit.competing.CV <- list()
leftout.prediction.competing <- leftout.performance.competing <- list()

for (i in 1:N.clust){
  
  data.in[[i]]$logtime <- log(data.in[[i]]$time)
  meth <- make.method(data.in[[i]])
  pred <- make.predictorMatrix(data.in[[i]])
  pred[,c("time")] <- 0
  imp.surv <- mice(data.in[[i]], pred = pred, meth = meth, m = n.impute)
  
  # get imputed datasets
  imputed <- list()
  impc <- complete(imp.surv, action="long")
  for(j in 1:n.impute){
    imputed[[j]] <- impc[impc$.imp==j, c(3:7, 13:15)]
  }
  
  leftout.performance.competing.each <- matrix(NA, nrow = n.impute, ncol = 4)
  for(j in 1:n.impute){
    fit.competing.CV <- CSC(
      formula = Hist(time, status) ~ x1 + x2 + x3 + x4,
      data = imputed[[j]]
    )
    
    leftout.performance.competing.each[j,] <- calculate_performance(
                          fit.model = fit.competing.CV,
                          time = data.leftout$time, 
                          status = data.leftout$status,
                          data.used = imputed[[j]],
                          time.horizon = time.horizon,
                          primary.event = primary.event)
  }
  leftout.performance.competing[[i]] <- apply(leftout.performance.competing.each, 2, mean)
}

# performance per cluster
leftout.performance.competing
int.ext.competing<-data.frame(cluster=1:N.clust,round(matrix(unlist(leftout.performance.competing),nrow=N.clust, byrow = T),2))
colnames(int.ext.competing)=c("cluster","calibration intercept", "calibration slope", "AUC", "c-index")
int.ext.competing



### calculate CIs for apparent performance ------- 
library(geepack) #geese
library(pec) #c-index

calculate_performance <- function(fitmodel = NULL, time = NULL, 
                                  status = NULL, data.used = NULL,
                                  time.horizon = NULL, primary.event = NULL, 
                                  bootstrap=100){
  
  score_vdata <- Score(
    list("csh_validation" = fitmodel),
    formula = Hist(time, status) ~ 1, 
    cens.model = "km", 
    data = data.used, 
    conf.int = TRUE, 
    times = time.horizon,
    metrics = c("auc", "brier"),
    summary = c("ipa"), 
    cause = primary.event,
    plots = "calibration" )
  
  pseudos <- data.frame(score_vdata$Calibration$plotframe)
  pseudos <- pseudos[order(pseudos$risk), ]
  pseudos$cll_pred <- log(-log(1 - pseudos$risk))  
  
  # Fit model for calibration intercept
  fit_cal_int <- geese(
    pseudovalue ~ offset(cll_pred), 
    data = pseudos,
    id = ID, 
    scale.fix = TRUE, 
    family = gaussian,
    mean.link = "cloglog",
    corstr = "independence", 
    jack = TRUE  )
  
  # Fit model for calibration slope
  fit_cal_slope <- geese(
    pseudovalue ~ offset(cll_pred) + cll_pred,
    data = pseudos,
    id = ID, 
    scale.fix = TRUE, 
    family = gaussian,
    mean.link = "cloglog",
    corstr = "independence", 
    jack = TRUE  )
  
  AUC <- score_vdata$AUC$score$AUC
  se.AUC <- score_vdata$AUC$score$se
  
  ### CI bootstraps
  boots_ls <- lapply(seq_len(bootstrap), function(b) {
    vdata_boot <- data.used[sample(nrow(data.used), replace = TRUE), ]
    cindex_boot <- pec::cindex(
      object = fitmodel, 
      formula = Hist(time, status) ~ 1,
      cause = primary.event, 
      eval.times = time.horizon, 
      data = vdata_boot,
      verbose = FALSE
    )$AppCindex$CauseSpecificCox
    
    
    cbind.data.frame("cindex" = cindex_boot)
  })
  
  cindex <- do.call(rbind.data.frame, boots_ls)
  
  #####
  returnVec <- c("cal.int" = summary(fit_cal_int)$mean$estimate,
                 "cal.int.se"=summary(fit_cal_int)$mean$san.se,
                 "cal.slope" = 1 + summary(fit_cal_slope)$mean["cll_pred",]$estimate,
                 "cal.slope.se"=summary(fit_cal_slope)$mean["cll_pred",]$san.se,
                 "AUC" = AUC,"se.AUC"=se.AUC,         
                 "c.index" = cindex)
  
  return(returnVec)
}  

apparent.competing.list <- list()
cindex=c(); auc=c(); se.auc=c(); cal.int=c();cal.int.se=c();cal.slope=c();cal.slope.se=c() 

for(i in 1:n.impute){
  apparent.competing.list[[i]]<- calculate_performance(fitmodel = fit.model[[i]],
                                                       time = complete.data$time,
                                                       status = complete.data$status,
                                                       data.used = complete.data,
                                                       time.horizon = time.horizon,
                                                       primary.event = 1, 
                                                       bootstrap=10)
  cindex=c(cindex, apparent.competing.list[[i]]$c.index)
  auc=c(auc,apparent.competing.list[[i]]$AUC )
  se.auc=c(se.auc, apparent.competing.list[[i]]$se.AUC)
  cal.int=c(cal.int,apparent.competing.list[[i]]$cal.int)
  cal.int.se=c(cal.int.se,apparent.competing.list[[i]]$cal.int.se)
  cal.slope=c(cal.slope,apparent.competing.list[[i]]$cal.slope)
  cal.slope.se=c(cal.slope.se,apparent.competing.list[[i]]$cal.slope.se)
}

apparent.cindex=round(c(mean(cindex), quantile(cindex, probs = c(0.025, 0.975)) ),2)
apparent.se.auc=sqrt(mean(se.auc^2)+(1+1/n.impute)/n.impute*(sum((mean(auc)-auc)^2)))
apparent.AUC=round(c(mean(auc), mean(auc)-1.96*apparent.se.auc, mean(auc)+1.96*apparent.se.auc),2)
names(apparent.AUC)=c("mean", "2.5%", "97.5%")

apparent.se.int=sqrt(mean(cal.int.se^2)+(1+1/n.impute)/n.impute*(sum((mean(cal.int)-cal.int)^2)))
apparent.se.slope=sqrt(mean(cal.slope.se^2)+(1+1/n.impute)/n.impute*(sum((mean(cal.slope)-cal.slope)^2)))

apparent.int=round(c(mean(cal.int), mean(cal.int)-1.96*apparent.se.int, mean(cal.int)+1.96*apparent.se.int),2)
names(apparent.int)=c("mean", "2.5%", "97.5%")

apparent.slope=round(c(mean(cal.slope), mean(cal.slope)-1.96*apparent.se.slope, mean(cal.slope)+1.96*apparent.se.slope),2)
names(apparent.slope)=c("mean", "2.5%", "97.5%")

# results
apparent.cindex
apparent.AUC
apparent.int
apparent.slope


### Decision curve analysis
primary.event<-1
time.of.interest<-15
predicted.for.dca1 <- prediction.competing(complete.data, multiple.fit = fit.model,
                                           time.horizon = time.of.interest, primary.event = 1)


for.dca=data.frame(time=complete.data$time, status=complete.data$status, pred=predicted.for.dca1)
for.dca$status[for.dca$status!=primary.event]=0
dca1<-dca(Surv(time, status) ~ pred, 
          data = for.dca,thresholds = 1:100 / 100,
          time = time.horizon) 
plot(dca1,smooth = TRUE)

setwd("G:/My Drive/PROJECT/practical guide to prediction models/R code/competing events")
#save.image(file="competing_example.RData")
load(file="competing_example.RData")
