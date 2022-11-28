remove(list=ls()) # empty memory
set.seed(42) # the answer to life the universe and everything
#useful functions
logit <- function(x){log(x/(1-x))}
expit <- function(x){exp(x)/(1+exp(x))}

# sample size calculations ------------------------------------------------
library(pmsampsize)
pmsampsize(type = "b", rsquared = 0.6, parameters = 10, prevalence = 0.2)

## simulate data
library(MASS)
N <- 200

Sigma <- outer(1:10, 1:10, function(x,y) 0.5^abs(x-y)) #variance covariance matrix for covariates

x <- mvrnorm(N, rep(0,10), Sigma)
x[,3] <- ifelse(x[,3] > 0.5, 1, 0)
x[,4] <- ifelse(x[,4] > 0, 1, 0)
x[,5] <- cut(x[,5], breaks=c(-Inf, -1, 0, 1, Inf), labels = FALSE)
x[,8] <- ifelse(x[,8] > 1, 0.5, 0)
x[,9] <- ifelse(x[,9] > 1.5, 1, 0)
x[,10] <- cut(x[,10], breaks=c(-Inf, -1, 0.5, 1, Inf), labels = FALSE)

data.bin.complete <- data.frame(x)
colnames(data.bin.complete) <- c(paste0("x", 1:5), paste0("z", 1:5))

logit.py <- with(data.bin.complete,-2+x1+0.2*x1^2+
                              0.3*x2+0.1*x2^2+0.2*(x3==2)+0.2*(x4==2)+
                             0.2*(x5==2)-0.1*(x5==3)+0.2*(x5==4)+rnorm(N,0,0.1))
py <- expit(logit.py)
data.bin.complete$y <- rbinom(N,1,py)
data.bin.complete[,c(3:5, 8:10, 11)] <- lapply(data.bin.complete[,c(3:5, 8:10, 11)], factor)

# introduce missing data
missing.matrix=matrix(0, nrow=nrow(data.bin.complete), ncol=ncol(data.bin.complete))
missing.matrix=matrix(rbinom(length(missing.matrix),1, p=0.1), nrow=nrow(data.bin.complete))
data.bin=data.bin.complete
data.bin[missing.matrix==1]=NA
table(data.bin$y)

#create a clustering variable
data.bin$clust <- factor(sample(1:5, size = N, replace = TRUE, prob = rep(0.2,5)))
data.bin <- data.bin[order(data.cont$clust),]
data.bin <- data.bin[order(data.bin$clust),]
head(data.bin)

# example of a new patient
new.patient <- data.frame(x1=-0.3, x2=-0.5, x3=1, x4=1, x5=2)
new.logit <- with(new.patient,-2 +0.5*x1+0.1*x1^2+0.2*x2-0.1*x2^2+0.1*(x3==2)+0.2*(x4==2)+
                 0.2*(x5==2)-0.1*(x5==3)+0.2*(x5==4))  # this is the true log odds
expit(new.logit)

# perform multiple imputations -----------------------------------------------------
sum(complete.cases(data.bin[,1:5])) # this is the number of patients with full observations on all predictors

# Imputing missing data
library(Hmisc)
n.impute <- 10
a <- aregImpute(data=data.bin, I(y)~x1+x2+I(x3)+I(x4)+I(x5)+z1+z2+z3+z4+z5+clust, n.impute=n.impute, nk=3, match='closest')

# get imputed datasets
imputed1 <- list()
for (i in 1:n.impute){
  imputed1[[i]] <- impute.transcan(a, imputation=i, data=data.bin, list.out=TRUE,
                                  pr=FALSE, check=FALSE)
}

###### 1. logistic regression with splines --------------
library(rms)
regression.splines <- list()
for (i in 1:n.impute){  
  regression.splines[[i]]<- lrm(y~rcs(x1,3)+rcs(x2,3)+x3+x4+x5, data=imputed1[[i]]) 
}

#plot the fitted splines
spl1<-ggplot(data.frame(x1=seq(-3,3, 0.1), 
            logit.py=Predict(regression.splines[[i]], x1=seq(-3,3,0.1), x2=0, x3=1, x4=1, x5=2)$yhat), 
            aes(x=x1, y=logit.py)) + geom_smooth()

spl2<-ggplot(data.frame(x2=seq(-3,3, 0.1), 
            logit.py=Predict(regression.splines[[i]], x1=1, x2=seq(-3,3,0.1), x3=1, x4=1, x5=2)$yhat), 
            aes(x=x2, y=logit.py)) + geom_smooth()

library(gridExtra)
grid.arrange(spl1,spl2,ncol=2)


prediction.lr <- function(new.patient, single.fit = NULL, multiple.fit = NULL){ 
  
  if(!is.null(multiple.fit)){
    mygrid <- expand.grid(k = 1:dim(new.patient)[1],i = 1:length(multiple.fit))  
  
    ff <- function(k,i){ 
      with(new.patient, Predict(multiple.fit[[i]], x1= x1[k], x2 = x2[k], x3= x3[k], x4= x4[k], x5= x5[k])$y)
    }
    
    prediction_matrix <- matrix(mapply(ff, mygrid$k, mygrid$i), 
                                nrow = dim(new.patient)[1], ncol = length(multiple.fit))
    prediction <- apply(prediction_matrix, 1, mean)
  } else if(!is.null(single.fit)){
    ff <- function(k){
      with(new.patient, Predict(single.fit, x1= x1[k], x2 = x2[k], x3= x3[k], x4= x4[k], x5= x5[k])$y)
    }
    prediction <- sapply(1:dim(new.patient)[1], ff)
  }
  return(prediction)
}

expit(prediction.lr(new.patient, multiple.fit = regression.splines))
complete.data <- data.bin[complete.cases(data.bin[,c("x1","x2","x3","x4","x5","y")]),]
predicted.lr <- prediction.lr(complete.data, multiple.fit = regression.splines) # these are patients with fully observed data
ggplot(data.frame(p=expit(predicted.lr)), aes(x=p)) + 
      geom_histogram() + geom_histogram(color="black", fill="white")

###### 2. penalized regression with smoothing splines  -------------------
library(mgcv)
fit.gam <- list()
for( i in 1:n.impute){
  fit.gam[[i]] <- gam(y~x3+x4+x5+s(x1)+s(x2), data = imputed1[[i]], family=binomial)
}

prediction.gam <- function(new.patient, single.fit = NULL, multiple.fit = NULL){ 
  
  if(!is.null(multiple.fit)){
    mygrid <- expand.grid(k = 1:dim(new.patient)[1],i = 1:length(multiple.fit))
    
    ff <- function(k,i){ 
      predict.gam(multiple.fit[[i]], newdata = new.patient[k,])  
    }
    
    prediction_matrix <- matrix(mapply(ff, mygrid$k, mygrid$i), 
                                nrow = dim(new.patient)[1], ncol = length(multiple.fit))
    prediction <- apply(prediction_matrix, 1, mean)  
  } else if(!is.null(single.fit)){
    
    ff <- function(k){ 
      predict.gam(single.fit, newdata = new.patient[k,])  
    }
    prediction <- sapply(1:dim(new.patient)[1], ff)
  }
  
  return(prediction)    
}
predicted.gam <- prediction.gam(complete.data, multiple.fit = fit.gam) 
expit(prediction.gam(new.patient, multiple.fit = fit.gam))

ggplot(data.frame(p=expit(predicted.gam)), aes(x=p)) + 
      geom_histogram() + geom_histogram(color="black", fill="white")

###### 3. Ridge regression    -------------------------------
library(glmnet)
lambdas <- 10^seq(2, -10, by = -0.3) 

fit.ridge <- list()

for( i in 1:n.impute){
  imp <- imputed1[[i]] 
  imp <- with(imp, data.frame(y, x1, x2, x3, x4, x5))
  data_glmnet <- model.matrix(y ~.,data = imp)
  data_glmnet <- data_glmnet[,-1]
  data_glmnet <- cbind(y = as.numeric(as.character(imp$y)), data_glmnet = data_glmnet)
  X <- as.matrix(data_glmnet[,-1])
  colnames(X)[3:4] <- c("x3", "x4")
  Y <- data_glmnet[,1]
  cvfit <- cv.glmnet(X,Y,family = "binomial",alpha=0,
                    lambda = lambdas, nfolds=10)
  lambda.min <- cvfit$lambda.min
  fit.ridge[[i]] <- glmnet(X,Y,family = "binomial",alpha=0, lambda = lambda.min)
}

# predict for new patients
prediction.ridge <- function(new.patient, single.fit = NULL, multiple.fit = NULL){

  if(!is.null(multiple.fit)){
    mygrid <- expand.grid(k = 1:dim(new.patient)[1],i = 1:length(multiple.fit))
    
    ff <- function(k,i){ 
      
      imp <- with(new.patient, data.frame(x1[k], x2[k], x3[k], x4[k], x52= x5[k]==2,
                                          x53= x5[k]==3, x54 = x5[k]==4))
      imp[,3:7] <- lapply(imp[,3:7], as.numeric)
      colnames(imp) <- c(paste0("x",1:4), "x52", "x53", "x54")
      predict(multiple.fit[[i]], newx = as.matrix(imp))
    }  
    prediction_matrix <- matrix(mapply(ff, mygrid$k, mygrid$i), 
                                nrow = dim(new.patient)[1], ncol = length(multiple.fit))
    prediction <- apply(prediction_matrix, 1, mean)
  } else if(!is.null(single.fit)){
    ff <- function(k){ 
    
      imp <- with(new.patient, data.frame(x1[k], x2[k], x3[k], x4[k], x52= x5[k]==2,
                                          x53= x5[k]==3, x54 = x5[k]==4))
      imp[,3:7] <- lapply(imp[,3:7], as.numeric)
      colnames(imp) <- c(paste0("x",1:4), "x52", "x53", "x54")
      predict(single.fit, newx = as.matrix(imp))
    }
    prediction <- sapply(1:dim(new.patient)[1], ff)
  }
  return(prediction)
}

predicted.ridge <- prediction.ridge(complete.data, multiple.fit = fit.ridge)

expit(prediction.ridge(new.patient, multiple.fit = fit.ridge))

ggplot(data.frame(p=expit(predicted.ridge)), aes(x=p)) + geom_histogram() + 
  geom_histogram(color="black", fill="white")

#compare predictions
axislim1=c(-5,3)
p1<-ggplot(data.frame(predicted.gam=predicted.gam, predicted.lr=predicted.lr),
          aes(x=predicted.lr, y=predicted.gam))+  geom_point(size=1)+ 
          geom_abline(intercept = 0, slope = 1, color="black",linetype="dashed", size=0.5)+
          xlim(axislim1)+ylim(axislim1)

p2<-ggplot(data.frame(predicted.lr=predicted.lr, predicted.ridge=predicted.ridge),
          aes(x=predicted.ridge, y=predicted.lr))+  geom_point(size=1)+ 
          geom_abline(intercept = 0, slope = 1, color="black",linetype="dashed", size=0.5)+
          xlim(axislim1)+ylim(axislim1)

p3<-ggplot(data.frame(predicted.gam=predicted.gam, predicted.ridge=predicted.ridge),
          aes(x=predicted.ridge, y=predicted.gam))+  geom_point(size=1)+ 
          geom_abline(intercept = 0, slope = 1, color="black",linetype="dashed", size=0.5)+
          xlim(axislim1)+ylim(axislim1)
grid.arrange(p1, p2, p3, ncol=3)


# Assess apparent performance (no optimism correction) -------------------
# discrimination

#calibration in the large
mean(complete.data$y==1)
mean(expit(predicted.lr))
mean(expit(predicted.gam))
mean(expit(predicted.ridge))

calculate_performance2 <- function(observed, predicted){
  
  auc <- auc(observed~predicted)
  
  glm.fit <- summary(glm(observed~predicted, family = binomial))
  calibration.intercept <- glm.fit$coef[1,1]
  calibration.slope <- glm.fit$coef[2,1]
  vec <- c(auc, calibration.intercept, calibration.slope)
  names(vec) <- c("auc", "calibration intercept", "calibration slope")
  return(vec)
}

library(pROC)
apparent.lr <- calculate_performance2(complete.data$y, predicted.lr)
apparent.gam <- calculate_performance2(complete.data$y, predicted.gam) 
apparent.ridge <- calculate_performance2(complete.data$y, predicted.ridge)
round(rbind(apparent.lr, apparent.gam, apparent.ridge),2)

#calibration plot
val.prob(y=as.numeric(complete.data$y)-1,p=expit(predicted.lr))
val.prob(y=as.numeric(complete.data$y)-1,p=expit(predicted.gam))
val.prob(y=as.numeric(complete.data$y)-1,p=expit(predicted.ridge))

# alternate calibration slope
df.calibration <- data.frame(y = complete.data$y, predicted.lr = expit(predicted.lr))
Ngroups <- 10
d1 <- quantile(df.calibration$predicted.lr, probs = seq(0, 1, 1/Ngroups))

g1<-list()
for (i in 1:Ngroups) {
  g1[[i]] <- df.calibration[df.calibration$predicted.lr >= d1[i] & df.calibration$predicted.lr < d1[i+1],] 
}

predicted <- observed <- vector(mode = "numeric", length = Ngroups)
for (i in 1:Ngroups) {
  predicted[i] <- mean(g1[[i]]$predicted.lr)
  observed[i] <- mean(g1[[i]]$y == 1)
}
  
dat1 <- data.frame(pred = predicted, obs = observed)
ggplot(dat1,aes(x=pred,y=obs))+geom_point(size=3,shape=20)+
  labs(x="Predicted  from LR", y="Observed") +
  geom_abline(intercept=0,slope=1,color="black",linetype="dashed",size=0.7) +
  geom_smooth(method="lm",colour="blue",size=0.7) + theme(aspect.ratio=1)

# Internal CV via bootstrapping -------------------
# bootstrap in each multiply imputed dataset
N.bootstrap <- 10
optimism.lr.eachbootstrap <- optimism.gam.eachbootstrap <- optimism.ridge.eachbootstrap <- matrix(NA, N.bootstrap, 3)
optimism.lr <- optimism.gam <- optimism.ridge <- matrix(NA, n.impute, 3)

for (i in 1:n.impute){
  for(j in 1:N.bootstrap){
    
    boot.sample <- sample(length(imputed1[[i]]$y),replace = T) 
    
    # create bootstrap sample
    imp.boot <- lapply(imputed1[[i]], function(x){x[boot.sample]}) 
    
    # fit spline model in bootstrap sample
    regression.splines.boot<- lrm(y~rcs(x1,3)+rcs(x2,3)+x3+x4+x5,data=imp.boot) 
    
    # fit gam model in bootstrap sample
    fit.gam.boot <- gam(y~x3+x4+x5+s(x1)+s(x2), data = imp.boot, family = binomial)

    # fit ridge model in bootstrap sample
    imp <- imp.boot
    imp <- with(imp, data.frame(y, x1, x2, x3, x4, x5))
    data_glmnet <- model.matrix(y ~.,data = imp)
    data_glmnet <- data_glmnet[,-1]
    data_glmnet <- cbind(y = as.numeric(as.character(imp$y)), data_glmnet = data_glmnet)
    X <- as.matrix(data_glmnet[,-1])
    colnames(X)[3:4] <- c("x3", "x4")
    Y <- data_glmnet[,1]
    cvfit <- cv.glmnet(X,Y,family = "binomial",alpha=0,
                       lambda = lambdas, nfolds=10)
    lambda.min <- cvfit$lambda.min
    fit.ridge.boot <- glmnet(X,Y,family = "binomial", alpha=0, lambda = lambda.min)
    
    # predict in bootstrap
    f1 <- as.data.frame(do.call(cbind, lapply(imp.boot, function(x) {as.numeric(as.character(x))})))
    
    boot.prediction.lr <- prediction.lr(f1, single.fit = regression.splines.boot)
    boot.prediction.gam <- prediction.gam(f1, single.fit = fit.gam.boot)
    boot.prediction.ridge <- prediction.ridge(f1, single.fit = fit.ridge.boot)
    
    lr.boot <- calculate_performance2(f1$y, boot.prediction.lr)
    gam.boot <- calculate_performance2(f1$y, boot.prediction.gam)
    ridge.boot <- calculate_performance2(f1$y, boot.prediction.ridge)
    
    # predict in test data
    f2 <- as.data.frame(do.call(cbind, lapply(imputed1[[i]], function(x) {as.numeric(as.character(x))})))
    test.prediction.lr <- prediction.lr(f2, single.fit = regression.splines.boot)
    test.prediction.gam <- prediction.gam(f2, single.fit = fit.gam.boot)
    test.prediction.ridge <- prediction.ridge(f2, single.fit = fit.ridge.boot)
    
    lr.test <- calculate_performance2(f2$y, test.prediction.lr)
    gam.test <- calculate_performance2(f2$y, test.prediction.gam)
    ridge.test <- calculate_performance2(f2$y, test.prediction.ridge)
    
    optimism.lr.eachbootstrap[j,] <- lr.boot - lr.test
    optimism.gam.eachbootstrap[j,] <- gam.boot - gam.test
    optimism.ridge.eachbootstrap[j,] <- ridge.boot - ridge.test
  }
  optimism.lr[i,] <- apply(optimism.lr.eachbootstrap, 2, mean)
  optimism.gam[i,] <- apply(optimism.gam.eachbootstrap, 2, mean)
  optimism.ridge[i,] <- apply(optimism.ridge.eachbootstrap, 2, mean)
  
  print(paste0("imputation done: ", i))
}
    
mean.optimism.lr <- apply(optimism.lr, 2, mean)
mean.optimism.gam <- apply(optimism.gam, 2, mean)
mean.optimism.ridge <- apply(optimism.ridge, 2, mean)

optimism.corrected.lr <- apparent.lr - mean.optimism.lr
optimism.corrected.gam <- apparent.gam - mean.optimism.gam
optimism.corrected.ridge <- apparent.ridge - mean.optimism.ridge
round(rbind(optimism.corrected.lr, optimism.corrected.gam, optimism.corrected.ridge),2)


# Internal-external CV -------------------
clusters <- unique(data.bin$clust)
N.clust <- length(clusters) # 5 clusters in this example
data.in <- data.leftout <- list()

#create the datasets
for(i in 1:N.clust){
  data.in[[i]]<-data.bin[data.bin$clust!=clusters[i],]
  data.leftout[[i]]<-data.bin[data.bin$clust==clusters[i],]
  complete.index <- complete.cases(data.leftout[[i]][,c(paste0("x", 1:5), "y")])
  data.leftout[[i]] <- data.leftout[[i]][complete.index,] 
}

n.impute <- 10

imputed <- regression.splines.CV <- fit.gam.CV <- fit.ridge.CV <- list()
leftout.prediction.lr <- leftout.prediction.gam <- leftout.prediction.ridge <- list()
leftout.performance.lr <- leftout.performance.gam <- leftout.performance.ridge <- list()

for (i in 1:N.clust){
  a <- aregImpute(data=data.in[[i]], I(y)~x1+x2+I(x3)+I(x4)+I(x5), n.impute=n.impute, nk=3, match='closest')
  
  for (j in 1:n.impute){
    imputed[[j]] <- impute.transcan(a, imputation=j, data=data.in[[i]], list.out=TRUE,
                                    pr=FALSE, check=FALSE)
    regression.splines.CV[[j]]<- lrm(y~rcs(x1,3)+rcs(x2,3)+x3+x4+x5,data=imputed[[j]]) 
    fit.gam.CV[[j]] <- gam(y ~ x3+x4+x5+s(x1)+s(x2), data = imputed[[j]], family = binomial)
    
    imp <- with(imputed[[j]], data.frame(y, x1, x2, x3, x4, x5))
    data_glmnet <- model.matrix(y ~.,data = imp)
    data_glmnet <- data_glmnet[,-1]
    data_glmnet <- cbind(y = as.numeric(as.character(imp$y)), data_glmnet = data_glmnet)
    X <- as.matrix(data_glmnet[,-1])
    colnames(X)[3:4] <- c("x3", "x4")
    Y <- data_glmnet[,1]
    cvfit <- cv.glmnet(X,Y,family = "binomial",alpha=0,
                       lambda = lambdas, nfolds=10)
    lambda.min <- cvfit$lambda.min
    fit.ridge.CV[[j]] <- glmnet(X,Y,family = "binomial", alpha=0, lambda = lambda.min)
  }
  
  leftout.prediction.lr[[i]] <- prediction.lr(data.leftout[[i]], multiple.fit = regression.splines.CV)
  leftout.prediction.gam[[i]] <- prediction.gam(data.leftout[[i]], multiple.fit = fit.gam.CV)
  leftout.prediction.ridge[[i]] <- prediction.ridge(data.leftout[[i]], multiple.fit = fit.ridge.CV)
  
  leftout.performance.lr[[i]] <- calculate_performance2(data.leftout[[i]]$y, leftout.prediction.lr[[i]])
  leftout.performance.gam[[i]] <- calculate_performance2(data.leftout[[i]]$y, leftout.prediction.gam[[i]])
  leftout.performance.ridge[[i]] <- calculate_performance2(data.leftout[[i]]$y, leftout.prediction.ridge[[i]])
}

#performance per cluster
#leftout.performance.lr
#leftout.performance.gam
#leftout.performance.ridge

# performance aggregating all cluster
leftout.prediction.lr.agg <- do.call(c, leftout.prediction.lr)
leftout.prediction.gam.agg <- do.call(c, leftout.prediction.gam)
leftout.prediction.ridge.agg <- do.call(c, leftout.prediction.ridge)
IECV.observed <- do.call(rbind, data.leftout)$y
IECV.cluster <- do.call(rbind, data.leftout)$clust

IECV.lr <- calculate_performance2(IECV.observed, leftout.prediction.lr.agg)
IECV.gam <- calculate_performance2(IECV.observed, leftout.prediction.gam.agg)
IECV.ridge <- calculate_performance2(IECV.observed, leftout.prediction.ridge.agg)
round(rbind(IECV.lr, IECV.gam, IECV.ridge),2)

#CI for ROC curve
#pROC::ci(roc(IECV.observed,leftout.prediction.lr.agg))
#pROC::ci(roc(IECV.observed,leftout.prediction.gam.agg))
#pROC::ci(roc(IECV.observed,leftout.prediction.ridge.agg))


# random effects meta-analysis of AUC
auc.clusters <- data.frame(auc=rep(NA,N.clust), SE=NA, cluster=NA)
for(i in 1:N.clust){
  d.cl <- data.leftout[[i]]
  roc1 <- roc(d.cl$y, leftout.prediction.ridge[[i]])
  auc.clusters$auc[i] <-auc(roc1)
  auc.clusters$SE[i] <-(pROC::ci(roc1)[3]-pROC::ci(roc1)[1])/3.92
  auc.clusters$cluster[i]<-clusters[i]
}
round(auc.clusters,digits=2)

#forest plot
library(meta)
meta.AUC<-metagen(TE=auc, seTE=SE, studlab = cluster, data=auc.clusters)
forestplot<-forest(meta.AUC, prediction = F, xlim=c(0.4,1),
           colgap.left="5mm", rightcols = c("effect", "ci"),
           leftlabs = c("Cluster", "AUC", "seTE"))


# plot a decision curve analysis
# see also https://cran.r-project.org/web/packages/dcurves/vignettes/dca.html
library(dcurves)
for.dca=data.frame(obs=as.numeric(complete.data$y)-1, 
                   predicted.gam=expit(predicted.gam),
                   predicted.lr=expit(predicted.lr), predicted.ridge=expit(predicted.ridge))
for.dca=data.frame(obs=as.numeric(complete.data$y)-1, 
                   predicted.gam=expit(predicted.gam), predicted.ridge=expit(predicted.ridge),
                   predicted.lr=expit(predicted.lr))
dca1<-dca(obs~predicted.lr+predicted.ridge+predicted.gam, data=for.dca)
plot(dca1, smooth = TRUE)


setwd("G:/My Drive/PROJECT/practical guide to prediction models/R code/binary outcome")
#save.image(file="binary_example.RData")
load(file="binary_example.RData")
