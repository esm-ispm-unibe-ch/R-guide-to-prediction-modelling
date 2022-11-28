remove(list=ls()) # empty memory
set.seed(42) # the answer to life the universe and everything

# sample size calculations ------------------------------------------------
library(pmsampsize)
pmsampsize(type = "c", rsquared = 0.6, parameters = 10, intercept = 0, sd = 1)

# simulating a toy example:
# generate dataset --------------------------------------------------
# The dataset will include 
#  - 5 predictors x1, x2, ... x5 (x1,x2: continuous, x3,x4: binary, x5 categorical)
#  - the outcome y we want to predict
#  - 5 auxilliary variables z1:z5. These will only be used in imputation
#  - a clustering variable clust: this will be used for internal-external CV
library(MASS) 
N <- 100
Sigma <- outer(1:10, 1:10, function(x,y) 0.5^abs(x-y)) #variance covariance matrix for covariates

x <- mvrnorm(N, rep(0,10), Sigma)
x[,3]<- ifelse(x[,3] > 0.5, 1, 0)
x[,4] <- ifelse(x[,4] > 0, 1, 0)
x[,5] <- cut(x[,5], breaks=c(-Inf, -1, 0, 1, Inf), labels = FALSE)
x[,8] <- ifelse(x[,8] > -0.5, 1, 0)
x[,9] <- ifelse(x[,9] > 0.5, 1, 0)
x[,10] <- cut(x[,10], breaks=c(-Inf, -1, 0, 1, Inf), labels = FALSE)

data.cont.complete <- data.frame(x)
colnames(data.cont.complete) <- c(paste0("x", 1:5), paste0("z", 1:5))

data.cont.complete$y<-with(data.cont.complete, x1+0.2*x1^2+0.5*x2-0.2*x2^2+0.3*x3+0.2*x4+
                             0.2*(x5==2)-0.1*(x5==3)+0.4*(x5==4)+rnorm(N,0,1))
data.cont.complete[,c(3:5, 8:10)] <- lapply(data.cont.complete[,c(3:5, 8:10)], factor)
head(data.cont.complete)

# introduce missing data
missing.matrix=matrix(0, nrow=nrow(data.cont.complete), ncol=ncol(data.cont.complete))
missing.matrix=matrix(rbinom(length(missing.matrix),1, p=0.1), nrow=nrow(data.cont.complete))
data.cont=data.cont.complete
data.cont[missing.matrix==1]=NA

#create a clustering variable
data.cont$clust <- factor(sample(1:5, size = N, replace = TRUE, prob = rep(0.2,5)))
data.cont <- data.cont[order(data.cont$clust),]
head(data.cont)

library(ggplot2)
ggplot(data.cont, aes(x=y)) + geom_histogram(color="black", fill="white", binwidth = 0.5)
summary(data.cont$y)

# example of a new patient
new.patient<-data.frame(x1=1.2, x2=-1.6, x3=0, x4=1, x5=4)

# this is the true expected outcome
with(new.patient,x1+0.2*x1^2-0.5*x2-0.2*x2^2+0.3*x3+0.2*x4+
       0.2*(x5==2)-0.1*(x5==3)+0.4*(x5==4))  

# perform multiple imputations -----------------------------------------------------
sum(complete.cases(data.cont[,1:5])) # this is the number of patients with full observations on y, x, z
sum(complete.cases(data.cont$y)) # this is the number of patients with full observations on y

library(mice)
md.pattern(data.cont, plot = FALSE)

 library(VIM) #for visualizations of missing data
 aggr_plot <- aggr(data.cont, col=c('navyblue','red'),
                  numbers=TRUE, sortVars=TRUE,
                  labels=names(data.cont),
                  cex.axis=.7, gap=3,
                  ylab=c("Histogram of missing data","Pattern"))


library(Hmisc)
n.impute <- 10
a <- aregImpute(data=data.cont, I(y)~x1+x2+I(x3)+I(x4)+I(x5)+
                z1+z2+I(z3)+I(z4)+I(z5)+clust, n.impute=n.impute, nk=3, match='closest')

# get imputed datasets
imputed1 <- list()
for (i in 1:n.impute){
  imputed1[[i]] <- impute.transcan(a, imputation=i, data=data.cont, list.out=TRUE, pr=FALSE, check=FALSE)
}

# fit prediction models ---------------------------------------------------
# You may want to first remove patients with imputed outcomes - check the paper!
# In that case first run:
# missing.y <- which(is.na(data.cont$y))
# for (i in 1:n.impute){
#   imputed1[[i]]<-as.data.frame(imputed1[[i]])
#   imputed1[[i]]<-imputed1[[i]][-missing.y,]
# }

###### 1. linear regression with splines -------
library(rms)
regression.splines<-list()
for (i in 1:n.impute){
  regression.splines[[i]]<- ols(y~rcs(x1,3)+rcs(x2,3)+x3+x4+x5,data=imputed1[[i]]) 
}


#plot the fitted splines
library(gridExtra) 
p.sp1<-ggplot(data.frame(x1=seq(-3,3, 0.1), 
              y=Predict(regression.splines[[i]], x1=seq(-3,3, 0.1), x2=0, x3=1, x4=1, x5=2)$yhat), 
              aes(x=x1, y=y)) + geom_smooth()
p.sp2<-ggplot(data.frame(x2=seq(-3,3, 0.1), 
              y=Predict(regression.splines[[i]], x1=1, x2=seq(-3,3, 0.1), x3=1, x4=1, x5=2)$yhat), 
              aes(x=x2, y=y)) + geom_smooth()
grid.arrange(p.sp1, p.sp2, ncol=2)

#predict for a new patient after a single fit or 
#averaging the predictions obtained from models developed in the multiply imputed datasets

prediction.ols <- function(new.patient, single.fit = NULL, multiple.fit = NULL){ 
  
  if(!is.null(multiple.fit)){
    mygrid <- expand.grid(k = 1:dim(new.patient)[1], i = 1:length(multiple.fit))  
    
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

prediction.ols(new.patient, multiple.fit = regression.splines)
complete.data <- data.cont[complete.cases(data.cont[,c(1:5,11)]),]
predicted.ols <- prediction.ols(complete.data, multiple.fit = regression.splines) # these are patients with fully observed data
ggplot(data.frame(observed=complete.data$y, predicted.ols=predicted.ols),
        aes(x=predicted.ols, y=observed)) +  
        geom_point(size=1) +
        geom_abline(intercept = 0, slope = 1, color="black",linetype="dashed", size=0.5) +
        xlim(-3,3) + ylim(-3,3)+geom_smooth(method=lm) 


###### Splines (alternative method) - 
# fit a model with splines using maximum likelihood
m2 <- fit.mult.impute(I(y)~ rcs(x1,3) + rcs(x2, 3) + x3 + x4 + x5, 
                      xtrans = a, data=data.cont, fitter=ols)
predict(m2, newdata=new.patient)
predicted.ols2 <- predict(m2, complete.data)


###### 2. penalized regression with smoothing splines ------ 
library(mgcv)
fit.gam=list()
for(i in 1:n.impute){
  fit.gam[[i]] <- gam(y ~ x3 + x4 + x5 + s(x1) + s(x2), data = imputed1[[i]])
}

#predict for a new patient after averaging the predictions 
#obtained from models developed in the multiply imputed datasets
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

prediction.gam(new.patient, multiple.fit = fit.gam)

new1 <- data.frame(x1=seq(-3,3, 0.2), x2=0, x3=1, x4=1, x5=2)
p.g1 <- ggplot(data.frame(x1=seq(-3,3, 0.2), 
              y=prediction.gam(new1, multiple.fit = fit.gam)), aes(x=x1, y=y)) + 
              geom_smooth()

new2 <- data.frame(x1=2, x2=seq(-3,3, 0.2), x3=1, x4=1, x5=2)
p.g2 <- ggplot(data.frame(x2=seq(-3,3, 0.2), 
              y=prediction.gam(new2, multiple.fit = fit.gam)), aes(x=x2, y=y)) + 
              geom_smooth()
grid.arrange(p.g1, p.g2, ncol=2)

predicted.gam <- prediction.gam(complete.data, multiple.fit = fit.gam)
ggplot(data.frame(observed=complete.data$y, 
      predicted.gam=predicted.gam),
      aes(x=predicted.gam, y=observed)) + geom_point(size=1) + 
      geom_abline(intercept = 0, slope = 1, color="black", linetype="dashed", size=0.5) +
      xlim(-3,3) + ylim(-3,3)+geom_smooth(method=lm) 



###### 3. Ridge regression and splines  ------ 
library(splines2)
library(glmnet)

lambdas <- 10^seq(2, -10, by = -0.3) 
#first find the position of the knots for the splines
bsMat.x1 <- bSpline(data.cont$x1[complete.cases(data.cont$x1)], 
                    knots = quantile(data.cont$x1, c(0.25, 0.5, 0.75), na.rm = TRUE))
bsMat.x2 <- bSpline(data.cont$x2[complete.cases(data.cont$x2)], 
                    knots = quantile(data.cont$x2, c(0.25, 0.5, 0.75), na.rm = TRUE))

fit.ridge <- list()
for( i in 1:n.impute){
  dfSplined.x1 <- as.data.frame(predict(bsMat.x1, imputed1[[i]]$x1))
  dfSplined.x2 <- as.data.frame(predict(bsMat.x2, imputed1[[i]]$x2))
  imp <- imputed1[[i]] 
  imp <- cbind(imp$y, imp$x3, imp$x4, imp$x5, dfSplined.x1, dfSplined.x2)
  colnames(imp) <- c("y", "x3", "x4", "x5", paste0("V", as.character(1:(length(colnames(imp))-4))))
  data_glmnet <- model.matrix(y ~., data = imp)
  X <- as.matrix(data_glmnet[,-1])
  colnames(X)[1:2] <- c("x3", "x4")
  Y <- imp$y  
  cvfit <- cv.glmnet(X,Y,family = "gaussian", alpha=0,
                     lambda = lambdas, nfolds=10)
  lambda.min <- cvfit$lambda.min
  fit.ridge[[i]] <- glmnet(X,Y,family = "gaussian", alpha=0, lambda = lambda.min)
}

# predict for new patients
prediction.ridge <- function(new.patient, single.fit = NULL, multiple.fit = NULL){
  
  if(!is.null(multiple.fit)){
    mygrid <- expand.grid(k = 1:dim(new.patient)[1],i = 1:length(multiple.fit))
    
    ff <- function(k,i){ 
      
      dfSplined.x1 <- as.data.frame(predict(bsMat.x1, new.patient$x1[k]))
      dfSplined.x2 <- as.data.frame(predict(bsMat.x2, new.patient$x2[k]))
      
      imp <- data.frame(x3=1*(new.patient$x3[k]==1),x4=1*(new.patient$x4[k]==1),x52=1*(new.patient$x5[k]==2),
                        x53=1*(new.patient$x5[k]==3), x54=1*(new.patient$x5[k]==4), dfSplined.x1, 
                        dfSplined.x2)
      colnames(imp)=c("x3","x4","x52","x53","x54", paste0("V", as.character(1:(length(colnames(imp))-5))))
      predict(multiple.fit[[i]], newx = as.matrix(imp))
    }
    prediction_matrix <- matrix(mapply(ff, mygrid$k, mygrid$i), 
                                nrow = dim(new.patient)[1], ncol = length(multiple.fit))
    prediction <- apply(prediction_matrix, 1, mean)
  } else if(!is.null(single.fit)){
    
    ff <- function(k){ 
      
      dfSplined.x1 <- as.data.frame(predict(bsMat.x1, new.patient$x1[k]))
      dfSplined.x2 <- as.data.frame(predict(bsMat.x2, new.patient$x2[k]))
      
      imp <- data.frame(x3=1*(new.patient$x3[k]==1),x4=1*(new.patient$x4[k]==1),x52=1*(new.patient$x5[k]==2),
                        x53=1*(new.patient$x5[k]==3), x54=1*(new.patient$x5[k]==4), dfSplined.x1, 
                        dfSplined.x2)
      colnames(imp)=c("x3","x4","x52","x53","x54", paste0("V", as.character(1:(length(colnames(imp))-5))))
      
      predict(single.fit, newx = as.matrix(imp))
    }
    prediction <- sapply(1:dim(new.patient)[1], ff)
  }
  return(prediction)    
}
  
prediction.ridge(new.patient, multiple.fit = fit.ridge)
predicted.ridge <- prediction.ridge(complete.data, multiple.fit = fit.ridge)

ggplot(data.frame(observed=complete.data$y, 
                  predicted.ridge=predicted.ridge),
       aes(x=predicted.ridge, y=observed))+ geom_point(size=1)+ 
       geom_abline(intercept = 0, slope = 1, color="black",
       linetype="dashed", size=0.5) +
       xlim(-3,3) + ylim(-3,3)+geom_smooth(method=lm) 

#compare predictions
p1<-ggplot(data.frame(predicted.gam=predicted.gam, predicted.ols=predicted.ols),
        aes(x=predicted.ols, y=predicted.gam))+  geom_point(size=1)+ 
        geom_abline(intercept = 0, slope = 1, color="black",linetype="dashed", size=0.5)+
        xlim(-4.,4.)+ylim(-4.,4.)

p2<-ggplot(data.frame(predicted.ols=predicted.ols, predicted.ridge=predicted.ridge),
        aes(x=predicted.ridge, y=predicted.ols))+  geom_point(size=1)+ 
        geom_abline(intercept = 0, slope = 1, color="black",linetype="dashed", size=0.5)+
        xlim(-4.,4.)+ylim(-4.,4.)

p3<-ggplot(data.frame(predicted.gam=predicted.gam, predicted.ridge=predicted.ridge),
        aes(x=predicted.ridge, y=predicted.gam))+  geom_point(size=1)+ 
        geom_abline(intercept = 0, slope = 1, color="black",linetype="dashed", size=0.5)+
        xlim(-4.,4.)+ylim(-4.,4.)

grid.arrange(p1, p2, p3, ncol=3)

# Assess apparent performance (no optimism correction) -------------------
calculate_performance <- function(observed, predicted){
  
  MAE <- mean(abs(observed - predicted))
  MSE <- mean((observed - predicted)^2)
  R2 <- summary(lm(observed~predicted))$r.squared
  
  vec <- c(MAE, MSE, R2)
  names(vec) <- c("MAE", "MSE", "R2")
  return(vec)
}

apparent.ols <- calculate_performance(complete.data$y, predicted.ols)
apparent.gam <- calculate_performance(complete.data$y, predicted.gam) 
apparent.ridge <- calculate_performance(complete.data$y, predicted.ridge)
round(rbind(apparent.ols, apparent.gam, apparent.ridge), digits=2)

# Internal CV via bootstrapping -------------------
# bootstrap in each multiply imputed dataset
N.bootstrap <- 5
optimism.ols.eachbootstrap <- optimism.gam.eachbootstrap <- optimism.ridge.eachbootstrap <- matrix(NA, N.bootstrap, 3)
optimism.ols <- optimism.gam <- optimism.ridge <- matrix(NA, n.impute, 3)

for (i in 1:n.impute){
  for(j in 1:N.bootstrap){
    
    boot.sample <- sample(length(imputed1[[i]]$y),replace = T) 
    
    # create bootstrap sample
    imp.boot <- lapply(imputed1[[i]], function(x){x[boot.sample]}) 
    
    regression.splines.boot <- ols(y~rcs(x1,3)+rcs(x2,3)+x3+x4+x5,data= imp.boot) 
    
    fit.gam.boot <- gam(y ~ x3 + x4 + x5 + s(x1) + s(x2), data = imp.boot)
    
    dfSplined.x1 <- as.data.frame(predict(bsMat.x1, imp.boot$x1))
    dfSplined.x2 <- as.data.frame(predict(bsMat.x2, imp.boot$x2))
    imp <- cbind(imp.boot$y, imp.boot$x3, imp.boot$x4, imp.boot$x5, dfSplined.x1, dfSplined.x2)
    colnames(imp)=c("y", "x3", "x4", "x5", paste0("V", as.character(1:(length(colnames(imp))-4))))
    data_glmnet <- model.matrix(y ~ ., data = imp)
    X <- as.matrix(data_glmnet[,-1])
    colnames(X)[1:2] <- c("x3", "x4")
    Y <- imp$y  
    cvfit <- cv.glmnet(X,Y,family = "gaussian", alpha=0,
                       lambda = lambdas, nfolds=10)
    lambda.min <- cvfit$lambda.min
    fit.ridge.boot <- glmnet(X,Y,family = "gaussian", alpha=0, lambda = lambda.min)
    
    # predict in bootstrap
    f1 <- as.data.frame(do.call(cbind, lapply(imp.boot, function(x) {as.numeric(as.character(x))})))
    
    boot.prediction.ols <- prediction.ols(f1, single.fit = regression.splines.boot)
    boot.prediction.gam <- prediction.gam(f1, single.fit = fit.gam.boot)
    boot.prediction.ridge <- prediction.ridge(f1, single.fit = fit.ridge.boot)
    
    ols.boot <- calculate_performance(f1$y, boot.prediction.ols)
    gam.boot <- calculate_performance(f1$y, boot.prediction.gam)
    ridge.boot <- calculate_performance(f1$y, boot.prediction.ridge)
    
    # predict in test data
    f2 <- as.data.frame(do.call(cbind, lapply(imputed1[[i]], function(x) {as.numeric(as.character(x))})))
    test.prediction.ols <- prediction.ols(f2, single.fit = regression.splines.boot)
    test.prediction.gam <- prediction.gam(f2, single.fit = fit.gam.boot)
    test.prediction.ridge <- prediction.ridge(f2, single.fit = fit.ridge.boot)
    
    ols.test <- calculate_performance(f2$y, test.prediction.ols)
    gam.test <- calculate_performance(f2$y, test.prediction.gam)
    ridge.test <- calculate_performance(f2$y, test.prediction.ridge)
    
    optimism.ols.eachbootstrap[j,] <- ols.boot - ols.test
    optimism.gam.eachbootstrap[j,] <- gam.boot - gam.test
    optimism.ridge.eachbootstrap[j,] <- ridge.boot - ridge.test
  }
  optimism.ols[i,] <- apply(optimism.ols.eachbootstrap, 2, mean)
  optimism.gam[i,] <- apply(optimism.gam.eachbootstrap, 2, mean)
  optimism.ridge[i,] <- apply(optimism.ridge.eachbootstrap, 2, mean)
  
  print(paste0("imputation done: ", i))
}

mean.optimism.ols <- apply(optimism.ols, 2, mean)
mean.optimism.gam <- apply(optimism.gam, 2, mean)
mean.optimism.ridge <- apply(optimism.ridge, 2, mean)

optimism.corrected.ols <- apparent.ols - mean.optimism.ols
optimism.corrected.gam <- apparent.gam - mean.optimism.gam
optimism.corrected.ridge <- apparent.ridge - mean.optimism.ridge

round(rbind(mean.optimism.ols, mean.optimism.gam, mean.optimism.ridge), digits=2)
round(rbind(optimism.corrected.ols, optimism.corrected.gam, optimism.corrected.ridge), digits=2)

# Internal-external CV -------------------
clusters <- unique(data.cont$clust)
N.clust <- length(clusters) # 5 clusters in this example
data.in <- data.leftout <- list()

#create the datasets
for(i in 1:N.clust){
  data.in[[i]]<-data.cont[data.cont$clust!=clusters[i],]
  data.leftout[[i]]<-data.cont[data.cont$clust==clusters[i],]
  complete.index <- complete.cases(data.leftout[[i]][,c(paste0("x", 1:5), "y")])
  data.leftout[[i]] <- data.leftout[[i]][complete.index,]
}

#impute the data and fit the model
n.impute <- 5

imputed <- regression.splines.CV <- fit.gam.CV <- fit.ridge.CV <- list()
leftout.prediction.ols <- leftout.prediction.gam <- leftout.prediction.ridge <- list()
leftout.performance.ols <- leftout.performance.gam <- leftout.performance.ridge <- list()

for (i in 1:N.clust){
  a <- aregImpute(data=data.in[[i]], I(y)~x1+x2+I(x3)+I(x4)+I(x5), n.impute=n.impute, nk=3, match='closest')

  for (j in 1:n.impute){
    imputed[[j]] <- impute.transcan(a, imputation=j, data=data.in[[i]], list.out=TRUE,
                                  pr=FALSE, check=FALSE)
    regression.splines.CV[[j]]<- ols(y~rcs(x1,3)+rcs(x2,3)+x3+x4+x5,data=imputed[[j]]) 
    fit.gam.CV[[j]] <- gam(y ~ x3+x4+x5+s(x1)+s(x2), data = imputed[[j]])

    dfSplined.x1 <- as.data.frame(predict(bsMat.x1, imputed[[j]]$x1))
    dfSplined.x2 <- as.data.frame(predict(bsMat.x2, imputed[[j]]$x2))
    imp <- cbind(imputed[[j]]$y, imputed[[j]]$x3, imputed[[j]]$x4, imputed[[j]]$x5, dfSplined.x1, dfSplined.x2)
    colnames(imp) <- c("y", "x3", "x4", "x5", paste0("V", as.character(1:(length(colnames(imp))-4))))
    data_glmnet <- model.matrix(y ~ ., data = imp)
    X <- as.matrix(data_glmnet[,-1])
    colnames(X)[1:2] <- c("x3", "x4")
    Y <- imp$y  
    cvfit <- cv.glmnet(X,Y,family = "gaussian", alpha=0,
                       lambda = lambdas, nfolds=10)
    lambda.min <- cvfit$lambda.min
    fit.ridge.CV[[j]] <- glmnet(X,Y,family = "gaussian", alpha=0, lambda = lambda.min)
  }
  
  leftout.prediction.ols[[i]] <- prediction.ols(data.leftout[[i]], multiple.fit = regression.splines.CV)
  leftout.prediction.gam[[i]] <- prediction.gam(data.leftout[[i]], multiple.fit = fit.gam.CV)
  leftout.prediction.ridge[[i]] <- prediction.ridge(data.leftout[[i]], multiple.fit = fit.ridge.CV)

  leftout.performance.ols[[i]] <- calculate_performance(data.leftout[[i]]$y, leftout.prediction.ols[[i]])
  leftout.performance.gam[[i]] <- calculate_performance(data.leftout[[i]]$y, leftout.prediction.gam[[i]])
  leftout.performance.ridge[[i]] <- calculate_performance(data.leftout[[i]]$y, leftout.prediction.ridge[[i]])
}

# performance per cluster
leftout.performance.ols
leftout.performance.gam
leftout.performance.ridge

per.cluster.ols=matrix(leftout.performance.ols[[1]], nrow=1)
for(i in 2: N.clust){
  per.cluster.ols=rbind(per.cluster.ols, matrix(leftout.performance.ols[[i]], nrow=1))
  }
per.cluster.ols=data.frame(per.cluster.ols)
colnames(per.cluster.ols)=c("MAE", "MSE", "R2")
rownames(per.cluster.ols)=paste("cluster",clusters,":")
round(per.cluster.ols, digits=2)

# performance aggregating all clusters
leftout.prediction.ols <- do.call(c, leftout.prediction.ols)
leftout.prediction.gam <- do.call(c, leftout.prediction.gam)
leftout.prediction.ridge <- do.call(c, leftout.prediction.ridge)
IECV.observed <- do.call(rbind, data.leftout)$y
IECV.cluster <- do.call(rbind, data.leftout)$clust

IECV.ols <- calculate_performance(IECV.observed, leftout.prediction.ols)
IECV.gam <- calculate_performance(IECV.observed, leftout.prediction.gam)
IECV.ridge <- calculate_performance(IECV.observed, leftout.prediction.ridge)
round(rbind(IECV.ols, IECV.gam, IECV.ridge), digits=2)

IECV <- data.frame(IECV.observed = IECV.observed, cluster = IECV.cluster, 
                   leftout.prediction.ols = leftout.prediction.ols, leftout.prediction.gam = leftout.prediction.gam, 
                   leftout.prediction.ridge = leftout.prediction.ridge)

p4<-ggplot(IECV,
       aes(x=leftout.prediction.ols, y=IECV.observed, color=cluster))+  geom_point(size=2)+ 
      geom_abline(intercept = 0, slope = 1, color="black",linetype="dashed", size=0.5)+
      xlim(-4,4)+ylim(-4,4) + theme(legend.position = "none")+
  geom_smooth(method=lm, se=F) +ylab("observed")

p5<-ggplot(IECV,
       aes(x=leftout.prediction.gam, y=IECV.observed, color=cluster))+  geom_point(size=2)+ 
  geom_abline(intercept = 0, slope = 1, color="black",linetype="dashed", size=0.5)+
  xlim(-4,4)+ylim(-4,4) + theme(legend.position = "none")+
  geom_smooth(method=lm, se=F) +ylab("observed")

p6<-ggplot(IECV,
       aes(x=leftout.prediction.ridge, y=IECV.observed, color=cluster))+  geom_point(size=2)+ 
        geom_abline(intercept = 0, slope = 1, color="black",linetype="dashed", size=0.5)+
  xlim(-4,4)+ylim(-4,4)+geom_smooth(method=lm, se=F) +ylab("observed")

library(ggpubr)
ggarrange(p4, p5, p6, ncol=3, common.legend = TRUE, legend="bottom")

setwd("G:/My Drive/PROJECT/practical guide to prediction models/R code/continuous outcome")
#save.image(file="cont_results")
load(file="cont_results")
