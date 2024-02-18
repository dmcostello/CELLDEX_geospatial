# Make training and test sets by 80-20 split 
library(caret)
set.seed(3)
# For the most rigorous test, we want to split 20% of the partners off for testing
Cdat_val <- cbind(Cdat,partnerid=C_stb$partnerid)
length(unique(Cdat_val$partnerid)) *0.2 #26 partners will be in the test dataset

test.partners <- sample(unique(Cdat_val$partnerid),26,replace = F)

#Split the data set 80-20 into train and test sets
train <- Cdat_val[!(Cdat_val$partnerid %in% test.partners),]
test  <- Cdat_val[Cdat_val$partnerid %in% test.partners,]

# Run boosted regression tree model with Gaussian error (e.g. linear regression)
Cgbm<- gbm(logk~., 
           data=train[,-103], 
           distribution="gaussian",
           n.trees=20000,
           shrinkage=0.001,
           interaction.depth=5,
           cv.folds=20)

# Check performance
(best.iter <- gbm.perf(Cgbm,method="cv"))

# Plot variable influence based on the estimated best number of trees
sum<-summary(Cgbm,n.trees=best.iter,method=permutation.test.gbm) 

head(sum,20)

# Calculate pseudo-R2
print(1-sum((test$logk - predict(Cgbm, newdata=test[,-103], n.trees=best.iter))^2)/
        sum((test$logk - mean(test$logk))^2))

# Calculate RMSE
library(Metrics)
print(rmse(test$logk, predict(Cgbm, newdata=test[,-103], n.trees=best.iter)))
print(rmse(train$logk, predict(Cgbm, newdata=train[,-103], n.trees=best.iter)))

# Leave-one-out cross validation
partners <- unique(Cdat_val$partnerid)

cv.mse <- rep(0,length(partners))
for(i in 1:length(partners)){
  subdata <- Cdat_val[!(Cdat_val$partnerid %in% partners[i]),]
  outdata <- Cdat_val[(Cdat_val$partnerid %in% partners[i]),]

  loo_mod<- gbm(logk~., 
             data=subdata[,-103], 
             distribution="gaussian",
             n.trees=20000,
             shrinkage=0.001,
             interaction.depth=5,
             cv.folds=20)
  
  loo_it <- gbm.perf(loo_mod,method="cv")
  loo.rmse[i] <- rmse(outdata$logk,predict(loo_mod, newdata=outdata[,-103], n.trees=loo_it))
  print(loo.rmse[i])
  print(paste("step =", i))
}

write.csv(cbind(loo.rmse,partners),"~/Desktop/loo.rmse.csv")
