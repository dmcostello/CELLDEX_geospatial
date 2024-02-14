# Make training and test sets by 80-20 split 
library(caret)
set.seed(3)
# Ideally, we want to split by locations, not actual cotton strips
# So we may want to wait to delete partnerid until running this split step
trainIndex <- createDataPartition(Cdat$partnerid, p = .8, list = FALSE, times = 1)

train <- Cdat[ trainIndex,]
test  <- Cdat[-trainIndex,]

# Run boosted regression tree model with Gaussian error (e.g. linear regression)
Cgbm<- gbm(logk~., 
           data=train, 
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
print(1-sum((test$logk - predict(Cgbm, newdata=test, n.trees=best.iter))^2)/
        sum((test$logk - mean(test$logk))^2))

# Calculate RMSE
library(Metrics)
print(rmse(test$logk, predict(Cgbm, newdata=test, n.trees=best.iter)))