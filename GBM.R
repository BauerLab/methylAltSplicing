
celltype.l <- c('IMR90', 'Gm12878',  'H1hesc')   # the celltype of the data (IMR90, Gm12878, H1hesc)

cellnames <- c('IMR-90', 'GM12878', 'H1-hESC')

#################################
## Load and install libraries  ##
#################################

# for the plots 
library(gbm)    # GBM package
library(data.table)
library(prodlim)      # for row matching
library(ggplot2) 
library(gplots)       # for heatmap.2
library(Rmisc)
library(ROCR)

#################################
## Function for the Training   ##
#################################

downSample = function(data.df){
  rows_inc = which(data.df[,1] == 1)
  rows_exc = which(data.df[,1] == 0)
  set.seed(123)
  downsample_inc = sample(rows_inc, length(rows_exc))
  downsampled_data.df = data.df[c(rows_exc,downsample_inc),]
  return(downsampled_data.df)
}

shuffle <- function(set){
  set.seed(123)
  newset <- set[sample(1:nrow(set), nrow(set)),]
  return(newset)
}

#################################################
# Split data in Train, Validation and Test set  #
#################################################

# print(paste('[NOTE] number of samples in the whole feature dataset =', nrow(data.df)))

# this function is an addition to the splitData function below
createDatasets <- function(ase, data.df, pTrain, pValid, setSeed, setBagging) {
  # take group of exclusion or inclusions
  rows.l <- which(data.df[,1] == ase)
  
  # length of group 
  nrows <- length(rows.l)
  
  # take samples for training 
  if ( setSeed ) {
    set.seed(123)
    rows.train.l <- sample(rows.l, round(pTrain*nrows), replace=setBagging)
  } else {
    rows.train.l <- sample(rows.l, round(pTrain*nrows), replace=setBagging)
  }
  
  # remove rows from whole set 
  red1.l <- rows.l[-which(rows.l %in% rows.train.l)]
  
  # look if validation set is needed 
  red2.l <- red1.l
  rows.valid.l <- as.numeric(c())
  if ( pValid != 0.0 ) {
    
    # take samples for validation set 
    if ( setSeed ) {
      set.seed(123)
      rows.valid.l <- sample(red1.l, round(pValid*nrows))
    } else {
      rows.valid.l <- sample(red1.l, round(pValid*nrows))
    }
    
    # remove rows from list 
    red2.l <- red1.l[-which(red1.l %in% rows.valid.l)]
  }
  
  # the rest is the test set 
  rows.test.l <- red2.l
  
  print(paste0('[NOTE] number of samples for train = ', length(rows.train.l)))
  print(paste0('[NOTE] number of samples for validation = ', length(rows.valid.l)))
  print(paste0('[NOTE] number of samples for test = ', length(rows.test.l)))
  
  # take sum for control 
  sumRows <- length(unique(rows.train.l)) + length(rows.valid.l) + length(rows.test.l)
  
  # test if the sum of the datasets equals the whole dataset 
  if ( sumRows != nrows ) {
    print(paste('[ERROR] something wrong with number of data sets',sumRows,'!=',nrows))
  }
  
  # test if the individual rows are unique for each data set 
  if ( length(intersect(intersect(rows.train.l,rows.valid.l),rows.test.l)) != 0 ) {
    print(paste('[ERROR] something wrong with the individual lines for the datasets'))
  }
  
  return(list(rows.train.l, rows.valid.l, rows.test.l))
}

# this function split the whole population in train,validation and test dataset
# for inclusion as well as exclusion
splitData <- function(data.df, pTrain, pValid, loadTest, setSeed, setBagging) {
  
  # split dataset (train, validation, test)
  print('[NOTE] split for exclusions')
  rows.excl.l <- createDatasets(0, data.df, pTrain, pValid, setSeed, setBagging)
  print('[NOTE] split for inclusions')
  rows.incl.l <- createDatasets(1, data.df, pTrain, pValid, setSeed, setBagging)
  
  #data.df[,1] = otherfeatures.df[,12]
  
  # collect training data
  features.train.df <- data.df[rows.excl.l[[1]],]
  features.train.df <- rbind(features.train.df, data.df[rows.incl.l[[1]],])
  
  # collect validation data
  features.validation.df <- c()
  if ( pValid != .0 ) {
    features.validation.df <- data.df[rows.excl.l[[2]],]
    features.validation.df <- rbind(features.validation.df, data.df[rows.incl.l[[2]],])
  }
  
  # collect test data
  features.test.df <- data.df[rows.excl.l[[3]],]
  features.test.df <- rbind(features.test.df, data.df[rows.incl.l[[3]],])
  
  # rownames changing if you sample with replacement
  if ( setBagging == TRUE ){
    rownames_training = gsub( '\\..*', '', rownames(features.train.df) )
  } else {
    rownames_training = rownames(features.train.df)
  }
  
  # get labels
  lables.train <- data.df[rownames_training, 1]
  lables.test <- data.df[rownames(features.test.df), 1]
  labels.valid <- data.df[rownames(features.validation.df), 1]
  
  train.df <- as.data.frame(features.train.df)
  test.df <- as.data.frame(features.test.df)
  valid.df <- as.data.frame(features.validation.df)
  
  y <- "target"
  x <- setdiff(names(train.df), y)
  
  train.df[,y] <- as.character(lables.train)
  test.df[,y] <- as.character(lables.test)
  valid.df[,y] <- as.character(labels.valid)
  
  return(list(train.df, valid.df, test.df, y, x,
              rownames(features.train.df), rownames(features.validation.df), rownames(features.test.df)))
}

# function used for k-fold crossvalidation
kfold.cv <- function(k, cv.df) {
  rownames.l <- rownames(cv.df)
  
  # get groups of exclusions and inclusions
  exc.l <- which(cv.df[,1] == '0')
  inc.l <- which(cv.df[,1] == '1')
  
  # get the length of each k-validation chunk
  lengthchunk <- round(length(exc.l)*1/k)
  
  chunks <- list()
  for ( i in 1:(k-1) ){
    # sample for exlsuions and inclusions
    rows1 <- sample( exc.l, lengthchunk ) 
    rows2 <- sample( inc.l, lengthchunk )
    
    # combine 
    all_rows <- c(rows1,rows2)
    
    # get chunk 
    chunks[[i]] <- rownames.l[all_rows]
    
    # remove rows 
    rownames.l <- rownames.l[-all_rows]
    
    # define new list 
    exc.l <- which(grepl('ex', rownames.l))
    inc.l <- which(grepl('in', rownames.l))
  }
  
  # last chunk 
  chunks[[k]] <- rownames.l
  
  # check if everthing add up
  if ( length(unique(unlist(chunks))) != nrow(cv.df) ){
    print('[NOTE] error something went wrong in the cross validation splitting')
  }
  return(chunks)
}

# function to capture for each epoch the parameter list and the 
# training and validation error 
mx.callback.params <- function(period, logger = NULL) {
  function (iteration, nbatch, env, verbose = TRUE) 
  {
    if (nbatch%%period == 0 && !is.null(env$metric)) {
      result <- env$metric$get(env$train.metric)
      if (nbatch != 0 & verbose) 
        cat(paste0("Batch [", nbatch, "] Train-", result$name, 
                   "=", result$value, "\n"))
      if (!is.null(logger)) {
        if (class(logger) != "mx.metric.logger") {
          stop("Invalid mx.metric.logger.")
        }
        logger$train <- c(logger$train, result$value)
        if (!is.null(env$eval.metric)) {
          result <- env$metric$get(env$eval.metric)
          if (nbatch != 0 & verbose) 
            cat(paste0("Batch [", nbatch, "] Validation-", 
                       result$name, "=", result$value, "\n"))
          logger$eval <- c(logger$eval, result$value)
        }
      }
    }
    if ( iteration %% period == 0 ){
      logList <<- c(logList, list(env$model$arg.params))
    }
    return(TRUE)
  }
}

model_ANN <- function(params, train.m, valid.m, rounds, batch, h){
  mx.set.seed(123)
  err <- mx.metric.logger$new()

  dropout <- 0.0

  data <- mx.symbol.Variable('data')
  fc1 <- mx.symbol.FullyConnected(data, name='fc1', num_hidden=h, dropout = dropout)
  act1 <- mx.symbol.Activation(fc1, name='rel1', act_type='relu')
  fc2 <- mx.symbol.FullyConnected(act1, name='fc2', num_hidden=2, dropout = dropout)
  softmax = mx.symbol.SoftmaxOutput(fc2, name='sm')
  model <- NULL

  lern_rate <- 0.01
  reg <- 0.001
  momentum <- 0.9

  if( length(params) == 0){
    if ( is.null(valid.m) ) {
      model <- mx.model.FeedForward.create(softmax,X=train.m[,-1], y=train.m[,1],
                                           array.batch.size=batch,
                                           learning.rate=lern_rate, wd = reg, momentum = momentum, num.round=rounds,
                                           initializer = mx.init.normal(1), optimizer = "sgd", eval.metric=mx.metric.accuracy,
                                           epoch.end.callback = mx.callback.params(1, err),
                                           ctx = mx.cpu(), verbose = TRUE)
    } else {
      model <- mx.model.FeedForward.create(softmax,X=train.m[,-1], y=train.m[,1],
                                           eval.data=list(data=valid.m[,-1], label=valid.m[,1]),array.batch.size=batch,
                                           learning.rate=lern_rate, wd = reg, momentum = momentum, num.round=rounds,
                                           initializer = mx.init.normal(1), optimizer = "sgd", eval.metric=mx.metric.accuracy,
                                           epoch.end.callback = mx.callback.params(1, err),
                                           ctx = mx.cpu(), verbose = TRUE)
    }
  } else {
    print('[NOTE] Use Defined Weight and Bias List')
    if ( is.null(valid.m) ) {
      model <- mx.model.FeedForward.create(softmax,X=train.m[,-1], y=train.m[,1],
                                           array.batch.size=batch,
                                           learning.rate=lern_rate, wd = reg, momentum = momentum, num.round=rounds,
                                           initializer = mx.init.normal(1), optimizer = "sgd", eval.metric=mx.metric.accuracy,
                                           epoch.end.callback = mx.callback.params(1, err),
                                           ctx = mx.cpu(), verbose = TRUE, arg.params = params)
    } else {
      model <- mx.model.FeedForward.create(softmax, X=train.m[,-1], y=train.m[,1],
                                           eval.data=list(data=valid.m[,-1], label=valid.m[,1]),array.batch.size=batch,
                                           learning.rate=lern_rate, wd = reg, momentum = momentum, num.round=rounds,
                                           initializer = mx.init.normal(1), optimizer = "sgd", eval.metric=mx.metric.accuracy,
                                           epoch.end.callback = mx.callback.params(1, err),
                                           ctx = mx.cpu(), verbose = TRUE, arg.params = params)
    }
  }

  if ( is.null(valid.m) == FALSE ){
    print(which.max(err$eval))
    print(length(logList))
    bestparams <- logList[[which.max(err$eval)]]
    model$arg.params <- bestparams
  }

  return(list(Model = model, Error = err))
}

# collect for all three cell types
c_all.data.df <- list()
all.data.df <- list()
all.positMatrix <- list()
correctSamples.l <- list()

# parameter set of the model
ann_parameters <- list()

# list which keep parametrs and errors of the model
listmodel <- NULL

binfeatures <- function(features, pre){
  listmeans <- matrix(nrow = nrow(features), ncol = 50)
  col <- 1
  for ( i in 1:50 ){
    buff <- features[,c(col:(col+9))]
    listmeans[,i] <- apply(buff,1, mean)
    col <- col + 10
  }
  colnames(listmeans) <- paste0(pre, c(1:50))
  return(listmeans)
}

# models list
models.l <- list()

# accuracy list
acc.l <- list()

################
## Load data  ##
################

for ( k in 1:3 ) {
  data.df <- as.data.frame(fread(paste0('data/network_data_',celltype.l[k],'.tsv'), sep =  "\t", header = TRUE))
  
  # first col is the rownames of data.df
  rownames(data.df) <- data.df[,1]
  data.df <- data.df[,-1]
  
  # get start and end of exon
  positMatrix <- matrix(nrow=nrow(data.df), ncol = 4)
  positMatrix[,1] <- data.df[,(ncol(data.df)-3)]
  positMatrix[,2] <- data.df[,(ncol(data.df)-2)]
  positMatrix[,3] <- data.df[,(ncol(data.df)-1)]
  positMatrix[,4] <- data.df[,ncol(data.df)]
  data.df <- data.df[,-c((ncol(data.df)-3):ncol(data.df))]
  colnames(positMatrix) <- c('Chr', 'Start', 'End', 'Strand')
  rownames(positMatrix) <- rownames(data.df)
  
  buff <- data.df
  data.df <- sapply(data.df, as.numeric)
  rownames(data.df) <- rownames(buff)
  
  features <- data.df[,paste0('up',c(1:500))]
  binned_Up <- binfeatures(features, 'p')
  features <- data.df[,paste0('do',c(1:500))]
  binned_Do <- binfeatures(features, 'd')
  
  pickedcols <- c('target', paste0('ex',c(2:49))))
  picked_data.df <- data.df[,pickedcols]
  colnames(picked_data.df)
  
  picked_data.df <- cbind(picked_data.df, binned_Up[,-50])
  picked_data.df <- cbind(picked_data.df, binned_Do[,-1])
  
  all.data.df[[k]] <- picked_data.df
  
  pickedcols <- c('target', paste0('up',c(1:490)), paste0('do',c(11:500)),
                            paste0('ex',c(2:49)))
  c_all.data.df[[k]] <- data.df[,pickedcols]
  print(colnames(c_all.data.df[[k]]))
  
  all.positMatrix[[k]] <- positMatrix
  
}

rm(data.df)
rm(buff)
rm(binned_Up)

###################
## Prepare data  ##
###################

# picking out samples which changes ASE
changingASE.l <- NULL

c1.df = all.positMatrix[[1]]
c2.df = all.positMatrix[[2]]
c3.df = all.positMatrix[[3]]

find_changedASE <- function(set1.df, set2.df) {
  row.matches12 <- row.match(as.data.frame(set1.df), as.data.frame(set2.df))
  d1 <- set1.df[which(!is.na(row.matches12)),]
  d2 <- set2.df[row.matches12[which(!is.na(row.matches12))],]
  s1 <- d1[grep('ex', rownames(d1)),]
  s2 <- d2[grep('ex', rownames(d2)),]
  row.matches12 <- row.match(as.data.frame(s1), as.data.frame(s2))
  e1 <- s1[which(is.na(row.matches12)),]
  e2 <- d2[grep('in', rownames(d2)),]
  row.matches12 <- row.match(as.data.frame(e1), as.data.frame(e2))
  if( length(which(is.na(row.matches12))) != 0 ){
    print('[ERROR] something went wrong with finding changed ASE')
  }
  e2 <- e2[row.matches12,]
  return( list(Set1 = rownames(e1),  Set2 = rownames(e2)) )
}

find_changedASE_v2 <- function(set1.df, set2.df) {
  row.matches12 <- row.match(as.data.frame(set1.df), as.data.frame(set2.df))
  d1 <- set1.df[which(!is.na(row.matches12)),]
  d2 <- set2.df[row.matches12[which(!is.na(row.matches12))],]
  
  # check first for exclusions
  s1 <- d1[grep('ex', rownames(d1)),]
  s2 <- d2[grep('ex', rownames(d2)),]
  row.matches12 <- row.match(as.data.frame(s1), as.data.frame(s2))
  e1 <- s1[which(is.na(row.matches12)),]
  e2 <- d2[grep('in', rownames(d2)),]
  row.matches12 <- row.match(as.data.frame(e1), as.data.frame(e2))
  if( length(which(is.na(row.matches12))) != 0 ){
    print('[ERROR] something went wrong with finding changed ASE exclusions')
  }
  e2 <- e2[row.matches12,]
  
  Set1 = rownames(e1)
  Set2 = rownames(e2)
  
  # check now for inclusions
  s1 <- d1[grep('in', rownames(d1)),]
  s2 <- d2[grep('in', rownames(d2)),]
  row.matches12 <- row.match(as.data.frame(s1), as.data.frame(s2))
  e1 <- s1[which(is.na(row.matches12)),]
  e2 <- d2[grep('ex', rownames(d2)),]
  row.matches12 <- row.match(as.data.frame(e1), as.data.frame(e2))
  if( length(which(is.na(row.matches12))) != 0 ){
    print('[ERROR] something went wrong with finding changed ASE inclusions')
  }
  e2 <- e2[row.matches12,]
  
  Set3 = rownames(e1)
  Set4 = rownames(e2)
  
  return( list(Set1 = Set1,  Set2 = Set2, Set3 = Set3, Set4 = Set4) )
}

set1v2 <- find_changedASE_v2(c1.df, c2.df) 
set1v3 <- find_changedASE_v2(c1.df, c3.df) 
set2v3 <- find_changedASE_v2(c2.df, c3.df) 

h1.df <- all.data.df[[1]]
h2.df <- all.data.df[[2]]
h3.df <- all.data.df[[3]]

changingASE.l[[1]] <- rbind(h1.df[set1v2$Set1,], h1.df[set1v3$Set1,])
changingASE.l[[2]] <- rbind(h2.df[set1v2$Set2,])
changingASE.l[[3]] <- rbind(h3.df[set1v3$Set2,])

for ( k in 1:3 ){
  set <- changingASE.l[[k]]
  changingASE.l[[k]] <- set[which(duplicated(rownames(set)) == FALSE),] 
}

h.df <- list()

h.df[[1]] <- h1.df[which(!rownames(h1.df) %in% c(set1v2$Set1,set1v3$Set1)),]
h.df[[2]] <- h2.df[which(!rownames(h2.df) %in% c(set1v2$Set2)),]
h.df[[3]] <- h3.df[which(!rownames(h3.df) %in% c(set1v3$Set2)),]

rm(h1.df)
rm(h2.df)
rm(h3.df)

all.train.l <- list()
all.pretrain.m <- list()
incorrect.l <- list()
auc.l <- numeric(3)
perf.l <- list()
logList <- NULL

plusChangedASE <- TRUE
params <-  c(1508, 0.001616662, 19, 36)

acc.l <- numeric(3)
for ( k in 1:3 ){

  ##################
  ## Train Model  ##
  ##################
  
  datasets.l <- splitData(downSample(h.df[[k]]), 0.7, 0.0, TRUE, TRUE, FALSE)
  
  train.m <- datasets.l[[1]]
  holdout.test.m <- datasets.l[[3]]
  
  train.m <- shuffle(train.m)
  holdout.test.m <- shuffle(holdout.test.m)
  
  if ( plusChangedASE ) {
    holdout.test.m <- rbind(holdout.test.m, changingASE.l[[k]])
  }

  print('Train Model')
  model <- gbm(target ~ ., distribution = "bernoulli",
               data = train.m, n.trees = params[1], shrinkage = params[2], 
               interaction.depth = params[3], n.minobsinnode = params[4])
  
  gbm.prediction <- predict(model, holdout.test.m[,-1], n.trees = params[1], type="response")
  print('FINISH')
  
  # round to get the classes
  finalprediction <- gbm.prediction
  finalprediction[which(gbm.prediction > 0.5)] <- 1
  finalprediction[which(gbm.prediction <= 0.5)] <- 0
  
  confusionMatrix <- table(holdout.test.m[,1], finalprediction)
  
  print(confusionMatrix)
  
  positiveRates <- numeric(ncol(confusionMatrix))
  for( i in 1:ncol(confusionMatrix)){
    row <- as.numeric(colnames(confusionMatrix))[i]+1
    positiveRates[i] <- confusionMatrix[row,i]/sum(confusionMatrix[row,])
  }
  
  print(positiveRates)
  
  if ( ncol(confusionMatrix) == nrow(confusionMatrix)  ) {
    acc.l[k] <- ((confusionMatrix[1,1] + confusionMatrix[2,2]) / 
            ( sum(confusionMatrix[1,]) + sum(confusionMatrix[2,]) ))
  } else {
    print('[NOTE] Confusion matrix misses some labels')
  }

  # save correct classified labels
  correctSamples.l[[k]] <- rownames(holdout.test.m)[which(holdout.test.m[,1] == finalprediction)]
  incorrect.l[[k]] <- rownames(holdout.test.m)[which(holdout.test.m[,1] != finalprediction)]

  roc.pred <- prediction(gbm.prediction, as.numeric(holdout.test.m[,1]), label.ordering = NULL)
  perf.l[[k]] <- performance(roc.pred, 'tpr', 'fpr')
  auc.l[k] <- performance(roc.pred, 'auc')@y.values[[1]]
  print(auc.l)
}

print(acc.l)
print(auc.l)

for ( k in 1:3 ) {
  perf <- perf.l[[k]]
  
  postscript(file = paste0('plots/GBM_performance_',celltype.l[k],'.eps'), 
      width=800, height=800, family="Times")
  par(cex = 2.0)
  plot(perf, xlab = 'False Positive Rate', ylab = 'True Positive Rate')
  acc.l <- unlist(acc.l)
  title(main = paste(cellnames[k],'\n Accuracy =',round(acc.l[k], digits = 3)
                     ,'AUC =', round(auc.l[k], digits = 3) ) )
  abline(a=0,b=1)
  dev.off()
}

#############################################################
## Get Samples Interesection of Correct Classified Samples ##
#############################################################

# get correct classified samples of test set of first cell 
all_correctSamples.l <- c(correctSamples.l[[1]], correctSamples.l[[2]], correctSamples.l[[3]])
positions.df <- rbind(rbind(all.positMatrix[[1]], all.positMatrix[[2]]), all.positMatrix[[3]])
cor_positons.df <- positions.df[all_correctSamples.l,]

# correctly classified samples which changes ASE
changingASE.m <- rbind(rbind(changingASE.l[[1]], changingASE.l[[2]]), changingASE.l[[3]])
correct_chASE.l <- rownames(cor_positons.df)[which(rownames(cor_positons.df) %in% rownames(changingASE.m))]
pos_correct_chASE.l <- cor_positons.df[correct_chASE.l,]
duplets <- as.data.frame(pos_correct_chASE.l[which(duplicated(pos_correct_chASE.l[,c(1:4)])),])
buff <- row.match(as.data.frame(pos_correct_chASE.l), duplets)
duplets <- pos_correct_chASE.l[which(is.na(buff) == FALSE), ]

# seperate correct classified samples for each cell type
namesRowAll <- rownames(cor_positons.df)

#######################################
## Cluster Correct Classified Sample ##
#######################################

# to calculate the AIC for the kmean clustering
kmeansAIC <- function(fit){
  m <- ncol(fit$centers)
  n <- length(fit$cluster)
  k <- nrow(fit$centers)
  D <- fit$tot.withinss
  return(D + 2*m*k)
}

# function to find optimal number of clusters
kopt <- function(distData.df, k.max, opt, cellname, flag){
  # number of clusters
  print(k.max) # Maximal number of clusters
  
  set.seed(123)
  km <- sapply(1:k.max,
               function(k){kmeans(distData.df, centers = k, nstart = 20)})
  
  # get total within-cluster sum of squares
  twss <- numeric(ncol(km))
  for ( i in 1:ncol(km) ){
    twss[i] <- km[,i]$tot.withinss
  }
  
  # get percentage of variance explained
  perc <- numeric(ncol(km))
  for ( i in 1:ncol(km) ){
    perc[i] <- (km[,i]$betweenss / km[,i]$totss) * 100
  }
  
  aic <- apply(km, 2, kmeansAIC)
  
  postscript(file = paste0('~/workspace/Methylation/plots/kmeans_optimalcluster_evalutaion_',
                        cellname,'_',flag,'.eps'), width=1000, height=700, family="Times")
  par(cex = 2.0, mfrow = c(2,2), cex = 1.5)
  
  plot(1:k.max, log(twss),
       type="b", pch = 19, frame = FALSE,
       xlab="Number of Clusters",
       ylab="Total Within Cluster Sum of Squares (log)",
       xlim=c(0,10))
  abline(v = opt, lty =2)
  
  plot(1:k.max, perc,
       type="b", pch = 19, frame = FALSE,
       xlab="Number of Clusters",
       ylab="Percent of Variance Explained",
       ylim = c(0,100),
       xlim=c(0,10))
  text( (opt - 0.8) , (perc[opt] + 10.0), round(perc[opt], digits = 1))
  abline(v = opt, lty =2)
  
  plot(1:k.max, log(aic),
       type="b", pch = 19, frame = FALSE,
       xlab="Number of Clusters",
       ylab="AIC (log)",
       xlim=c(0,10))
  abline(v = opt, lty =2)
  dev.off()
}

colfunc <- function(x){
  ramp <- colorRampPalette(c('blue', 'red'))
  uniquevalues <- unique(x)
  colors <- ramp(length(uniquevalues))
  sorted_uniques <- sort(uniquevalues)
  for ( i in 1:length(sorted_uniques) ){
    x[which(x == sorted_uniques[i])] = colors[i]
  }
  return(x)
}

optc.l <- list(c())

binning <- function(array, bins){
  pos <- 1
  binwidth <- length(array)/bins
  meanarray <- numeric(bins)
  for ( i in 1:bins ){
    buff <- array[c(pos:(pos+(binwidth-1)))]
    meanarray[i] <- mean(buff)
    pos <- pos + binwidth
  }
  return(meanarray)
}

bins <- 500
width <- 1
# get data for exclusions and inclusions
for ( i in 1:3 ){
  cluster_data.df <- all.data.df[[i]]
  cluster_data.df <- as.data.frame(cluster_data.df[correctSamples.l[[i]],])
  
  exc.df <- cluster_data.df[which(cluster_data.df[,1] == 0), -1]
  inc.df <- cluster_data.df[which(cluster_data.df[,1] == 1), -1]
  
  exc.dist.df <- dist(exc.df, method = "euclidean")
  optc <- 3
  kopt(exc.dist.df, 10, optc, celltype.l[i], 'exc_end')
  set.seed(123)
  clus_exc <- kmeans(exc.dist.df, centers = optc, nstart = 20)
  
  inc.dist.df <- dist(inc.df, method = "euclidean")
  optc2 <- 3
  kopt(inc.dist.df, 10, optc2, celltype.l[i], 'inc_end')
  set.seed(123)
  clus_inc <- kmeans(inc.dist.df, centers = optc2, nstart = 20)
  
  # plot PCA with clusters in it
  data.pca <- prcomp(cluster_data.df[,-1])
  
  # paste cluster names to dinstguish between exclusion and inclusion cluster
  clusters.l <- c(paste0(clus_exc$cluster, '_exc'), paste0(clus_inc$cluster, '_inc'))
  
  postscript(file = paste0('plots/pca_',
                        celltype.l[i],'_exclusions_clusters_end.eps'), width=800, height=800, family="Times")
  p1 <- qplot(data.pca$x[,1], data.pca$x[,2], 
              data = cluster_data.df, color = as.character(clusters.l), 
              xlab = 'PC1', ylab = 'PC2', shape = as.character(clusters.l)) + 
    labs(color = 'Cluster') + 
    scale_shape_manual('Cluster', values = c(1:(optc+optc2))) + 
    theme(text = element_text(size=15, family="Times"))
  p2 <- qplot(data.pca$x[,1], data.pca$x[,3], 
              data = cluster_data.df, color = as.character(clusters.l), 
              xlab = 'PC1', ylab = 'PC3', shape = as.character(clusters.l)) +
    labs(color = 'Cluster') + 
    scale_shape_manual('Cluster', values = c(1:(optc+optc2))) +
    theme(text = element_text(size=15, family="Times"))
  p3 <- qplot(data.pca$x[,2], data.pca$x[,3],
              data = cluster_data.df, color = as.character(clusters.l), 
              xlab = 'PC2', ylab = 'PC3', shape = as.character(clusters.l)) +
    labs(color = 'Cluster') + 
    scale_shape_manual('Cluster', values = c(1:(optc+optc2))) +
    theme(text = element_text(size=15, family="Times"))
  p4 <- qplot(data.pca$x[,1], data.pca$x[,4], 
              data = cluster_data.df, color = as.character(clusters.l), 
              xlab = 'PC1', ylab = 'PC4', shape = as.character(clusters.l)) +
    labs(color = 'Cluster') + 
    scale_shape_manual('Cluster', values = c(1:(optc+optc2))) +
    theme(text = element_text(size=15, family="Times"))
  multiplot(plotlist = list(p1,p2,p3,p4), cols = 2)
  dev.off()
  
  postscript(file = paste0('plots/clustering_',
                        celltype.l[i],'_exclusions_cluster_up_introns.eps'), width=1000, height=1300, family="Times")
  par( mar = c(4,4,2,2), mfrow = c( 3, 1 ), cex = 2.0 )
  for ( j in 1:optc ){
    # sample the specified number of profiles
    c <- rownames(exc.df)[which(clus_exc$cluster == j)]
    cad <- c_all.data.df[[i]]
    c <- cad[c,]
    c <- apply(c, 2, mean)
    
    p1 <- c(paste0('up', c(1:500)))
    c1 <- c[p1] + 1
    c1 <- binning(c1, bins)
    c1 <- c1[-which(is.na(c1))]
    
    colors_values <- colfunc(c1)

    barplot(c1, axes = TRUE, horiz = FALSE, 
            xlab = '', ylab = 'Methylation Level' ,main = paste('Cluster', j), 
            names.arg = '', space = 0, width = width, col = colors_values, 
            border = colors_values,
            xlim = c(0.0, 500) )
    axis(side = 1, at=c(0,100,200,300,400, 490), labels = c(1,100,200,300,400,490))
    title(xlab = 'Relative Position (bp)', line = 2)
  }
  dev.off()
  
  postscript(file = paste0('plots/clustering_',
                        celltype.l[i],'_exclusions_cluster_do_introns.eps'), width=1000, height=1300, family="Times")
  par( mar = c(4,4,2,2), mfrow = c( 3, 1 ), cex = 2.0)
  for ( j in 1:optc ){
    # sample the specified number of profiles
    c <- rownames(exc.df)[which(clus_exc$cluster == j)]
    cad <- c_all.data.df[[i]]
    c <- cad[c,]
    c <- apply(c, 2, mean)
    
    p1 <- c(paste0('do', c(1:500)))
    c1 <- c[p1] + 1 
    c1 <- binning(c1, bins)
    c1 <- c1[-which(is.na(c1))]
    
    colors_values <- colfunc(c1)
   
    barplot(c1, axes = TRUE, horiz = FALSE, 
            xlab = '', ylab = 'Methylation Level' ,main = paste('Cluster', j), 
            names.arg = '', space = 0, width = width, col = colors_values, 
            border = colors_values,
            xlim = c(0.0, 500) )
    axis(side = 1, at=c(0,100,200,300,400, 490), labels = c(11,100,200,300,400,500))
    title(xlab = 'Relative Position (bp)', line = 2)
  }
  dev.off()
  
  postscript(file = paste0('plots/clustering_',
                        celltype.l[i],'_inclusions_cluster_up_introns.eps'), width=1000, height=1300, family="Times")
  par( mar = c(4,4,2,2), mfrow = c( 3, 1 ), cex = 2.0)
  for ( j in 1:optc2 ){
    # sample the specified number of profiles
    c <- rownames(inc.df)[which(clus_inc$cluster == j)]
    cad <- c_all.data.df[[i]]
    c <- cad[c,]
    c <- apply(c, 2, mean)
    
    p1 <- c(paste0('up', c(1:500)))
    c1 <- c[p1] + 1
    c1 <- binning(c1, bins)
    c1 <- c1[-which(is.na(c1))]
    
    colors_values <- colfunc(c1)
    
    barplot(c1, axes = TRUE, horiz = FALSE, 
            xlab = '', ylab = 'Methylation Level' ,main = paste('Cluster', j), 
            names.arg = '', space = 0, width = width, col = colors_values, 
            border = colors_values,
            xlim = c(0.0, 500) )
    axis(side = 1, at=c(0,100,200,300,400, 490), labels = c(1,100,200,300,400,490))
    title(xlab = 'Relative Position (bp)', line = 2)
  }
  dev.off()
  
  
  postscript(file = paste0('plots/clustering_',
                        celltype.l[i],'_inclusions_cluster_do_introns.eps'), width=1000, height=1300, family="Times")
  par( mar = c(4,4,2,2), mfrow = c( 3, 1 ), cex = 2.0)
  for ( j in 1:optc2 ){
    # sample the specified number of profiles
    c <- rownames(inc.df)[which(clus_inc$cluster == j)]
    cad <- c_all.data.df[[i]]
    c <- cad[c,]
    c <- apply(c, 2, mean)
    
    p1 <- c(paste0('do', c(1:500)))
    c1 <- c[p1] + 1
    c1 <- binning(c1, bins)
    c1 <- c1[-which(is.na(c1))]
    
    colors_values <- colfunc(c1)
    
    barplot(c1, axes = TRUE, horiz = FALSE, 
            xlab = '', ylab = 'Methylation Level' ,main = paste('Cluster', j), 
            names.arg = '', space = 0, width = width, col = colors_values, 
            border = colors_values,
            xlim = c(0.0, 500) )
    axis(side = 1, at=c(0,100,200,300,400, 490), labels = c(11,100,200,300,400,500))
    title(xlab = 'Relative Position (bp)', line = 2)
  }
  dev.off()
  
  
  postscript(file = paste0('plots/clustering_',
                        celltype.l[i],'_exclusions_cluster_exons.eps'), width=1000, height=1300, family="Times")
  par( mar = c(4,4,2,2), mfrow = c( 3, 1 ), cex = 2.0 )
  for ( j in 1:optc ){
    # sample the specified number of profiles
    c <- rownames(exc.df)[which(clus_exc$cluster == j)]
    cad <- c_all.data.df[[i]]
    c <- cad[c,]
    c <- apply(c, 2, mean)
    
    p2 <- c(paste0('ex', c(2:49)))
    c2 <- c[p2]
    
    colors_values <- colfunc(c2)
    
    barplot(c2, axes = TRUE, horiz = FALSE, 
            xlab = '', ylab = 'Methylation Level' ,main = paste('Cluster', j), 
            names.arg = '', space = 0, width = 10, col = colors_values,
            border = colors_values,
            xlim = c(0.0, 500) )
    axis(side = 1, at=c(0,100,200,300,400, 480), labels = c(2,10,20,30,40,49))
    title(xlab = 'Relative Position (bin)', line = 2)
  }
  dev.off()
  
  postscript(file = paste0('plots/clustering_',
                        celltype.l[i],'_inclusions_cluster_exons.eps'), width=1000, height=1300, family="Times")
  par( mar = c(4,4,2,2), mfrow = c( 3, 1 ), cex = 2.0 )
  for ( j in 1:optc2 ){
    # sample the specified number of profiles
    c <- rownames(inc.df)[which(clus_inc$cluster == j)]
    cad <- c_all.data.df[[i]]
    c <- cad[c,]
    c <- apply(c, 2, mean)
    
    p2 <- c(paste0('ex', c(2:49)))
    c2 <- c[p2]
    
    colors_values <- colfunc(c2)
    
    barplot(c2, axes = TRUE, horiz = FALSE, 
            xlab = '', ylab = 'Methylation Level' ,main = paste('Cluster', j), 
            names.arg = '', space = 0, width = 10, col = colors_values, 
            border = colors_values,
            xlim = c(0.0, 500) )
    axis(side = 1, at=c(0,100,200,300,400, 480), labels = c(2,10,20,30,40,49))
    title(xlab = 'Relative Position (bin)', line = 2)
  }
  dev.off()
  
}

###########################################
## PCA for correct and incorrect Samples ##
###########################################

PCA <- TRUE

if ( PCA ) {

  # incorrect samples
  incorrect.m <- NULL
  for ( i in 1:3 ){
    set <- all.data.df[[i]]
    if ( i == 1 ){
      incorrect.m <- set[incorrect.l[[i]],]
    } else {
      incorrect.m <- rbind(incorrect.m, set[incorrect.l[[i]],])
    }
  }
  
  # correct samples
  correct.m <- NULL
  for ( i in 1:3 ){
    set <- all.data.df[[i]]
    if ( i == 1 ){
      correct.m <- set[correctSamples.l[[i]],]
    } else {
      correct.m <- rbind(correct.m, set[correctSamples.l[[i]],])
    }
  }
  
  pick <- c(paste0('p', c(1:49)))
  put_inPCA.df <- as.data.frame(rbind(incorrect.m[,pick], correct.m[,pick]))
  data.pca <- prcomp(put_inPCA.df)
  labels.l <- c(rep('Incorrect', nrow(incorrect.m)), rep('Correct', nrow(correct.m)))
  
  postscript(file = paste0('plots/pca_incor_cor_intron.eps'), width=800, height=800, family="Times")
  p1 <- qplot(data.pca$x[,1], data.pca$x[,2], 
              data = put_inPCA.df, color = as.character(labels.l), 
              xlab = 'PC1', ylab = 'PC2', shape = as.character(labels.l)) + 
    labs(color = 'Cluster') + 
    scale_shape_manual('Cluster', values = c(1:2)) + 
    theme(text = element_text(size=15, family="Times"))
  p2 <- qplot(data.pca$x[,1], data.pca$x[,3], 
              data = put_inPCA.df, color = as.character(labels.l), 
              xlab = 'PC1', ylab = 'PC3', shape = as.character(labels.l)) +
    labs(color = 'Cluster') + 
    scale_shape_manual('Cluster', values = c(1:2)) +
    theme(text = element_text(size=15, family="Times"))
  p3 <- qplot(data.pca$x[,2], data.pca$x[,3],
              data = put_inPCA.df, color = as.character(labels.l), 
              xlab = 'PC2', ylab = 'PC3', shape = as.character(labels.l)) +
    labs(color = 'Cluster') + 
    scale_shape_manual('Cluster', values = c(1:2)) +
    theme(text = element_text(size=15, family="Times"))
  p4 <- qplot(data.pca$x[,1], data.pca$x[,4], 
              data = put_inPCA.df, color = as.character(labels.l), 
              xlab = 'PC1', ylab = 'PC4', shape = as.character(labels.l)) +
    labs(color = 'Cluster') + 
    scale_shape_manual('Cluster', values = c(1:2)) +
    theme(text = element_text(size=15, family="Times"))
  multiplot(plotlist = list(p1,p2,p3,p4), cols = 2)
  dev.off()
  
  
  pick <- c(paste0('ex', c(2:49)))
  put_inPCA.df <- as.data.frame(rbind(incorrect.m[,pick], correct.m[,pick]))
  data.pca <- prcomp(put_inPCA.df)
  labels.l <- c(rep('Incorrect', nrow(incorrect.m)), rep('Correct', nrow(correct.m)))
  
  postscript(file = paste0('plots/pca_incor_cor_exons.eps'), width=800, height=800, family="Times")
  p1 <- qplot(data.pca$x[,1], data.pca$x[,2], 
              data = put_inPCA.df, color = as.character(labels.l), 
              xlab = 'PC1', ylab = 'PC2', shape = as.character(labels.l)) + 
    labs(color = 'Cluster') + 
    scale_shape_manual('Cluster', values = c(1:2)) + 
    theme(text = element_text(size=15, family="Times"))
  p2 <- qplot(data.pca$x[,1], data.pca$x[,3], 
              data = put_inPCA.df, color = as.character(labels.l), 
              xlab = 'PC1', ylab = 'PC3', shape = as.character(labels.l)) +
    labs(color = 'Cluster') + 
    scale_shape_manual('Cluster', values = c(1:2)) +
    theme(text = element_text(size=15, family="Times"))
  p3 <- qplot(data.pca$x[,2], data.pca$x[,3],
              data = put_inPCA.df, color = as.character(labels.l), 
              xlab = 'PC2', ylab = 'PC3', shape = as.character(labels.l)) +
    labs(color = 'Cluster') + 
    scale_shape_manual('Cluster', values = c(1:2)) +
    theme(text = element_text(size=15, family="Times"))
  p4 <- qplot(data.pca$x[,1], data.pca$x[,4], 
              data = put_inPCA.df, color = as.character(labels.l), 
              xlab = 'PC1', ylab = 'PC4', shape = as.character(labels.l)) +
    labs(color = 'Cluster') + 
    scale_shape_manual('Cluster', values = c(1:2)) +
    theme(text = element_text(size=15, family="Times"))
  multiplot(plotlist = list(p1,p2,p3,p4), cols = 2)
  dev.off()
  
  
  pick <- c(paste0('p', c(1:49)))
  put_inPCA.df <- as.data.frame(correct.m[,pick])
  data.pca <- prcomp(put_inPCA.df)
  labels.l <- rep('', nrow(correct.m))
  labels.l[grep('ex', rownames(correct.m))] <- 'Exclusions'
  labels.l[grep('in', rownames(correct.m))] <- 'Inclusions'
  
  postscript(file = paste0('plots/pca_cor_intron.eps'), width=800, height=800, family="Times")
  p1 <- qplot(data.pca$x[,1], data.pca$x[,2], 
              data = put_inPCA.df, color = as.character(labels.l), 
              xlab = 'PC1', ylab = 'PC2', shape = as.character(labels.l)) + 
    labs(color = 'Cluster') + 
    scale_shape_manual('Cluster', values = c(1:2)) + 
    theme(text = element_text(size=15, family="Times"))
  p2 <- qplot(data.pca$x[,1], data.pca$x[,3], 
              data = put_inPCA.df, color = as.character(labels.l), 
              xlab = 'PC1', ylab = 'PC3', shape = as.character(labels.l)) +
    labs(color = 'Cluster') + 
    scale_shape_manual('Cluster', values = c(1:2)) +
    theme(text = element_text(size=15, family="Times"))
  p3 <- qplot(data.pca$x[,2], data.pca$x[,3],
              data = put_inPCA.df, color = as.character(labels.l), 
              xlab = 'PC2', ylab = 'PC3', shape = as.character(labels.l)) +
    labs(color = 'Cluster') + 
    scale_shape_manual('Cluster', values = c(1:2)) +
    theme(text = element_text(size=15, family="Times"))
  p4 <- qplot(data.pca$x[,1], data.pca$x[,4], 
              data = put_inPCA.df, color = as.character(labels.l), 
              xlab = 'PC1', ylab = 'PC4', shape = as.character(labels.l)) +
    labs(color = 'Cluster') + 
    scale_shape_manual('Cluster', values = c(1:2)) +
    theme(text = element_text(size=15, family="Times"))
  multiplot(plotlist = list(p1,p2,p3,p4), cols = 2)
  dev.off()
  
  
  
  pick <- c(paste0('ex', c(2:49)))
  put_inPCA.df <- as.data.frame(correct.m[,pick])
  data.pca <- prcomp(put_inPCA.df)
  labels.l <- rep('', nrow(correct.m))
  labels.l[grep('ex', rownames(correct.m))] <- 'Exclusions'
  labels.l[grep('in', rownames(correct.m))] <- 'Inclusions'
  
  postscript(file = paste0('plots/pca_cor_exons.eps'), width=800, height=800, family="Times")
  p1 <- qplot(data.pca$x[,1], data.pca$x[,2], 
              data = put_inPCA.df, color = as.character(labels.l), 
              xlab = 'PC1', ylab = 'PC2', shape = as.character(labels.l)) + 
    labs(color = 'Cluster') + 
    scale_shape_manual('Cluster', values = c(1:2)) + 
    theme(text = element_text(size=15, family="Times"))
  p2 <- qplot(data.pca$x[,1], data.pca$x[,3], 
              data = put_inPCA.df, color = as.character(labels.l), 
              xlab = 'PC1', ylab = 'PC3', shape = as.character(labels.l)) +
    labs(color = 'Cluster') + 
    scale_shape_manual('Cluster', values = c(1:2)) +
    theme(text = element_text(size=15, family="Times"))
  p3 <- qplot(data.pca$x[,2], data.pca$x[,3],
              data = put_inPCA.df, color = as.character(labels.l), 
              xlab = 'PC2', ylab = 'PC3', shape = as.character(labels.l)) +
    labs(color = 'Cluster') + 
    scale_shape_manual('Cluster', values = c(1:2)) +
    theme(text = element_text(size=15, family="Times"))
  p4 <- qplot(data.pca$x[,1], data.pca$x[,4], 
              data = put_inPCA.df, color = as.character(labels.l), 
              xlab = 'PC1', ylab = 'PC4', shape = as.character(labels.l)) +
    labs(color = 'Cluster') + 
    scale_shape_manual('Cluster', values = c(1:2)) +
    theme(text = element_text(size=15, family="Times"))
  multiplot(plotlist = list(p1,p2,p3,p4), cols = 2)
  dev.off()
  
}

#############################
## Hierarchical Clustering ##
#############################

my_palette <- colorRampPalette(c('white', 'red', 'black'))(n = 500)

pick <- c(paste0('p', c(1:49)))
all.m <- correct.m[, pick]
sepExcl <- numeric(nrow(correct.m))
sepExcl[which(correct.m[,1] == 0)] <- 'blue'
sepExcl[which(correct.m[,1] == 1)] <- 'red'

sums <- apply(all.m, 1, sum)
all.m <- all.m[which(sums != -49),]
sepExcl <- sepExcl[which(sums != -49)]

postscript(file = paste0('plots/heatmap_cor_up_introns_all.eps'), width=800, height=800, family="Times")
par(cex.axis = 1.5, cex.main = 1.5, cex.lab = 1.5)
heatmap.2(all.m, Rowv = TRUE, Colv = FALSE,  dendrogram = 'row',
          distfun = function(d) dist(d, method = 'euclidean'), 
          hclustfun = function(d) hclust(d, method = 'average'),
          density.info = 'none', trace = 'none', col = my_palette, labRow = FALSE, 
          ylab = 'Samples', xlab = 'Features', cexCol = 1.2, margins = c(5,5), RowSideColors = sepExcl)
dev.off()


pick <- c(paste0('d', c(2:50)))
all.m <- correct.m[, pick]
sepExcl <- numeric(nrow(correct.m))
sepExcl[which(correct.m[,1] == 0)] <- 'blue'
sepExcl[which(correct.m[,1] == 1)] <- 'red'

sums <- apply(all.m, 1, sum)
all.m <- all.m[which(sums != -49),]
sepExcl <- sepExcl[which(sums != -49)]

postscript(file = paste0('plots/heatmap_cor_do_introns_all.eps'), width=800, height=800, family="Times")
par(cex.axis = 1.5, cex.main = 1.5, cex.lab = 1.5)
heatmap.2(all.m, Rowv = TRUE, Colv = FALSE,  dendrogram = 'row',
          distfun = function(d) dist(d, method = 'euclidean'), 
          hclustfun = function(d) hclust(d, method = 'average'),
          density.info = 'none', trace = 'none', col = my_palette, labRow = FALSE, 
          ylab = 'Samples', xlab = 'Features', cexCol = 1.2, margins = c(5,5), RowSideColors = sepExcl)
dev.off()

pick <- c(paste0('ex', c(2:49)))
all.m <- correct.m[, pick]
sepExcl <- numeric(nrow(correct.m))
sepExcl[which(correct.m[,1] == 0)] <- 'blue'
sepExcl[which(correct.m[,1] == 1)] <- 'red'

sums <- apply(all.m, 1, sum)
all.m <- all.m[which(sums != 0),]
sepExcl <- sepExcl[which(sums != 0)]

postscript(file = paste0('plots/heatmap_cor_exons_all.eps'), width=800, height=800, family="Times")
par(cex.axis = 1.5, cex.main = 1.5, cex.lab = 1.5)
heatmap.2(all.m, Rowv = TRUE, Colv = FALSE,  dendrogram = 'row',
          distfun = function(d) dist(d, method = 'euclidean'), 
          hclustfun = function(d) hclust(d, method = 'average'),
          density.info = 'none', trace = 'none', col = my_palette, labRow = FALSE, 
          ylab = 'Samples', xlab = 'Features', cexCol = 1.0, margins = c(5,5), RowSideColors = sepExcl)
dev.off()

picks <- grep('ex', rownames(duplets))
picks.l <- NULL
for ( i in 1:length(picks) ){
  one <- duplets[picks[i],]
  two <- row.match(as.data.frame(duplets[,c(1:4)]), matrix(one[c(1:4)], nrow = 1, ncol = 4))
  two <- which(!is.na(two))
  n_one <- rownames(duplets)[picks[i]]
  two <- rownames(duplets)[two[which(!rownames(duplets)[two] %in% n_one)]]
  two <- two[grep('in', two)]
  picks.l <- append(picks.l, c(n_one, two[1]))
}

# incorrect samples
incorrect.m <- NULL
for ( i in 1:3 ){
  set <- c_all.data.df[[i]]
  if ( i == 1 ){
    incorrect.m <- set[incorrect.l[[i]],]
  } else {
    incorrect.m <- rbind(incorrect.m, set[incorrect.l[[i]],])
  }
}

# correct samples
correct.m <- NULL
for ( i in 1:3 ){
  set <- c_all.data.df[[i]]
  if ( i == 1 ){
    correct.m <- set[correctSamples.l[[i]],]
  } else {
    correct.m <- rbind(correct.m, set[correctSamples.l[[i]],])
  }
}

my_palette <- colorRampPalette(c('white', 'red', 'black'))(n = 3)
pick <- c(paste0('up', c(1:490)))
all.m <- correct.m[picks.l, pick]
sepExcl <- numeric(nrow(all.m))
sepExcl[grep('ex',rownames(all.m))] <- 'blue'
sepExcl[grep('in',rownames(all.m))] <- 'red'

postscript(file = paste0('plots/heatmap_ASEall_up_introns.eps'), width=800, height=800, family="Times")
par(cex.axis = 1.5, cex.main = 1.5, cex.lab = 1.5)
heatmap.2(all.m, Rowv = FALSE, Colv = FALSE,  dendrogram = 'none',
          density.info = 'none', trace = 'none', col = my_palette,
          ylab = 'Samples', xlab = 'Features', cexRow = 1.0, cexCol = 0.1, margins = c(5,10),
          RowSideColors = sepExcl, 
          key.xtickfun=function() {
            breaks <- parent.frame()$breaks
            return(list(
              at=parent.frame()$scale01(c(-1, 0.0, 1.0)),
              labels=c(as.character(-1), as.character(0.0), as.character(1))
            ))
          })
dev.off()


pick <- c(paste0('do', c(11:500)))
all.m <- correct.m[picks.l, pick]
sepExcl <- numeric(nrow(all.m))
sepExcl[grep('ex',rownames(all.m))] <- 'blue'
sepExcl[grep('in',rownames(all.m))] <- 'red'

postscript(file = paste0('plots/heatmap_ASEall_do_introns.eps'), width=800, height=800, family="Times")
par(cex.axis = 1.5, cex.main = 1.5, cex.lab = 1.5)
heatmap.2(all.m, Rowv = FALSE, Colv = FALSE,  dendrogram = 'none',
          density.info = 'none', trace = 'none', col = my_palette,
          ylab = 'Samples', xlab = 'Features', cexRow = 1.0, cexCol = 0.1, margins = c(5,10),
          RowSideColors = sepExcl, 
          key.xtickfun=function() {
            breaks <- parent.frame()$breaks
            return(list(
              at=parent.frame()$scale01(c(-1, 0.0, 1.0)),
              labels=c(as.character(-1), as.character(0.0), as.character(1))
            ))
          })
dev.off()


my_palette <- colorRampPalette(c('white', 'red', 'black'))(n = 500)
pick <- c(paste0('ex', c(2:49)))
all.m <- correct.m[picks.l, pick]
sepExcl <- numeric(nrow(all.m))
sepExcl[grep('ex',rownames(all.m))] <- 'blue'
sepExcl[grep('in',rownames(all.m))] <- 'red'

postscript(file = paste0('plots/heatmap_ASEall_exons.eps'), width=800, height=800, family="Times")
par(cex.axis = 1.5, cex.main = 1.5, cex.lab = 1.5)
heatmap.2(all.m, Rowv = FALSE, Colv = FALSE,  dendrogram = 'none',
          density.info = 'none', trace = 'none', col = my_palette,
          ylab = 'Samples', xlab = 'Features', cexCol = 1.0, margins = c(5,10),
          RowSideColors = sepExcl)
dev.off()