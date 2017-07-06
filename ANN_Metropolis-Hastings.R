
args = commandArgs(trailingOnly = TRUE)
params <- as.numeric(args)         # for setting default parameters or random parameters
                       		   # have to be an array defining all parameters
                       		   # (units, layer, batch size, dropout, 
                       		   #   rate, momentum, regularization, )

celltype <- 'IMR90'    # the celltype of the data (IMR90, Gm12878, H1hesc)


print(as.vector(params))
print(params[-1])

#################################
## Load and install libraries  ##
#################################

# for the plots 
library(ggplot2)
library(Rmisc)
library(splines)
library(MASS)
library(mxnet)
library(doParallel) 
library(parallel)
library(data.table)

################
## Load data  ##
################

data.df <- as.data.frame(fread(paste0('data/network_data_',
                                      paste0(celltype, collapse = '_'),'.tsv'), sep =  "\t", header = TRUE))

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

# change to numeric type 
buff.df <- data.df
data.df <- sapply(data.df, as.numeric)
rownames(data.df) <- rownames(buff.df)

# testing 
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

features <- data.df[,paste0('up',c(1:500))]
binned_Up <- binfeatures(features, 'u')

# pick cols
pickedcols <- c('target', paste0('ex',c(2:49)))
picked_data.df <- data.df[,pickedcols]
colnames(picked_data.df) 

picked_data.df <- cbind(picked_data.df, binned_Up[,-50])

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

print(paste('[NOTE] number of samples in the whole feature dataset =', nrow(data.df)))

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

##################################
## first test run other models  ##
##################################

downsampled_data.df = downSample(picked_data.df)

# function(data.df, pTrain, pValid, loadTest, setSeed, setBagging)
datasets.l <- splitData(downsampled_data.df, 0.6, 0.1, TRUE, TRUE, FALSE)
train.m <- data.matrix(datasets.l[[1]])
valid.m <- data.matrix(datasets.l[[2]])
test.m <- data.matrix(datasets.l[[3]])

combination.train.valid.m <- rbind(train.m, valid.m, test.m)

###################################
## Metroplois Hasting Algorithm  ##
###################################

# Prior distribution
prior <- function(p){
  first_prior = dunif(p[1], min=1, max=2000, log = T)
  second_prior = dunif(p[2], min=0.0001, max=1.0, log = T)
  third_prior = dunif(p[3], min=0.0, max=0.9, log = T)
  fourth_prior = dunif(p[4], min=0.00001, max=1.0, log = T)
  fifth_prior = dunif(p[5], min=0.001, max=0.95, log = T)
  return(first_prior + second_prior + third_prior + fourth_prior + fifth_prior)
}

likelihoodTarget <- function(arch,p) {
  
  mx.set.seed(Sys.time())
  
  model <- mx.mlp(data=train.m[,-1], label=train.m[,1], hidden_node=arch, 
                  out_node=2, num.round=10, array.batch.size=as.integer(p[1]), 
                  learning.rate=p[2],
                  initializer = mx.init.uniform(0.01), optimizer = "sgd", dropout=p[3],  
                  activation = "tanh", out_activation="softmax", eval.metric=mx.metric.accuracy,
                  eval.data=list(data=valid.m[,-1], label=valid.m[,1]), wd = p[4], 
                  momentum=p[5])
  
  # make prediction
  pred <- predict(model, valid.m[,-1])

  c0 <- which(as.data.frame(valid.m[,1]) == 0)
  c1 <- which(as.data.frame(valid.m[,1]) == 1)
  
  t.pred <- t(pred)
  
  p0 <- t.pred[c0]
  p1 <- t.pred[c1]

  return( sum( c(log(p0), log(p1)) ) )
}

posterior <- function(arch,param){
  return (likelihoodTarget(arch,param) + prior(param))
}

rwmetro <- function(N, arch, p, var, burnin=0) {
  samples <- p
  
  for (i in 2:(burnin+N))
  {   
    # proposal function
    prop <- abs(mvrnorm(n = 1, p, var*diag(length(p)) ))
    
    max.l <- c(2000, 1.0, 0.9, 1.0, 0.95)
    
    for ( i in 1:length(prop) ){
      if ( prop[i] > max.l[i] ){
        prop[i] <- runif(1, min=0.0, max = max.l[i])
      }
    }
    
    print(prop)
    # I set min for the proposal to take a look what comes out
    # The min of the proposal is either 1 if 
    # posterior(arch,prop) > posterior(arch,p)
    # of smaller 1 if 
    # posterior(arch,prop) < posterior(arch,p)
    proposal <- min(1, exp(posterior(arch,prop)-posterior(arch,p)))
    print(proposal)
    
    if ( !is.na(proposal) ) {
      # acceptance function 
      if (runif(1) < min(1, proposal )) {
        p <- prop
      }
    }
    
    # bind parameter setting to the matrix
    samples <- rbind(samples, c(as.integer(p[1]), p[2], p[3], p[4], p[5] ) )
    
    if (i %% 100 == 0 ){print(i)}
    
  }
  samples[(burnin+1):(N+burnin),]
}

# keep units and lazers constant 

print(params[1])
print(params[-1])

# params = c(units, batch size, rate, dropout, regularization, momentum)
MH.results.df <- rwmetro(1000, params[1],
                         params[-1], c(10.0, 0.1, 0.1, 0.1, 0.1), 100)

acceptance = 1-mean(as.numeric(duplicated(MH.results.df)))

print(MH.results.df)
print(acceptance)

log.file = 'data/ANN_MH_log.txt'

if ( file.exists(log.file) == FALSE ){
  file.create(log.file)
} 

sink(file=log.file, append = TRUE)
print(MH.results.df[nrow(MH.results.df),])
print(acceptance)
sink()
closeAllConnections()



