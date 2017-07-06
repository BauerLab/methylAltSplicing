
intronLength <- 500    # same number as intron length in the data 
degradation <- 0.1     # for the methylation density profile
filterB <- FALSE       # to apply the filter based on the quantiles 
quart <- 3             # which quartile you want to filter on (2=first,3=median,4=third)
noise1 <- FALSE        # to add noise to the data before you filtered based on the quantiles 
noise2 <- FALSE        # to add noise to the data after you filtered based on the quantiles 
densityData <- FALSE   # to use the density methylation pattern data set this TRUE
regimes <- FALSE       # defines new labels for the samples
celltype <- c('IMR90')    # the celltype of the data (IMR90, Gm12878, H1hesc)

# add additional features 
otherfeatures.b <- c(TRUE, # methylation content of exons and introns
                     FALSE, # GC content of exons and introns
                     TRUE, # exon methylation pattern
                     FALSE, # cytosine content of upstream introns, exon and downstream introns
                     FALSE, # runstest p-value for upstream and downstream introns
                     TRUE, # include 10 bp around splice sites
                     FALSE, # remove the whole methylation profile for exons and introns
                     TRUE, # mCpG/CpG ratio for exons and introns 
                     TRUE, # add methyaltion density profile of the introns
                     TRUE, # add mCpG/CpG ratio profile of exons 
                     FALSE)  # remove downstream intron pattern

#################################
## Load and install libraries  ##
#################################

# for the plots 
library(ggplot2)
library(Rmisc)
library(splines)
library(MASS)

alldata.df <- NULL
for ( c in 1:length(celltype) ) {

  ################
  ## Load data  ##
  ################
  
  # load exclusions
  features.exclusion.df <- read.table(file=paste0('data/features/features.FALSE.FALSE.0.',intronLength,'.',celltype[c],'.tsv'), 
                                     header=TRUE, colClasses='character')
  row.names(features.exclusion.df) <- NULL
  
  # load inclusions
  features.inclusion.df <- read.delim(file=paste0('data/features/features.TRUE.FALSE.0.',intronLength,'.',celltype[c],'.tsv'), 
                                     header=TRUE, colClasses='character')
  row.names(features.inclusion.df) <- NULL
  
  print(paste0('[NOTE] number of exclusion = ', nrow(features.exclusion.df)))
  print(paste0('[NOTE] number of inclusion = ', nrow(features.inclusion.df)))
  
  ###################
  ## Prepare data  ##
  ###################
  
  # function to create matrix out of the bitwise vector
  converBitvector <- function(features, matrix, splitchar) {
    list <- lapply(features, function(x) return( as.numeric( unlist(strsplit(x, split=splitchar)) ) ))
    for(i in 1:length(list) ){
      matrix[i,] <- list[[i]]
    }
    return(matrix)
  }
  
  # prepare inclusions features
  # splite updtream and downstream bitwise vectors into columns
  up.inc <- matrix(nrow=nrow(features.inclusion.df), ncol=intronLength)
  up.inc <- converBitvector(features.inclusion.df[,6], up.inc, '\\|')
  
  do.inc <- matrix(nrow=nrow(features.inclusion.df), ncol=intronLength)
  do.inc <- converBitvector(features.inclusion.df[,7], do.inc, '\\|')
  
  # create matrix out of the target and bitwise vector for upstream and downstream intron
  pre_features.inclusions.df <- matrix(nrow=nrow(features.inclusion.df), ncol=((2*intronLength)+1))
  
  colnames(pre_features.inclusions.df) <- c('target', 
                                           paste0('up',as.character(1:intronLength)), 
                                           paste0('do',as.character(1:intronLength)))
  
  pre_features.inclusions.df[,1] <- as.numeric(features.inclusion.df[,4])
  pre_features.inclusions.df[,2:(intronLength+1)] <- up.inc[,1:intronLength]
  pre_features.inclusions.df[,(intronLength+2):((2*intronLength)+1)] <- do.inc[,1:intronLength]
  
  # prepare exclusion features
  # splite updtream and downstream bitwise vectors into columns
  up.ex <- matrix(nrow=nrow(features.exclusion.df), ncol=intronLength)
  up.ex <- converBitvector(features.exclusion.df[,6], up.ex, '\\|')
  
  do.ex <- matrix(nrow=nrow(features.exclusion.df), ncol=intronLength)
  do.ex <- converBitvector(features.exclusion.df[,7], do.ex, '\\|')
  
  # create matrix out of the target and bitwise vector for upstream and downstream intron
  pre_features.exclusions.df <- matrix(nrow=nrow(features.exclusion.df), ncol=((2*intronLength)+1))
  
  colnames(pre_features.exclusions.df) <- c('target', 
                                     paste0('up',as.character(1:intronLength)), 
                                     paste0('do',as.character(1:intronLength)))
  
  pre_features.exclusions.df[,1] <- as.numeric(features.exclusion.df[,4])
  pre_features.exclusions.df[,2:(intronLength+1)] <- up.ex[,1:intronLength]
  pre_features.exclusions.df[,(intronLength+2):((2*intronLength)+1)] <- do.ex[,1:intronLength]
  
  # join exclusions and inclusion features together 
  pre_features.df <- rbind(pre_features.exclusions.df, pre_features.inclusions.df)
  
  # set row names of feature matrix 
  rownames(pre_features.df) <- c(paste0('ex',features.exclusion.df[,1]), paste0('in',features.inclusion.df[,1]))
  
  # matrix which holds the PSI, MAP_PSI and average methyaltion rate of the exons  
  PSI.df <- matrix(nrow=(nrow(features.exclusion.df)+nrow(features.inclusion.df)), ncol=17)
  PSI.df[,1] <- c(features.exclusion.df[,2], features.inclusion.df[,2])
  PSI.df[,2] <- c(features.exclusion.df[,3], features.inclusion.df[,3])
  PSI.df[,3] <- c(features.exclusion.df[,5], features.inclusion.df[,5])
  PSI.df[,6] <- c(features.exclusion.df[,8], features.inclusion.df[,8])
  PSI.df[,7] <- c(features.exclusion.df[,9], features.inclusion.df[,9])
  PSI.df[,8] <- c(features.exclusion.df[,10], features.inclusion.df[,10])
  PSI.df[,9] <- c(features.exclusion.df[,11], features.inclusion.df[,11])
  PSI.df[,10] <- c(features.exclusion.df[,12], features.inclusion.df[,12])
  PSI.df[,11] <- c(features.exclusion.df[,13], features.inclusion.df[,13])
  PSI.df[,13] <- c(features.exclusion.df[,15], features.inclusion.df[,15])
  PSI.df[,14] <- c(features.exclusion.df[,16], features.inclusion.df[,16])
  PSI.df[,15] <- c(features.exclusion.df[,17], features.inclusion.df[,17])
  PSI.df[,16] <- c(features.exclusion.df[,18], features.inclusion.df[,18])
  PSI.df[,17] <- c(features.exclusion.df[,19], features.inclusion.df[,19])
  colnames(PSI.df) <- c('PSI','MAP_PSI', 'methEx', 'methUp', 'methDo','UpGC',
                       'DoGC','ExGC','UpC','DoC','ExC','M_Labels','UpRun','DoRun',
                       'UpmCpG/Cpg', 'ExmCp/Cpg', 'DomCpg/Cpg')
  rownames(PSI.df) <- c(paste0('ex',features.exclusion.df[,1]), paste0('in',features.inclusion.df[,1]))
  
  # set data
  data.df <- pre_features.df
  
  ################################
  ## Get Exons Feature Vectors  ##
  ################################
  
  # get the number of bins for the exon profile 
  exBins <- length(unlist(strsplit(features.exclusion.df[1,14], split='\\|')))
  
  # convert methylation pattern into seperated features
  ex.inc <- matrix(nrow=nrow(features.inclusion.df), ncol=exBins)
  ex.inc <- converBitvector(features.inclusion.df[,14], ex.inc, '\\|')
  
  ex.exc <- matrix(nrow=nrow(features.exclusion.df), ncol=exBins)
  ex.exc <- converBitvector(features.exclusion.df[,14], ex.exc, '\\|')
  
  # put data together for inclusion and exclusion
  exon.features.df <- ex.exc
  exon.features.df <- rbind(exon.features.df, ex.inc)
  
  colnames(exon.features.df) <- c(paste0('ex', as.character(1:exBins)))
  rownames(exon.features.df) <- c(paste0('ex', features.exclusion.df[,1]), paste0('in',features.inclusion.df[,1]))
  
  # remove -1 for the methylation profile 
  # else the exon profile is distorted towards -1
  # the bins have a value part of the 
  exon.features.df[which(exon.features.df == -1)] = 0.0
  
  #######################################
  ## Get Exons mCpG/CpG Ratio Vectors  ##
  #######################################
  
  # convert methylation pattern into seperated features
  ratio_ex.inc <- matrix(nrow=nrow(features.inclusion.df), ncol=exBins)
  ratio_ex.inc <- converBitvector(features.inclusion.df[,20], ratio_ex.inc, '\\|')
  
  ratio_ex.exc <- matrix(nrow=nrow(features.exclusion.df), ncol=exBins)
  ratio_ex.exc <- converBitvector(features.exclusion.df[,20], ratio_ex.exc, '\\|')
  
  # put data together for inclusion and exclusion
  ratio_exon.features.df <- ratio_ex.exc
  ratio_exon.features.df <- rbind(ratio_exon.features.df, ratio_ex.inc)
  
  colnames(ratio_exon.features.df) <- c(paste0('exR',as.character(1:exBins)))
  rownames(ratio_exon.features.df) <- c(paste0('ex',features.exclusion.df[,1]), paste0('in',features.inclusion.df[,1]))
  
  ###################################
  ## Get Average Methylation Rate  ##
  ###################################
  
  # get target of features 
  sum.features.df <- as.data.frame(pre_features.df[,1])
  
  # change -1 to 0
  changed_pre_features.df <- pre_features.df
  changed_pre_features.df[which(changed_pre_features.df == -1)] <- 0
  
  # get the sum for each row for upstream region 
  sum.features.df <- cbind(sum.features.df, apply(changed_pre_features.df[,2:(intronLength+1)], 1, sum))
  
  # get the sum for each row for downstream region 
  sum.features.df <- cbind(sum.features.df, apply(changed_pre_features.df[,(intronLength+2):((intronLength*2)+1)], 1, sum))
  
  # create df
  colnames(sum.features.df) <- c('target', 'sumUp', 'sumDo')
  rownames(sum.features.df) <- rownames(pre_features.df)
  
  # save the methylation rates of upstream and downstream intron 
  PSI.df[,4] <- sum.features.df[,2] / intronLength
  PSI.df[,5] <- sum.features.df[,3] / intronLength
  
  ########################
  ## Create new Labels  ##
  ########################
  
  # get data from introns
  # get data from introns
  exc.l <- which(data.df[,1] == 0)
  inc.l <- which(data.df[,1] == 1)
  
  exc.ex.mC.q <- quantile(as.numeric(PSI.df[exc.l, 16]))
  lowmeth_exons.exlc <- exc.l[which(as.numeric(PSI.df[exc.l, 16]) <= 0.0)] 
  highmeth_exons.exlc <- exc.l[which(as.numeric(PSI.df[exc.l, 16]) == 1.0)] 
  
  labels <- c()
  if(regimes == F) {
    labels <- c(rep('0',3), rep('1',3))
  } else {
    labels <- c('M0', 'L0', 'H0', 'M1', 'L1', "H1")
  }
  
  PSI.df[exc.l, 12] <- labels[1]
  PSI.df[lowmeth_exons.exlc, 12] <- labels[2]
  PSI.df[highmeth_exons.exlc, 12] <- labels[3]
  
  inc.ex.mC.q <- quantile(as.numeric(PSI.df[inc.l, 16]))
  lowmeth_exons.incl <- inc.l[which(as.numeric(PSI.df[inc.l, 16]) <= 0.0)]
  highmeth_exons.incl <- inc.l[which(as.numeric(PSI.df[inc.l, 16]) == 1.0)]
  
  PSI.df[inc.l, 12] <- labels[4]
  PSI.df[lowmeth_exons.incl, 12] <- labels[5]
  PSI.df[highmeth_exons.incl, 12] <- labels[6]
  
  #############################################
  ## Create pseudocount methylation profile  ##
  #############################################
  
  peakMeth <- function(x) {
    counter <- 0
    y <- numeric(intronLength)
    for(i in 1:intronLength){
      if (x[i] == 1) {
        counter <- counter + 1
      } else {
        counter <- counter - degradation
      }
  
      if (counter < 0) {
        counter <- 0
      }
  
      y[i] <- counter
    }
    return(y)
  }
  
  if ( densityData == TRUE || otherfeatures.b[9] == TRUE ) {
  
    # apply !!!! you have to transpose the matrix again 
    peaks_up.df <- apply(pre_features.df[,2:(intronLength+1)], 1, peakMeth)
    peaks_up.df <- t(peaks_up.df)
    
    peaks_do.df <- apply(pre_features.df[,(intronLength+2):((intronLength*2)+1)], 1, peakMeth)
    peaks_do.df <- t(peaks_do.df)
    
    peaks.df <- pre_features.df
    
    # take the over all methylation rate into account 
    peaks_up.df <- peaks_up.df/sum.features.df[rownames(peaks_up.df),2]
    peaks_do.df <- peaks_do.df/sum.features.df[rownames(peaks_up.df),3]
    
    peaks_up.df[which(is.na(peaks_up.df))] <- 0.0
    peaks_do.df[which(is.na(peaks_do.df))] <- 0.0
    
    # add gaussian noise if true
    noise1 <- FALSE
    if ( noise1 == TRUE ) {
      print('[NOTE] your turned on the noise1')
      noise_peaks_up.df <- peaks_up.df + rnorm(nrow(peaks_up.df) * ncol(peaks_up.df))
      noise_peaks_do.df <- peaks_do.df + rnorm(nrow(peaks_do.df) * ncol(peaks_do.df))
      peaks.df[,2:(intronLength+1)] <- noise_peaks_up.df
      peaks.df[,(intronLength+2):((intronLength*2)+1)] <- noise_peaks_do.df
    } else {
      peaks.df[,2:(intronLength+1)] <- peaks_up.df
      peaks.df[,(intronLength+2):((intronLength*2)+1)] <- peaks_do.df
    }
  
  
    ###############
    ## Filter B  ##
    ###############
  
    if ( filterB == TRUE ) {
      # the function looks at the methylation profile
      # of the upstream and downstream intron and picks a group of
      # up/do-introns based on their maximal peak in the methylation profile
      peaksQuantil <- function(peaks, ase, q){
        ase.l <- which(peaks[,1] == ase)
        max.up.l <- apply(peaks[ase.l,2:(intronLength+1)], 1, max)
        max.do.l <- apply(peaks[ase.l,(intronLength+2):((intronLength*2)+1)], 1, max)
      
        # keep upstream and downstream intron if both are above the quantile q
        keep.l <- ase.l[which(max.up.l > quantile(max.up.l)[[q]] )]
        keep.l <- append(keep.l, ase.l[which(max.do.l > quantile(max.do.l)[[q]] )])
        keep.l <- unique(keep.l)
        return(keep.l)
      }
      
      # reduce peaks to significant introns
      keep0.l <- peaksQuantil(peaks.df,0,quart)
      keep1.l <- peaksQuantil(peaks.df,1,quart)
      
      data.df <- peaks.df[keep0.l,]
      data.df <- rbind(data.df, peaks.df[keep1.l,])
      
      rownames(data.df) <- c(rownames(peaks.df)[keep0.l], rownames(peaks.df)[keep1.l])
      
      # reduce PSI matrix to significant introns
      sign_PSI.df <- PSI.df[keep0.l,]
      sign_PSI.df <- rbind(sign_PSI.df, PSI.df[keep1.l,])
      
      # reduce exon feature matrix to significant introns
      sign_exon.features.df <- exon.features.df[keep0.l,]
      sign_exon.features.df <- rbind(sign_exon.features.df, exon.features.df[keep1.l,])
      
    } else {
      data.df <- peaks.df
    }
      
    if ( noise2 == TRUE ) {
      print('[NOTE] your turned on the noise2')
      # data.df[,-1] = data.df[,-1] + rnorm(nrow(data.df) * (ncol(data.df)-1))
      data.df[,-1] <- jitter(data.df[,-1])
    }
  }
  
  ###############################
  ## Featur selection process  ##
  ###############################
  
  # check which data the user has chosen (methylation pattern or methylation density pattern)
  otherfeatures.df <- PSI.df
  exon.data.df <- exon.features.df
  if ( densityData == FALSE ) {
    print('[NOTE] you use the normal methylation pattern data')
    data.df <- pre_features.df
  } else {
    print('[NOTE] you use the density methylation pattern data')
    if ( filterB == TRUE ) {
      otherfeatures.df <- sign_PSI.df
      exon.data.df <- sign_exon.features.df
    } else {
      data.df <- peaks.df
    }
  }
  # look if the main features and the other features have the same rows 
  if ( nrow(data.df) != nrow(otherfeatures.df) ) {
    print('[ERROR] something went wrong between the main and subfeatures')
    print(paste(nrow(data.df),'!=',nrow(otherfeatures.df)))
  }
  
  # look if the main features and the exon features have the same rows 
  if ( nrow(data.df) != nrow(exon.data.df) ) {
    print('[ERROR] something went wrong between the main and exon features')
    print(paste(nrow(data.df),'!=',nrow(exon.data.df)))
  }
  
  # add features from otherfeatures to the feature vector
  if ( otherfeatures.b[3] == TRUE ) {
    print('[NOTE] you added the exon features')
    data.df <- cbind(data.df, exon.data.df)
  }
  
  additionalCols = c()
  if ( otherfeatures.b[1] == TRUE ) {
    print('[NOTE] you added the methylation ratio of the exons and introns as a feature')
    data.df <- cbind(data.df, as.numeric(otherfeatures.df[,3]) )
    colnames(data.df)[ncol(data.df)] <- 'methEx'
    data.df <- cbind(data.df, as.numeric(otherfeatures.df[,4]) )
    colnames(data.df)[ncol(data.df)] <- 'methUp'
    data.df <- cbind(data.df, as.numeric(otherfeatures.df[,5]) )
    colnames(data.df)[ncol(data.df)] <- 'methDo'
    additionalCols <- append(additionalCols, c('methEx', 'methUp', 'methDo'))
  }
  
  if ( otherfeatures.b[2] == TRUE ) {
    print('[NOTE] you added the gc ratio of the exons and introns as a feature')
    data.df <- cbind(data.df, as.numeric(otherfeatures.df[,6]))
    colnames(data.df)[ncol(data.df)] <- 'UpGC'
    data.df <- cbind(data.df, as.numeric(otherfeatures.df[,7]))
    colnames(data.df)[ncol(data.df)] <- 'DoGC'
    data.df <- cbind(data.df, as.numeric(otherfeatures.df[,8]))
    colnames(data.df)[ncol(data.df)] <- 'ExGC'
    additionalCols <- append(additionalCols, c('UpGC', 'DoGC', 'ExGC'))
  }
  
  if ( otherfeatures.b[4] == TRUE ) {
    print('[NOTE] you added the cytosine content of the exon and up-/downstream intron as a feature')
    data.df <- cbind(data.df, as.numeric(otherfeatures.df[,9]))
    colnames(data.df)[ncol(data.df)] <- 'UpC'
    data.df <- cbind(data.df, as.numeric(otherfeatures.df[,10]))
    colnames(data.df)[ncol(data.df)] <- 'DoC'
    data.df <- cbind(data.df, as.numeric(otherfeatures.df[,11]))
    colnames(data.df)[ncol(data.df)] <- 'ExC'
    additionalCols <- append(additionalCols, c('UpC', 'DoC', 'ExC'))
  }
  
  if ( otherfeatures.b[5] == TRUE ) {
    print('[NOTE] you added the Run Test p-value of the up-/downstream intron as a feature')
    
    # changes NAs to 0.0 (NAs appears for pattern with only -1, I treat them as a significant
    # not randomly distributed pattern)
    otherfeatures.df[which(is.na(otherfeatures.df[,13])), 13] <- 0.0
    otherfeatures.df[which(is.na(otherfeatures.df[,14])), 14] <- 0.0
    
    data.df <- cbind(data.df, as.numeric(otherfeatures.df[,13]))
    colnames(data.df)[ncol(data.df)] <- 'UpRun'
    data.df <- cbind(data.df, as.numeric(otherfeatures.df[,14]))
    colnames(data.df)[ncol(data.df)] <- 'DoRun'
    additionalCols <- append(additionalCols, c('UpRun', 'DoRun'))
  }
  
  # TODO change order
  if ( otherfeatures.b[10] == TRUE ) {
    print('[NOTE] you added the mCpG/CpG profile of the exons as additional features')
    data.df <- cbind(data.df, ratio_exon.features.df)
  }
  
  # take out splice sites if chosen 
  if ( otherfeatures.b[6] == FALSE ) {
    print('[NOTE] you removed 10 bp from up- and downstream intron around splice site and 1 bin up and down of the exon')
    pickcols.l <- c('target', paste0('up',1:(intronLength-10)), paste0('do',11:intronLength))
    if ( otherfeatures.b[3] == TRUE ) {
      pickcols.l <- append(pickcols.l, paste0('ex',2:(exBins-1)))
    }
    if ( otherfeatures.b[10] == TRUE ) {
      pickcols.l <- append(pickcols.l, paste0('exR',2:(exBins-1)))
    }
    pickcols.l <- append(pickcols.l, additionalCols)
    data.df <- data.df[,pickcols.l]
  }
  
  # take out methylation profiles 
  if ( otherfeatures.b[7] == TRUE ) {
    print('[NOTE] you removed the whole methylation profile for exons and introns')
    pickcols.l <- c('target')
    pickcols.l <- append(pickcols.l, additionalCols)
    data.df <- data.df[,pickcols.l]
  }
  
  if ( otherfeatures.b[8] == TRUE ) {
    print('[NOTE] you added the cytosine content of the exon and up-/downstream intron as a feature')
    data.df <- cbind(data.df, as.numeric(otherfeatures.df[,15]))
    colnames(data.df)[ncol(data.df)] <- 'UpmCpG/Cpg'
    data.df <- cbind(data.df, as.numeric(otherfeatures.df[,16]))
    colnames(data.df)[ncol(data.df)] <- 'ExmCp/Cpg'
    data.df <- cbind(data.df, as.numeric(otherfeatures.df[,17]))
    colnames(data.df)[ncol(data.df)] <- 'DomCpg/Cpg'
  }
  
  if ( otherfeatures.b[9] == TRUE ) {
    print('[NOTE] you added the methyaltion density profile of the introns as additional features')
    data_peaks.df <- peaks.df[,-1]               # strip off the target col
    pre <- ncol(data.df)+1
    data.df <- cbind(data.df, data_peaks.df)
    colnames(data.df)[pre:ncol(data.df)] <- c(paste0('d', 1:ncol(data_peaks.df)))
  }
  
  if ( otherfeatures.b[11] == TRUE ) {
    print('[NOTE] you removed the downstream methylation profile')
    pickcols.l <- c('target', paste0('up',1:intronLength))
    data.df <- data.df[,pickcols.l]
  }
  
  # get start and end position 
  chr <- c(features.exclusion.df[,21], features.inclusion.df[,21])
  startExon <- c(features.exclusion.df[,22], features.inclusion.df[,22])
  endExon <- c(features.exclusion.df[,23], features.inclusion.df[,23])
  strand <- c(features.exclusion.df[,24], features.inclusion.df[,24])
  
  # bind start and end position of exon 
  data.df <- cbind(data.df, 'ChrExon' = chr, 'StartExon' = startExon, 
                   'EndExon' = endExon, 'StrandExon' = strand)
  
  if ( c == 1 ){
    alldata.df <- data.df
    rownames(alldata.df) <- paste0(rownames(alldata.df), '_', celltype[c])
  } else {
    rownames(data.df) <- paste0(rownames(data.df), '_', celltype[c])
    alldata.df <- rbind(alldata.df, data.df)
  }
  
  print(celltype[c])
}

write.table(alldata.df, paste0('data/network_data_',
                               paste0(celltype, collapse = '_'),'.tsv'), sep =  "\t", append = FALSE, col.names=NA)



