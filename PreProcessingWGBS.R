
library(bsseq)
library(prodlim)     # for row matching

# celltype = 'H1hesc'
args = commandArgs(trailingOnly = TRUE)
celltype = args[1]

# set empirical thresholds
t1 = args[2]
t2 = args[3]

print(celltype)
print(t1)
print(t2)

####################
## Get WGBS data  ##
####################

# load methylation data 
# WGBS

if ( celltype == 'IMR90' ) {
  data <- read.delim("data/IMR90.Bisulfite-Seq.combined.bed", header=FALSE)
} else if ( celltype == 'Gm12878' ) {
  data <- read.delim("data/wgEncodeHaibMethylWgbsGm12878CpGSites.bed", header=FALSE)
} else if ( celltype == 'H1hesc' ) {
  data <- read.delim("data/wgEncodeHaibMethylWgbsH1hescCpGSites.bed", header=FALSE)
} else {
  stop('[ERROR] wrong celltype, use IMR90, Gm12878 or H1hesc')
}

# load annotation data for comparison (this is necessary due to werid 
# CpGs in the celltypes of the encode data see example below)

# V1  V2  V3   V4 V5 V6  V7  V8        V9 V10 V11
# 18037612 chr17  91  92 WGBS  2  -  91  92   0,255,0   2   0    <- CpG
# 18037613 chr17 141 142 WGBS  1  + 141 142   0,255,0   1   0    <- not a CpG
# 18037614 chr17 142 143 WGBS  6  - 142 143 255,105,0   6  83    <- CpG
# 18037615 chr17 172 173 WGBS  3  + 172 173 155,255,0   3  33    <- not a CpG
# 18037616 chr17 173 174 WGBS  9  - 173 174 255,155,0   9  77    <- CpG
# 18037617 chr17 184 185 WGBS  7  + 184 185  55,255,0   7  14    <- not a CpG

compare_data <- read.delim("data/CpG_whole_hg19.bed", header=FALSE)

# first three line are useless
head(compare_data)
compare_data = compare_data[-c(1:3),]
head(compare_data)

######################################################
## Compare between Annotation data and Encode data  ##
######################################################

# remove all rows which do not match with the annotation data
match.with.annotaton.l = row.match(data[,1:3],compare_data[,1:3])

# CpGs which immediately match
realCpGs.l = which(is.na(match.with.annotaton.l) == FALSE)

# CpG which do not match 
unrealCpG.l = which(is.na(match.with.annotaton.l))

# due to bed file it can happen that the position has to be shifted be 1 
unrealCpG_data = data[unrealCpG.l,]
unrealCpG_data[,2] = unrealCpG_data[,2] + 1
unrealCpG_data[,3] = unrealCpG_data[,3] + 1

# rowmatch between the modified CpG which do not match in the first round
match.unreal.with.annotaton.l = row.match(unrealCpG_data[,1:3],compare_data[,1:3])

# get now CpGS which match the annotation 
round2.realCpGs.l = which(is.na(match.unreal.with.annotaton.l) == FALSE)

# look again if some CpGs do not match
round2.unreal.CpG = which(is.na(match.unreal.with.annotaton.l))

if ( celltype == 'IMR90' ){
  unrealCpG_r3_data = data[round2.unreal.CpG,]
  unrealCpG_r3_data[,2] = unrealCpG_r3_data[,2] + 1
  
  # rowmatch between the modified CpG which do not match in the second round
  match.unreal.r3.with.annotaton.l = row.match(unrealCpG_r3_data[,1:3],compare_data[,1:3])
  
  # get now CpGS which match the annotation 
  round3.realCpGs.l = which(is.na(match.unreal.r3.with.annotaton.l) == FALSE)
  
  # look again if some CpGs do not match
  round3.unreal.CpG = which(is.na(match.unreal.r3.with.annotaton.l))
}

# get from first round the real CpGS
if ( length(realCpGs.l) != 0 ){
  new_data = data[realCpGs.l,]

  # bind from second round the unreal modified CpGs
  new_data = rbind(new_data, unrealCpG_data[round2.realCpGs.l,])
} else {
  
  if ( length(round2.realCpGs.l) != 0 ){
    new_data = unrealCpG_data[round2.realCpGs.l,]
  } else {
    new_data = unrealCpG_r3_data[round3.realCpGs.l,]
  }
}

# look for unique CpGs (just in case)
duplicated_entries = which(duplicated(new_data[,1:3]))

# remove duplicated entries
if (  length(duplicated_entries) != 0  ){
  print(paste('[NOTE] remove',length(duplicated_entries),'duplicated entries'))
  new_data = new_data[-duplicated_entries,]
} 

write.table(new_data, file = paste0('tmp/CpG_',celltype,'.tsv'), sep='\t', append=FALSE, row.names=FALSE)

################################
## Change matrix of CpG sites ##
################################

ch_data = new_data
plotcellname <- ''

if ( celltype == 'IMR90' ){
  colnames(ch_data) = c('chrom','start','end','probability of methylation')
  plotcellname <- 'IMR-90'
}

if ( celltype == 'Gm12878' ){
  colnames(ch_data) = c('chrom','start','end','name',
                                'cov','strand','start','end',
                                'RBG','cov','probability of methylation')
  ch_data[,11] = ch_data[,11] / 100
  plotcellname <- 'GM12878'
}

if ( celltype == 'H1hesc' ){
  colnames(ch_data) = c('chrom','start','end','name',
                                'cov','strand','start','end',
                                'RBG','cov','probability of methylation')
  ch_data[,11] = ch_data[,11] / 100
  plotcellname <- 'H1-hESC'
}

##############
## Filter A ##
##############

CpGs = ch_data

# throw out all CpGs related to chromosome X and Y
chromosomes = CpGs[,1]
removeXY = c(which(chromosomes == 'chrX'), which(chromosomes == 'chrY'))

# from CpGs
CpGs = CpGs[-removeXY,]

###############
## Smoothing ##
###############

if ( celltype != 'IMR90' ) {

  # look if smoothed data already exist 
  filename = paste0('/tmp/BismoothedValues_',celltype,'.tsv')
  
  if ( file.exists(filename) == FALSE ) {
  
    if ( celltype == 'IMR90' ) {
      #M = matrix(CpGs[,4] * CpGs[,11])
      #cov = matrix(CpGs[,4])
    } else {
      M = matrix(CpGs[,5] * CpGs[,11])
      cov = matrix(CpGs[,5])
    }
    
    BStmp = BSseq( chr=CpGs[,1] , pos=CpGs[,2], M=M, Cov=cov )
    
    # 6 cores to speed up the process
    Bsm = BSmooth(BStmp, verbose = TRUE, mc.cores=10)
    
    # Better use coefficient from bsmooth
    # With smoothing (using the 'coef' slot to locally weight)
    SmoothedValues = getMeth(Bsm, type='smooth', what='perBase')
    
    print('[NOTE] smoothing finished')
    
    write.table(as.data.frame(SmoothedValues), 
                file = filename, sep='\t', append=FALSE)
    
    print('[NOTE] wrote smoothed matrix')
  } else {
    SmoothedValues = read.delim(filename, header=TRUE)
  }
}
####################################################
## Plot distribution of Methylation Probabilities ##
####################################################

if ( celltype == 'Gm12878' || celltype == 'H1hesc'  ) {
  CpGs = CpGs[,c(1,2,3,11)]
  
  # add smoothed values
  CpGs = cbind(CpGs, SmoothedValues)
  colnames(CpGs)[5] = 'smoothed Values'
}

# throw out all CpGs related to chromosome M
chromosomes = CpGs[,1]
removeM = c(which(chromosomes == 'chrM'))

# from CpGs
CpGs = CpGs[-removeM,]

# list of all chromosomes and "all" = whole genome (except chromosome X)
listplots.l = c('Whole Genome (Except X, Y and M)', unique(as.character(CpGs[,1])))

for ( i in 1:length(listplots.l) ) {
  
  print(i)
  
  # change which CpGs you want to analyse and change file name of png
  if ( i != 1  ){
    name = listplots.l[i]
    lookatCpGs = CpGs[which(CpGs[,1] %in% listplots.l[i]),]
  } else {
    name = 'all'
    lookatCpGs = CpGs
  }

  # get rows for mCpGs
  rows_mCpGs = intersect( which(lookatCpGs[,4] >= t1), which(lookatCpGs[,4] <= t2))

  h = 700
  if (celltype == 'IMR90') {
	h = 350
  } 
  
  # create plots
  png(filename = paste0('plots/',celltype,'/histogramWGBS_',
                        name,'.png'), width = 800, height = h)
  par(family = 'serif')
  if ( celltype != 'IMR90' ){
  
    # get rows for mCpGs based on smoothed values 
    sm_rows_mCpGs = intersect( which(lookatCpGs[,5] >= t1), which(lookatCpGs[,5] <= t2))
    
    # create plots
    par(mfrow=c(2, 1))
    
    plot(density(lookatCpGs[,5], adjust=0.2, na.rm=TRUE),
         main=paste('mCpG Smoothed Probability of \n', nrow(lookatCpGs),
                    'CpGs',listplots.l[i],'for',plotcellname),
         xlab='Pobability of Methylation',
         ylab='Density',
         ylim=c(0.0,15.0), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
    abline(v=t1, col='red')
    abline(v=t2, col='red')
    text(x=0.4,y=10.0,paste('#CpGs',length(sm_rows_mCpGs)),col='red', cex = 1.5)
  }
  
  plot(density(lookatCpGs[,4], adjust=0.2, na.rm=TRUE),
       main=paste('mCpG Probability of \n', nrow(lookatCpGs),
                  'CpGs',listplots.l[i],'for',plotcellname),
       xlab='Pobability of Methylation',
       ylab='Density',
       ylim=c(0.0,15.0), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
  abline(v=t1, col='red')
  abline(v=t2, col='red')
  text(x=0.4,y=10.0,paste('#CpGs',length(rows_mCpGs)),col='red', cex = 1.5)
  
  dev.off()
}

###########################
## Write Methylated CpGs ##
###########################

# write out GenomeRanges as a table for mCpgs
if (  celltype == 'IMR90' ) {
  # use the unsmoothed values to get methylated CpGs
  mCpGs = CpGs[intersect(which(CpGs[,4] >= t1), which(CpGs[,4] <= t2)),]
  
  print(nrow(mCpGs)/nrow(CpGs))

  write.table(mCpGs, file = 'tmp/mCpG_IMR90.tsv', sep='\t', append=FALSE, row.names=FALSE)
}


if (  celltype == 'Gm12878' ) {
  # use the unsmoothed values to get methylated CpGs
  mCpGs = CpGs[intersect(which(CpGs[,4] >= t1), which(CpGs[,4] <= t2)),]

  print(nrow(mCpGs))
  print(nrow(mCpGs)/nrow(CpGs))
  
  # take out pobability of methylation
  mCpGs = mCpGs[,-5]
  mCpGs = mCpGs[,-4]

  write.table(mCpGs, file = 'tmp/mCpG_Gm12878.tsv', sep='\t', append=FALSE, row.names=FALSE)
}

if (  celltype == 'H1hesc' ) {
  # use the unsmoothed values to get methylated
  mCpGs = CpGs[intersect(which(CpGs[,4] >= t1), which(CpGs[,4] <= t2)),]

  print(nrow(mCpGs))
  print(nrow(mCpGs)/nrow(CpGs))
  
  # take out pobability of methylation
  mCpGs = mCpGs[,-5]
  mCpGs = mCpGs[,-4]

  write.table(mCpGs, file = 'tmp/mCpG_H1hesc.tsv', sep='\t', append=FALSE, row.names=FALSE)
}

print('[FINISH]')
