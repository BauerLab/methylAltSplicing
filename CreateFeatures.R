
#####################
## Load libraries  ##
#####################

library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Rmisc)
library(bsseq)
library(tseries)  # for runs test of randomness
library(plyr)
library(data.table)   # for quick import of data

####################
## Get meth data  ##
####################

celltype = 'H1hesc'

if ( celltype == 'IMR90' ) { 
  # load methylation data 
  # WGBS for methylated CpGs
  mCpGs.df = fread('tmp/mCpG_IMR90.tsv', sep='\t', header=TRUE)
  
  # remove chromosome M CpGs 
  keepCpGs = which(as.character(mCpGs.df[[1]]) != 'chrM')
  mCpGs.df = mCpGs.df[keepCpGs,]
  
  # create GRanges object
  mCpGs.gr = GRanges(seqnames = mCpGs.df[[1]], 
                     ranges = IRanges(start = as.numeric(mCpGs.df[[2]]), 
                                      end = as.numeric(mCpGs.df[[2]])), 
                     strand = rep('*',nrow(mCpGs.df)) )
  
  # remove big data from memory
  rm(mCpGs.df)
  
  # load methylation data 
  # WGBS (keep in mind that this is bed format see BismoothData.v*.R)
  CpGs.df <- fread('tmp/CpG_IMR90.tsv', sep='\t', header=TRUE)
  
  # remove chromosome M CpGs 
  keepCpGs = which(as.character(CpGs.df[[1]]) != 'chrM')
  CpGs.df = CpGs.df[keepCpGs,]
  
  # create GRanges object
  CpGs.gr = GRanges(seqnames = CpGs.df[[1]], 
                    ranges = IRanges(start = as.numeric(CpGs.df[[2]]), # changed start site
                                     end = as.numeric(CpGs.df[[2]])),  # same end position
                    strand = rep('*',nrow(CpGs.df)) )
  
  # remove big data from memory
  rm(CpGs.df)
  
  # read data from PSI.R
  psi.df = read.delim("data/PSI.10.IMR90.tsv", header=T)
  row.names(psi.df) = NULL
}

if ( celltype == 'Gm12878' ) {
  # load methylation data 
  # WGBS for methylated CpGs
  mCpGs.df = fread('tmp/mCpG_Gm12878.tsv', sep='\t', header=TRUE)
  
  # remove chromosome M CpGs 
  keepCpGs = which(as.character(mCpGs.df[[1]]) != 'chrM')
  mCpGs.df = mCpGs.df[keepCpGs,]
  
  # create GRanges object
  mCpGs.gr = GRanges(seqnames = mCpGs.df[[1]], 
                     ranges = IRanges(start = as.numeric(mCpGs.df[[2]]), 
                                      end = as.numeric(mCpGs.df[[2]])), 
                     strand = rep('*',nrow(mCpGs.df)) )
  
  # remove big data from memory
  rm(mCpGs.df)
  
  # load methylation data 
  # WGBS (keep in mind that this is bed format see BismoothData.v*.R)
  CpGs.df <- fread('tmp/CpG_Gm12878.tsv', sep='\t', header=TRUE)
  
  # remove chromosome M CpGs 
  keepCpGs = which(as.character(CpGs.df[[1]]) != 'chrM')
  CpGs.df = CpGs.df[keepCpGs,]
  
  # create GRanges Object
  CpGs.gr = GRanges(seqnames = CpGs.df[[1]], 
                    ranges = IRanges(start = as.numeric(CpGs.df[[2]]), # changed start site
                                     end = as.numeric(CpGs.df[[2]])),  # same end position
                    strand = rep('*',nrow(CpGs.df)) )
  
  # remove big data from memory
  rm(CpGs.df)
  
  # read data from PSI.R
  psi.df = read.delim("data/PSI.10.Gm12878.tsv", header=T)
  row.names(psi.df) = NULL
}

if ( celltype == 'H1hesc' ) {
  # load methylation data 
  # WGBS for methylated CpGs
  mCpGs.df = fread('tmp/mCpG_H1hesc.tsv', sep='\t', header=TRUE)
  
  # remove chromosome M CpGs 
  keepCpGs = which(as.character(mCpGs.df[[1]]) != 'chrM')
  mCpGs.df = mCpGs.df[keepCpGs,]
  
  # create GRanges object
  mCpGs.gr = GRanges(seqnames = mCpGs.df[[1]], 
                     ranges = IRanges(start = as.numeric(mCpGs.df[[2]]), 
                                      end = as.numeric(mCpGs.df[[2]])), 
                     strand = rep('*',nrow(mCpGs.df)) )
  
  # remove big data from memory
  rm(mCpGs.df)
  
  # load methylation data 
  # WGBS (keep in mind that this is bed format see BismoothData.v*.R)
  CpGs.df <- fread('tmp/CpG_H1hesc.tsv', sep='\t', header=TRUE)
  
  # remove chromosome M CpGs 
  keepCpGs = which(as.character(CpGs.df[[1]]) != 'chrM')
  CpGs.df = CpGs.df[keepCpGs,]
  
  # create GRanges object 
  CpGs.gr = GRanges(seqnames = CpGs.df[[1]], 
                    ranges = IRanges(start = as.numeric(CpGs.df[[2]]), # changed start site
                                     end = as.numeric(CpGs.df[[2]])),  # same end position
                    strand = rep('*',nrow(CpGs.df)) )
  
  # remove big data from memory
  rm(CpGs.df)
  
  # read data from PSI.R
  psi.df = read.delim("data/PSI.10.H1hesc.tsv", header=T)
  row.names(psi.df) = NULL
}

#####################
## Get other data  ##
#####################

# get exons for hg19 
all.exons.df = fread('tmp/sorted_exons.bed', sep='\t', header=FALSE)

#####################
## Get other data  ##
#####################

# get exons for hg19 
all.exons.df = fread('tmp/sorted_exons.bed', sep='\t', header=FALSE)

###################
## Get PSI data  ##
###################

# read data from PSI.R
PSI.exclusion.df = read.delim(file=paste0('data/unfiltered/PSI.unfiltered.chosenExons.FALSE.',celltype,'.tsv'), header=T)
row.names(PSI.exclusion.df) = NULL

PSI.inclusion.df = read.delim(file=paste0('data/unfiltered/PSI.unfiltered.chosenExons.TRUE.',celltype,'.tsv'), header=T)
row.names(PSI.inclusion.df) = NULL


########################
## Start of Pipeline  ##
########################

pipe = function(cs, shortening, shift, inlength, exonBins) {

  # read controls
  PSI.exclusion.control.df = read.delim(file=paste0('data/filtered/PSI.filtered.chosenExons.FALSE.0.',inlength,'.',celltype,'.tsv'), header=T)
  PSI.inclusion.control.df = read.delim(file=paste0('data/filtered/PSI.filtered.chosenExons.TRUE.0.',inlength,'.',celltype,'.tsv'), header=T)
  
  ###############
  ## Pipeline  ##
  ###############
  
  chosenExons = PSI.exclusion.df
  control = PSI.exclusion.control.df
  if ( cs ) {
    chosenExons = PSI.inclusion.df
    control = PSI.inclusion.control.df
  }
  
  #################
  ## Get Introns ##
  #################
  
  if ( length(which(chosenExons[,4] == '*')) != 0 ){
    stop('[Error] found undefined strand in one exon entry')
  }
  
  print(paste0('[NOTE] number of exons for testing = ',nrow(chosenExons)))
  
  # replicate the exons 
  Upintrons.df = Dointrons.df = chosenExons
  rownames(Upintrons.df) = rownames(Dointrons.df) = NULL
  
  intronLength = inlength
  s1 = 1 + shift
  s2 = 1 + shift
  if(shortening) {
    s1 = 20
    s2 = 6
  }
  
  # first end pos, second start pos (+ upstream introns)
  Upintrons.df[which(Upintrons.df[,4] == "+"),3] = as.numeric(Upintrons.df[which(Upintrons.df[,4] == "+"),2]) - s1
  Upintrons.df[which(Upintrons.df[,4] == "+"),2] = as.numeric(Upintrons.df[which(Upintrons.df[,4] == "+"),2]) - intronLength - (s1 - 1)
  
  # first start pos, second end pos (- downstream introns)
  Dointrons.df[which(Dointrons.df[,4] == "-"),3] = as.numeric(Dointrons.df[which(Dointrons.df[,4] == "-"),2]) - s2
  Dointrons.df[which(Dointrons.df[,4] == "-"),2] = as.numeric(Dointrons.df[which(Dointrons.df[,4] == "-"),2]) - intronLength - (s2 - 1)
  
  # first start pos, second end pos (+ downstream introns)
  Dointrons.df[which(Dointrons.df[,4] == "+"),2] = as.numeric(Dointrons.df[which(Dointrons.df[,4] == "+"),3]) + s2
  Dointrons.df[which(Dointrons.df[,4] == "+"),3] = as.numeric(Dointrons.df[which(Dointrons.df[,4] == "+"),3]) + intronLength + (s2 - 1)
  
  # first end pos, second start pos (- upstream introns)
  Upintrons.df[which(Upintrons.df[,4] == "-"),2] = as.numeric(Upintrons.df[which(Upintrons.df[,4] == "-"),3]) + s1
  Upintrons.df[which(Upintrons.df[,4] == "-"),3] = as.numeric(Upintrons.df[which(Upintrons.df[,4] == "-"),3]) + intronLength + (s1 - 1)
  
  # create GRanges objects
  all.exons.gr = GRanges(seqnames=all.exons.df[[1]],
                         ranges=IRanges(start=as.numeric(all.exons.df[[2]]),
                                        end=as.numeric(all.exons.df[[3]])),
                         strand=as.character(all.exons.df[[6]]))
  
  Upintrons.gr = GRanges(seqnames=Upintrons.df[,1],
                         ranges=IRanges(start=as.numeric(Upintrons.df[,2]),
                                        end=as.numeric(Upintrons.df[,3])),
                         strand=Upintrons.df[,4])
  
  Dointrons.gr = GRanges(seqnames=Dointrons.df[,1],
                         ranges=IRanges(start=as.numeric(Dointrons.df[,2]),
                                        end=as.numeric(Dointrons.df[,3])),
                         strand=Dointrons.df[,4])
  
  
  # remove introns which overlaps with exons
  overlaps1 = findOverlaps(all.exons.gr, Upintrons.gr)
  removeIntrons1 = unique(subjectHits(overlaps1))
  
  overlaps2 = findOverlaps(all.exons.gr, Dointrons.gr)
  
  removeIntrons2 = unique(c(removeIntrons1,unique(subjectHits(overlaps2))))
  filtered_Upintrons.df = Upintrons.df[-removeIntrons2,]
  filtered_Dointrons.df = Dointrons.df[-removeIntrons2,]
  chosenExons = chosenExons[-removeIntrons2,]
  
  print(paste0('[NOTE] number of upstream introns for testing = ',nrow(filtered_Upintrons.df)))
  print(paste0('[NOTE] number of exons for testing = ',nrow(chosenExons)))
  print(paste0('[NOTE] number of downstream introns for testing = ',nrow(filtered_Dointrons.df)))
  
  if ( nrow(chosenExons) != nrow(control) ){
    print('[INFO] attention you used a different shift or shortening')
    print(paste('control with',nrow(control),'exons'))
  }
  
  #########################
  ## Methylation content ##
  #########################
  
  # create GRanges object out of the exons
  grExons <- GRanges(seqnames = chosenExons[,1], 
                     ranges = IRanges(start = as.numeric(chosenExons[,2]), 
                                      end = as.numeric(chosenExons[,3])), 
                     strand = chosenExons[,4])
  
  # count methylation loci per exon
  overlapcountsExons <- countOverlaps(grExons,mCpGs.gr)
  
  # create datatable
  Ematrix <- matrix(ncol=6, nrow=length(grExons))
  colnames(Ematrix) <- c("ID", "#Methylathions", "seqlength", 
                         'lengthNorm', 'CpG', 'mCpG/CpG')
  
  Ematrix[,1] <- 1:nrow(Ematrix)
  Ematrix[,2] <- overlapcountsExons
  Ematrix[,3] <- width(grExons)
  Ematrix[,4] <- (as.numeric(Ematrix[,2]) / as.numeric(Ematrix[,3]))
  
  # calculate mCpG/CpG-ratio for exons
  # count overlap with all CpGs
  overlapcountsExons <- countOverlaps(grExons,CpGs.gr)
  Ematrix[,5] = overlapcountsExons
  Ematrix[,6] = (as.numeric(Ematrix[,2]) / as.numeric(Ematrix[,5]))
  
  # check for NaNs
  Ematrix[which(is.na(Ematrix[,6])),6] = -1
  
  # print how many exons are not 100% or 0% methylated 
  print(paste('[DATA] number of exons with mCpG/CpG ratio != 0 or 1 :', length(which(!Ematrix[,6] %in% c(-1,1,0))) ) )
  
  
  # check for Inf
  if ( length(which(is.infinite(Ematrix[,6]))) != 0 ) {
    print('[ERROR] something went wrong with the calculation of the mCpG/CpG Ratio exons, infinite valeus')
  }
  
  # check if mCpG/CpG makes sense
  if( length(which(Ematrix[,6] > 1.0)) != 0 ){
    print('[ERROR] something went wrong with the calculation of the mCpG/CpG Ratio exons, values bigger than 1')
  }
  
  ################
  ### UPSTREAM ### 
  ################
  
  # create GRanges object out of the upstream introns
  grUpIntrons <- GRanges(seqnames = filtered_Upintrons.df[,1], 
                         ranges = IRanges(start = as.numeric(filtered_Upintrons.df[,2]), 
                                          end = as.numeric(filtered_Upintrons.df[,3])), 
                         strand = filtered_Upintrons.df[,4])
  
  # count methylation loci per exon
  overlapcountsUpIntrons <- countOverlaps(grUpIntrons,mCpGs.gr)
  
  # create datatable
  UpImatrix <- matrix(ncol=6, nrow=length(grUpIntrons))
  colnames(UpImatrix) <- c("ID", "#Methylathions", "seqlength", 
                           'lengthNorm', 'CpG', 'mCpG/CpG')
  
  UpImatrix[,1] <- 1:nrow(UpImatrix)
  UpImatrix[,2] <- overlapcountsUpIntrons
  UpImatrix[,3] = width(grUpIntrons)
  UpImatrix[,4] = (as.numeric(UpImatrix[,2]) / as.numeric(UpImatrix[,3]))
  
  # calculate mCpG/CpG-ratio for upstream introns
  # count overlap with all CpGs
  overlapcountsUpIntrons <- countOverlaps(grUpIntrons,CpGs.gr)
  UpImatrix[,5] = overlapcountsUpIntrons
  UpImatrix[,6] = (as.numeric(UpImatrix[,2]) / as.numeric(UpImatrix[,5]))
  
  # check for NaNs
  UpImatrix[which(is.na(UpImatrix[,6])),6] = -1
  
  # print how many upstream introns are not 100% or 0% methylated 
  print(paste('[DATA] number of upstream introns with mCpG/CpG ratio != 0 or 1 :', length(which(!UpImatrix[,6] %in% c(-1,1,0))) ) )
  
  # check for Inf
  if ( length(which(is.infinite(UpImatrix[,6]))) != 0 ) {
    print('[ERROR] something went wrong with the calculation of the mCpG/CpG Ratio exons, infinite valeus')
  }
  
  # check if mCpG/CpG makes sense
  if( length(which(UpImatrix[,6] > 1.0)) != 0 ){
    print('[ERROR] something went wrong with the calculation of the mCpG/CpG Ratio up introns, values bigger than 1')
  }
  
  ################
  ### DOSTREAM ### 
  ################
  
  # create GRanges object out of the downstream introns
  grDoIntrons <- GRanges(seqnames = filtered_Dointrons.df[,1], 
                         ranges = IRanges(start = as.numeric(filtered_Dointrons.df[,2]), 
                                          end = as.numeric(filtered_Dointrons.df[,3])), 
                         strand = filtered_Dointrons.df[,4])
  
  # count methylation loci per exon
  overlapcountsDoIntrons <- countOverlaps(grDoIntrons,mCpGs.gr)
  
  # create datatable
  DoImatrix <- matrix(ncol=6, nrow=length(grDoIntrons))
  colnames(DoImatrix) <- c("ID", "#Methylathions", "seqlength", 
                           'lengthNorm', 'CpG', 'mCpG/CpG')
  
  DoImatrix[,1] <- 1:nrow(DoImatrix)
  DoImatrix[,2] <- overlapcountsDoIntrons
  DoImatrix[,3] = width(grDoIntrons)
  DoImatrix[,4] = (as.numeric(DoImatrix[,2]) / as.numeric(DoImatrix[,3]))
  
  
  # calculate mCpG/CpG-ratio for downstream introns
  # count overlap with all CpGs
  overlapcountsDoIntrons <- countOverlaps(grDoIntrons,CpGs.gr)
  DoImatrix[,5] = overlapcountsDoIntrons
  DoImatrix[,6] = (as.numeric(DoImatrix[,2]) / as.numeric(DoImatrix[,5]))
  
  # check for NaNs
  DoImatrix[which(is.na(DoImatrix[,6])),6] = -1
  
  # print how many downstream introns are not 100% or 0% methylated 
  print(paste('[DATA] number of downstream introns with mCpG/CpG ratio != 0 or 1 :', length(which(!DoImatrix[,6] %in% c(-1,1,0))) ) )
  
  # check for Inf
  if ( length(which(is.infinite(DoImatrix[,6]))) != 0 ) {
    print('[ERROR] something went wrong with the calculation of the mCpG/CpG Ratio exons, infinite valeus')
  }
  
  # check if mCpG/CpG makes sense
  if( length(which(DoImatrix[,6] > 1.0)) != 0 ){
    print('[ERROR] something went wrong with the calculation of the mCpG/CpG Ratio do introns, values bigger than 1')
  }
  
  #################################
  ## randomness of methyaltions  ##
  #################################
  
  runsTest = function(x) {
    
    # check if array is not only a zero vector
    if ( length(which(x != 0)) != 0 ) {
      return(runs.test(factor(x))[[3]])
    } else {
      return(0.0)
    }
    
  }
  
  ###########################################
  ## Binning Methylation content for Exons ##
  ###########################################
  
  if ( exonBins != 0 ) {
    
    # merge overlaps 
    Exmerge = mergeByOverlaps(grExons, mCpGs.gr)
    
    # e.g.  
    # grExons       lean.smooth.rowRanges
    # <GRanges>                   <GRanges>
    #   1 chrX:-:[99885755, 99885863] chrX:*:[99885829, 99885829]
    # 2 chrX:-:[99885755, 99885863] chrX:*:[99885845, 99885845]
    # 3 chrX:-:[99888401, 99888536] chrX:*:[99888515, 99888515]
    # 4 chrX:-:[99890554, 99890743] chrX:*:[99890580, 99890580]
    # 5 chrX:-:[99890554, 99890743] chrX:*:[99890623, 99890623]
    # 6 chrX:-:[99890554, 99890743] chrX:*:[99890644, 99890644]
    
    # for exons chunk the exons width in bins 
    num.bins = exonBins
    
    # get exons width
    width.exons = width(Exmerge[,1])
    
    # get bin width for individual exon depending on the exon width 
    width.bin = width.exons / num.bins
    
    # get relative position inside of the exon regardin exon width 
    rel_positions = abs(start(Exmerge[,1]) - start(Exmerge[,2]))
    
    # negative strand
    neg = which(strand(Exmerge[,1]) == '-')
    rel_positions[neg] = end(Exmerge[neg,1]) - start(Exmerge[neg,2])
    
    # get the number of the bin the methylation is in 
    rel_bin_position = ceiling(rel_positions / width.bin)
    
    # matrix for bin-exon-length normalization 
    norm.df = matrix(nrow=length(rel_positions), ncol=4)
    norm.df[,1] = as.character(Exmerge[,1])
    norm.df[,2] = rel_bin_position
    norm.df[,3] = width.bin
    norm.df[,4] = 1/width.bin
    colnames(norm.df) = c('ID', 'bin', 'width_bin','inverse_width_bin')
    
    # because there is a zero bin I have to change the bin position by one
    norm.df[,2] = as.numeric(norm.df[,2]) + 1
    
    # because some values overshoot a bit and land in the bin which is 
    # (length + 1) I have to put them to the last bin 
    norm.df[which(norm.df[,2] == as.character(num.bins + 1)),2] = as.character(num.bins)
    
    # unique exon matches found by GenomicRanges
    uniqueEntries.ex.mCpG = unique(Exmerge[,1])
    
    # matrix containing bit-wise-vector for exon regions 
    x.ex.mCpG.m = matrix(0,nrow=length(uniqueEntries.ex.mCpG), ncol=exonBins)
    
    # array for a later check up
    doubleEntryIds.l = c()
    
    # got over all bins  
    for(i in 1:exonBins){
      posFori = which(norm.df[,2] == i)
      
      # for an exon the bin could be matched more than one time with a methylation
      IDs.df = as.data.frame(norm.df[posFori,])
      count_posFori = count(IDs.df[,1])  

      # 
      #     1   chr1:16254585-16262761:+    1
      #     2 chr1:204437997-204439014:-    1
      #     3   chr1:39852858-39854330:+    1
      #     4  chr10:13043196-13043697:-    1
      #     5  chr11:62283376-62301546:-    2
      #     6  chr11:67933181-67934645:-    1
      #     which(norm.df[,1] == 'chr11:62283376-62301546:-')
      #     norm.df[38799 ,4]
      #     "0.0275163722414837"
      
      # fill array for a later check up
      doubleEntryIds.l = append(doubleEntryIds.l, 
                                as.character(count_posFori[which(count_posFori[,2] != 1),1]))
      
      # get values for the bins
      uq_IDs.df = IDs.df[which(duplicated(IDs.df[,1]) == F),]
  
      # match ids between the count dataframe and the reduced Id dataframe 
      matching = match(as.character(uq_IDs.df[,1]), as.character(count_posFori[,1]))
  
      # multiply the count to the inverse bin width 
      uq_IDs.df[,4] = as.numeric(as.character(uq_IDs.df[,4])) * as.numeric(count_posFori[matching, 2])
  
      # collect values 
      rows = which(uniqueEntries.ex.mCpG %in% Exmerge[posFori,1])
      x.ex.mCpG.m[rows,i] = uq_IDs.df[,4]
      if(i%%10 == 0) print(i)
    }
  }
  
  # list of matches between exons and methylation 
  # some exons are not listed 
  exonMatches.mCpGs.l = match(uniqueEntries.ex.mCpG, grExons)
  
  ######################################
  ## Binning mCpG/CpG Ratio for Exons ##
  ######################################
  
  if ( exonBins != 0 ) {
    
    # merge overlaps with CpGs
    Exmerge.CpG = mergeByOverlaps(grExons, CpGs.gr)

    # get exons width
    width.exons.CpG = width(Exmerge.CpG[,1])
    
    # get bin width for individual exon depending on the exon width 
    width.bin.CpG = width.exons.CpG / num.bins
    
    # get relative position inside of the exon regardin exon width 
    rel_positions.CpG = abs(start(Exmerge.CpG[,1]) - start(Exmerge.CpG[,2]))
    
    # negative strand
    neg.CpG = which(strand(Exmerge.CpG[,1]) == '-')
    rel_positions.CpG[neg.CpG] = end(Exmerge.CpG[neg.CpG,1]) - start(Exmerge.CpG[neg.CpG,2])
    
    # get the number of the bin the methylation is in 
    rel_positions.CpG = ceiling(rel_positions.CpG / width.bin.CpG)
    
    # matrix for bin-exon-length normalization 
    norm.CpG.df = matrix(nrow=length(rel_positions.CpG), ncol=4)
    norm.CpG.df[,1] = as.character(Exmerge.CpG[,1])
    norm.CpG.df[,2] = rel_positions.CpG
    norm.CpG.df[,3] = width.bin.CpG
    norm.CpG.df[,4] = 1/width.bin.CpG
    colnames(norm.CpG.df) = c('ID', 'bin', 'width_bin','inverse_width_bin')
    
    # because there is a zero bin I have to change the bin position by one
    norm.CpG.df[,2] = as.numeric(norm.CpG.df[,2]) + 1
    
    # because some values overshoot a bit and land in the bin which is 
    # (length + 1) I have to put them to the last bin 
    norm.CpG.df[which(norm.CpG.df[,2] == as.character(num.bins + 1)),2] = as.character(num.bins)
    
    # unique exon matches found by GenomicRanges
    uniqueEntries.ex.CpG = unique(Exmerge.CpG[,1])
    
    # matrix containing bit-wise-vector for exon regions 
    x.ex.CpG.m = matrix(0,nrow=length(uniqueEntries.ex.CpG), ncol=exonBins)
    
    # got over all bins 
    for(i in 1:exonBins){
      posFori.CpG = which(norm.CpG.df[,2] == i)
      
      # for an exon the bin could be matched more than one time with a methylation
      IDs.CpG.df = as.data.frame(norm.CpG.df[posFori.CpG,])
      count_posFori.CpG = count(IDs.CpG.df[,1])  
      
      # 
      #     1   chr1:16254585-16262761:+    1
      #     2 chr1:204437997-204439014:-    1
      #     3   chr1:39852858-39854330:+    1
      #     4  chr10:13043196-13043697:-    1
      #     5  chr11:62283376-62301546:-    2
      #     6  chr11:67933181-67934645:-    1
      #     which(norm.df[,1] == 'chr11:62283376-62301546:-')
      #     norm.df[38799 ,4]
      #     "0.0275163722414837"
      
      # fill array for a later check up
      doubleEntryIds.l = append(doubleEntryIds.l, 
                                as.character(count_posFori.CpG[which(count_posFori.CpG[,2] != 1),1]))
      
      # get values for the bins
      uq_IDs.CpG.df = IDs.CpG.df[which(duplicated(IDs.CpG.df[,1]) == F),]
      
      # match ids between the count dataframe and the reduced Id dataframe 
      matching.CpG = match(as.character(uq_IDs.CpG.df[,1]), as.character(count_posFori.CpG[,1]))
      
      # multiply the count to the inverse bin width 
      uq_IDs.CpG.df[,4] = as.numeric(as.character(uq_IDs.CpG.df[,4])) * 
                          as.numeric(count_posFori.CpG[matching.CpG, 2])
      
      # collect values 
      rows.CpG = which(uniqueEntries.ex.CpG %in% Exmerge.CpG[posFori.CpG,1])
      x.ex.CpG.m[rows.CpG,i] = uq_IDs.CpG.df[,4]
      if(i%%10 == 0) print(i)
    }
  }
  
  # list of matches between mCpGs and CpGs exons
  # some exons are not listed 
  exonMatches.CpG.l = match(uniqueEntries.ex.mCpG, uniqueEntries.ex.CpG)
  
  # calculate the mCpG/CpG ratio for each bin
  x.ex.ratio.m = x.ex.CpG.m
  
  # take elementwise inverse 
  inverse_x.ex.ratio.m = 1/x.ex.ratio.m

  # change all x.ex.mCpG.m zeros to -1 to capture site which cant be methylated
  # min1_x.ex.mCpG.m = x.ex.mCpG.m
  # min1_x.ex.mCpG.m[which(min1_x.ex.mCpG.m == 0 )] = -1
  
  # calulate new ratio profile
  x.ex.ratio.m = inverse_x.ex.ratio.m
  x.ex.ratio.m[exonMatches.CpG.l,] = inverse_x.ex.ratio.m[exonMatches.CpG.l,] * x.ex.mCpG.m

  # for rows where there is only CpG sites but no mCpGs
  no_mCpGS.l = which(is.na(match(uniqueEntries.ex.CpG, uniqueEntries.ex.mCpG)))
  x.ex.ratio.m[no_mCpGS.l,] = x.ex.ratio.m[no_mCpGS.l,] * rep(0.0, ncol(x.ex.ratio.m))
  
  # every entry which is NaN can not be methylated ---> -1
  x.ex.ratio.m[which(is.na(x.ex.ratio.m))] = -1
  
  if ( length(which(is.infinite(x.ex.ratio.m))) != 0  ) {
    print('[ERROR] something went wrong with the ratio matrix, there are infinite values')
  }

  # list of matches between exons and CpGs 
  # some exons are not listed 
  exonMatches.CpGs.l = match(uniqueEntries.ex.CpG, grExons)
  
  # match array with exon which have double entires for a later check up
  doubleEntrySeq = matrix(unlist(strsplit(doubleEntryIds.l, split='\\:|\\-')), ncol=4, byrow=TRUE)
  doubleEntrySeq[which(doubleEntrySeq[,4] == ''),4] = '-'
  doubleEntry.gr = GRanges(seqnames = doubleEntrySeq[,1], 
                           ranges = IRanges(start = as.numeric(doubleEntrySeq[,2]), 
                                            end = as.numeric(doubleEntrySeq[,3])), 
                           strand = doubleEntrySeq[,4])
  doubleEntryMatchIDs = match(doubleEntry.gr, grExons)
  
  ##################
  ### DOWNSTREAM ### 
  ##################
  
  ### CPG
  Domerge.CpG = mergeByOverlaps(grDoIntrons,CpGs.gr)
  
  # get relative position inside of the exon regardin exon width 
  rel_positions_CpG.do = abs(start(Domerge.CpG[,1]) - start(Domerge.CpG[,2]))
  
  # get negative strand 
  neg.do = which(strand(Domerge.CpG[,1]) == '-')
  rel_positions_CpG.do[neg.do] = end(Domerge.CpG[neg.do,1]) - start(Domerge.CpG[neg.do,2])
  
  # +1 because else there is a 0 position which is not possible
  rel_positions_CpG.do = rel_positions_CpG.do + 1
  
  # unique introns for matches found by GenomicRanges
  uniqueEntries.CpG.do = unique(Domerge.CpG[,1])
  
  # matrix containing bit-wise-vector for downstream intron regions 
  x.do.m = matrix(-1,nrow=length(uniqueEntries.CpG.do), ncol=inlength)
  
  # set 0 in position where you have a CpG 
  for(i in 1:inlength){
    posFori = which(rel_positions_CpG.do == i)
    rows = which(uniqueEntries.CpG.do %in% Domerge.CpG[posFori,1])
    x.do.m[rows,i] = 0
    if(i%%10 == 0) print(i)
  }
  
  ### mCPG
  Domerge.mCpG = mergeByOverlaps(grDoIntrons,mCpGs.gr)
  
  # get relative position inside of the exon regardin exon width 
  rel_positions_mCpG.do = abs(start(Domerge.mCpG[,1]) - start(Domerge.mCpG[,2]))
  
  # get negative strand 
  neg.do = which(strand(Domerge.mCpG[,1]) == '-')
  rel_positions_mCpG.do[neg.do] = end(Domerge.mCpG[neg.do,1]) - start(Domerge.mCpG[neg.do,2])
  
  # +1 because else there is a 0 position which is not possible
  rel_positions_mCpG.do = rel_positions_mCpG.do + 1
  
  # unique introns for matches found by GenomicRanges
  uniqueEntries.mCpG.do = unique(Domerge.mCpG[,1])
  
  # set 1 in position where you have a methylation 
  for(i in 1:inlength){
    posFori = which(rel_positions_mCpG.do == i)
    
    # other match before due to the structure of the matrix
    rows = match(Domerge.mCpG[posFori,1],uniqueEntries.CpG.do)
    x.do.m[rows,i] = 1
    if(i%%10 == 0) print(i)
  }
  
  # list of matches between upstream introns and methylation 
  # some introns are not listed 
  uniqueEntries = unique(c(uniqueEntries.CpG.do, uniqueEntries.mCpG.do))
  intronMatches.do.l = match(uniqueEntries, grDoIntrons)
  
  # run test of randomness for upstream introns 
  # to do the runs.test you need dichotomous data 
  runs_x.do.m = x.do.m
  runs_x.do.m[which(x.do.m == -1)] = 0
  runsTest.do.l = apply(runs_x.do.m, 1, runsTest)
  
  ################
  ### UPSTREAM ### 
  ################
  
  ### CPG
  Upmerge.CpG = mergeByOverlaps(grUpIntrons,CpGs.gr)
  
  # get relative position inside of the exon regardin exon width 
  rel_positions_CpG.up = abs(start(Upmerge.CpG[,1]) - start(Upmerge.CpG[,2]))
  
  # get negative strand 
  neg.up = which(strand(Upmerge.CpG[,1]) == '-')
  rel_positions_CpG.up[neg.up] = end(Upmerge.CpG[neg.up,1]) - start(Upmerge.CpG[neg.up,2])
  
  # +1 because else there is a 0 position which is not possible
  rel_positions_CpG.up = rel_positions_CpG.up + 1
  
  # unique introns for matches found by GenomicRanges
  uniqueEntries.CpG = unique(Upmerge.CpG[,1])
  
  # matrix containing bit-wise-vector for upstream intron regions 
  x.up.m = matrix(-1,nrow=length(uniqueEntries.CpG), ncol=inlength)
  
  # set 0 in position where you have a CpG 
  for(i in 1:inlength){
    posFori = which(rel_positions_CpG.up == i)
    rows = which(uniqueEntries.CpG %in% Upmerge.CpG[posFori,1])
    x.up.m[rows,i] = 0
    if(i%%10 == 0) print(i)
  }
  
  ### mCPG
  Upmerge.mCpG = mergeByOverlaps(grUpIntrons,mCpGs.gr)
  
  # get relative position inside of the exon regardin exon width 
  rel_positions_mCpG.up = abs(start(Upmerge.mCpG[,1]) - start(Upmerge.mCpG[,2]))
  
  # get negative strand 
  neg.up = which(strand(Upmerge.mCpG[,1]) == '-')
  rel_positions_mCpG.up[neg.up] = end(Upmerge.mCpG[neg.up,1]) - start(Upmerge.mCpG[neg.up,2])
  
  # +1 because else there is a 0 position which is not possible
  rel_positions_mCpG.up = rel_positions_mCpG.up + 1
  
  # unique introns for matches found by GenomicRanges
  uniqueEntries.mCpG = unique(Upmerge.mCpG[,1])
  
  # set 1 in position where you have a methylation 
  for(i in 1:inlength){
    posFori = which(rel_positions_mCpG.up == i)
    
    # other match before due to the structure of the matrix
    rows = match(Upmerge.mCpG[posFori,1],uniqueEntries.CpG)
    x.up.m[rows,i] = 1
    if(i%%10 == 0) print(i)
  }
  
  # list of matches between upstream introns and methylation 
  # some introns are not listed 
  uniqueEntries = unique(c(uniqueEntries.CpG, uniqueEntries.mCpG))
  intronMatches.up.l = match(uniqueEntries, grUpIntrons)
  
  # run test of randomness for upstream introns 
  # to do the runs.test you need dichotomous data 
  runs_x.up.m = x.up.m
  runs_x.up.m[which(x.up.m == -1)] = 0
  runsTest.up.l = apply(runs_x.up.m, 1, runsTest)
  
  #########################################
  ## create bp-wise methylation profile  ##
  #########################################
  
  profile.up.l = apply( x.up.m, 1, function(x) 
                                    return( paste0(as.character(x), collapse="|") ) )
  
  
  profile.do.l = apply( x.do.m, 1, function(x) 
                                    return( paste0(as.character(x), collapse="|") ) )
  
  profile.ex.l = apply( x.ex.mCpG.m, 1, function(x) 
                                    return( paste0(as.character(x), collapse="|") ) )
  
  profile.ex.ratio.l = apply( x.ex.ratio.m, 1, function(x) 
                                    return( paste0(as.character(x), collapse="|") ) )
  
  
  #################################
  ## calculate GC and C content  ##
  #################################
  
  GCandC_Content <-function(x) {
    x = DNAString(x)
    alf <- alphabetFrequency(x, as.prob=TRUE)
    return( c(sum(alf[c("G", "C")]), sum(alf[c("C")])) )
  }
  
  reference = BSgenome.Hsapiens.UCSC.hg19
  
  ################
  ### UPSTREAM ### 
  ################ 
  
  sequences = getSeq(reference, names=filtered_Upintrons.df[,1],   
                     start=as.numeric(filtered_Upintrons.df[,2]), 
                     end=as.numeric(filtered_Upintrons.df[,3]),
                     strand=as.character(filtered_Upintrons.df[,4]),
                     as.character=T)
  
  sequences.df = data.frame(sequences)
  UpGCCvalues = apply(sequences.df, 1, GCandC_Content)
  UpGCCvalues = matrix(unlist(UpGCCvalues), ncol=2, byrow=TRUE)
  print('GC and C uptream introns finished')
  
  ##################
  ### DOWNSTREAM ### 
  ##################
  
  sequences = getSeq(reference, names=filtered_Dointrons.df[,1],   
                     start=as.numeric(filtered_Dointrons.df[,2]), 
                     end=as.numeric(filtered_Dointrons.df[,3]),
                     strand=as.character(filtered_Dointrons.df[,4]),
                     as.character=T)
  
  sequences.df = data.frame(sequences)
  DoGCCvalues = apply(sequences.df, 1, GCandC_Content)
  DoGCCvalues = matrix(unlist(DoGCCvalues), ncol=2, byrow=TRUE)
  print('GC and C donwstream introns finished')
  
  #############
  ### Exons ### 
  #############
  
  sequences = getSeq(reference, names=chosenExons[,1],   
                     start=as.numeric(chosenExons[,2]), 
                     end=as.numeric(chosenExons[,3]),
                     strand=as.character(chosenExons[,4]),
                     as.character=T)
  
  sequences.df = data.frame(sequences)
  ExGCCvalues = apply(sequences.df, 1, GCandC_Content)
  ExGCCvalues = matrix(unlist(ExGCCvalues), ncol=2, byrow=TRUE)
  print('GC and C exons finished')
  
  ############################
  ## create feature matrix  ##
  ############################
  
  # offset profile for Introns
  min1Vector = rep('-1|', (inlength-1))
  min1Vector = paste0(min1Vector, collapse='')
  min1Vector = paste0(min1Vector, '-1', collapse='')
  
  # offset profile for Exons
  min1VectorEx = rep('-1|', (exonBins-1))
  min1VectorEx = paste0(min1VectorEx, collapse='')
  min1VectorEx = paste0(min1VectorEx, '-1', collapse='')

  # more for deep learning 
  features.df = matrix(nrow=nrow(Ematrix), ncol=24)
  features.df[,1] = Ematrix[,1]
  features.df[,2] = chosenExons[,9]
  features.df[,3] = chosenExons[,10]
  features.df[,4] = rep(cs, nrow(features.df))
  features.df[,5] = Ematrix[,4]
  features.df[,6] = rep(min1Vector, nrow(features.df))
  features.df[,7] = rep(min1Vector, nrow(features.df))
  features.df[intronMatches.up.l,6] = profile.up.l
  features.df[intronMatches.do.l,7] = profile.do.l
  
  features.df[,8] = UpGCCvalues[,1]
  features.df[,9] = DoGCCvalues[,1]
  features.df[,10] = ExGCCvalues[,1]
  features.df[,11] = UpGCCvalues[,2]
  features.df[,12] = DoGCCvalues[,2]
  features.df[,13] = ExGCCvalues[,2]
  features.df[,14] = rep(min1VectorEx, nrow(features.df))
  features.df[exonMatches.mCpGs.l,14] = profile.ex.l
  features.df[intronMatches.up.l,15] = runsTest.up.l
  features.df[intronMatches.do.l,16] = runsTest.do.l
  features.df[,17] = UpImatrix[,6]
  features.df[,18] = Ematrix[,6]
  features.df[,19] = DoImatrix[,6]
  features.df[,20] = rep(min1VectorEx, nrow(features.df))
  features.df[exonMatches.CpGs.l,20] = profile.ex.ratio.l
  features.df[,21] = as.character(chosenExons[,1])
  features.df[,22] = as.numeric(chosenExons[,2])
  features.df[,23] = as.numeric(chosenExons[,3])
  features.df[,24] = as.character(chosenExons[,4])
  colnames(features.df) = c("ID", 'PSI', 'MAP_PSI', 'In/Out', 'methEx', 
                            'UpProfile', 'DoProfile', 'UpGC', 'DoGC', 'ExGC',
                            'UpC', 'DoC', 'ExC', 'ExProfile', 'UpRun', 'DoUp',
                            'UpmCpg/Cpg', 'ExmCpG/CpG', 'DomCpg/Cpg', 'ExRatioProf',
                            'ChrExon', 'StartExon', 'EndExon', 'StrandExon')
  
  # check if bitwise vectors have the same number of methylation as in datamatrix before
  num.CpGsUp = nchar(gsub('-1|\\|','',features.df[,6]))
  num.CpGsDo = nchar(gsub('-1|\\|','',features.df[,7]))
  
  num.mCpGsUp = nchar(gsub('-1|\\||0','',features.df[,6])) 
  num.mCpGsDo = nchar(gsub('-1|\\||0','',features.df[,7])) 
  
  if ( length(which(num.mCpGsUp != UpImatrix[,2])) != 0 ) {
    print('[ERROR] number of methylations for the upstream region not correct in the feature vector')
  }
  
  if ( length(which(num.mCpGsDo != DoImatrix[,2])) != 0 ) {
    print('[ERROR] number of methylations for the downstream region not correct in the feature vector')
  }
  
  # in the script methPSI I changed NAs to -1, 
  # for the checkup I need to change them to a ration of 0 
  changedNAsUp = UpImatrix[,6]
  changedNAsUp[which(changedNAsUp == -1)] = 0
  
  if ( length(which((num.mCpGsUp/num.CpGsUp) != changedNAsUp)) != 0 ) {
    print('[ERROR] number of mCpG/CpG for the upstream region not correct in the feature vector')
  }
  
  # in the script methPSI I changed NAs to -1, 
  # for the checkup I need to change them to a ration of 0 
  changedNAsDo = DoImatrix[,6]
  changedNAsDo[which(changedNAsDo == -1)] = 0
  
  if ( length(which((num.mCpGsDo/num.CpGsDo) != changedNAsDo)) != 0 ) {
    print('[ERROR] number of mCpG/CpG for the downstream region not correct in the feature vector')
  }
  
  num.featuresEx = lapply(features.df[,14], function(x) {
                                             v = strsplit(x, split='[|]')[[1]]
                                             v = as.numeric(unlist(v))
                                             return(sum(v) / exonBins)
                                             })
  num.featuresEx = unlist(num.featuresEx)
  num.featuresEx[which(num.featuresEx == -1)] = 0
  
  differences = which(round(num.featuresEx, digits=6) != round(Ematrix[,4], digits=6) )
  
  if ( length(differences) != 0 ) {
    print('[ERROR] number of methylations for the exons not correct in the feature vector')
    print(num.featuresEx[differences])
    print(Ematrix[differences,4])
  }
  
  # check if mCpg/Cpg profile matches mCpG/Cpg ratio from matrix before
  check.Ex = lapply(features.df[,20], function(x) {
    v = strsplit(x, split='[|]')[[1]]
    v = as.numeric(unlist(v))
    v = v[which(v != -1)]
    w = 1/v
    w[which(is.infinite(w))] = 1
    return(c(length(which(v != 0)), sum(w)))
  })
  
  check.Ex.df = matrix(unlist(check.Ex), ncol=2, byrow = TRUE)
  check.ratio.l = (check.Ex.df[,1]/check.Ex.df[,2])
  check.ratio.l[which(is.na(check.ratio.l))] = 0
  
  # in the script methPSI I changed NAs to -1, 
  # for the checkup I need to change them to a ration of 0 
  changedNAsEx = Ematrix[,6]
  changedNAsEx[which(changedNAsEx == -1)] = 0
  
  if ( length(which(check.ratio.l != changedNAsEx)) != 0 ) {
    
    # pick wrong entries 
    wrongEntires.l = which(check.ratio.l != changedNAsEx)
    
    # look if wrong entries matches with IDs with has double or more CpGs or mCpGs
    # in a bin ---> this means everything is right but it can be hard to check 
    if ( length(which(!wrongEntires.l %in% doubleEntryMatchIDs)) != 0) {
      print('[ERROR] mCpG/CpG ratio for the exon region not correct')
    }
  }
  
  # checkpoint
  checkpoint1 = paste0('data/features/features.',cs,'.',shortening,'.',shift,'.',inlength,'.',celltype,'.tsv')
  if ( file.exists(checkpoint1) ) {
    file.remove(checkpoint1)
  }
  write.table(features.df, file = checkpoint1, sep='\t', append=F)
  
  print('[FINISH]')
}

#######################
## Execute Pipeline  ##
#######################

# cs, shortening, shift, inlength, exonBins
pipe(TRUE, FALSE, 0.0, 500, 50)
pipe(FALSE, FALSE, 0.0, 500, 50)
