
library(GenomicRanges)
library(plyr)
library(stringr)      # library two split one column into multiple cols

args = commandArgs(trailingOnly = TRUE)
celltype = args[1]

# read read counts for exons
if ( celltype == 'IMR90' ) {
  tablecounts = read.delim("/tmp/new_counts_IMR90.tsv", header=F)
}

if ( celltype == 'Gm12878' ) {
  tablecounts = read.delim("tmp/new_counts_Gm12878.tsv", header=F)
}

if ( celltype == 'H1hesc' ) {
  tablecounts = read.delim("tmp/new_counts_H1hesc.tsv", header=F)
}

################
##  Filter A  ##
################

# order matrix
ordered_tablecounts = tablecounts[order(tablecounts[,5], as.numeric(tablecounts[,2])),]

# Filter out Genes which has less than 3 exons (no triplett formable).
filterA_count = ordered_tablecounts

# get uniqe starts for each gene and count occurrence of gene per start position
uniqueStart = unique(filterA_count[,c(2,5)])
genStartsCount = count(uniqueStart[,2])
genOnlyOneStart = genStartsCount[which(genStartsCount[,2] < 3),1]

# remove genes with less than 3 starting position
filterA_count = filterA_count[-which(filterA_count[,5] %in% genOnlyOneStart),]

# get uniqe ends for each gene and count occurrence of gene per end position
uniqueEnd = unique(filterA_count[,c(3,5)])
genEndsCount = count(uniqueEnd[,2])
genOnlyOneEnd = genEndsCount[which(genEndsCount[,2] < 3),1]

# remove gene with less than 3 end position
filterA_count = filterA_count[-which(filterA_count[,5] %in% genOnlyOneEnd),]

################
##  Filter B  ##
################

# filter out exons which has no read coverage whats so ever
filterB_count = filterA_count

if ( length(which(filterB_count[,8] == 0)) != 0 ) { 
  # filterB_count = filterB_count[-which(filterB_count[,8] == 0),]
  print('[ERROR] there are still zero covered reads in data')
}

# and filter out exons which are two small (<10bp)
filterB_count = filterB_count[-which( (filterB_count[,3] - filterB_count[,2]) < 50 ),]

################
##  Filter C  ##
################

# filter out entries for exons which are exactly the same (same chr,start,end,strand,gene)
filterC_count = filterB_count
filterC_count = filterC_count[-which(duplicated(filterC_count[,1:5])),]

################
##  Filter D  ##
################

# filter out exons which overlaps with other genes

# collect first and last exons
firstexons = which(duplicated(ordered_tablecounts[,5]) == F )
lastexons = which(duplicated(ordered_tablecounts[,5], fromLast=T) == F )

genes = matrix(nrow=length(firstexons), ncol=5)
colnames(genes) = c("chr", "start", "end", "strand", "ID")
genes[,1] = as.character(ordered_tablecounts[firstexons,1])
genes[,2] = as.numeric(ordered_tablecounts[firstexons,2])
genes[,3] = as.numeric(ordered_tablecounts[lastexons,3])
genes[,4] = as.character(ordered_tablecounts[firstexons,4])
genes[,5] = as.character(ordered_tablecounts[firstexons,5])

# create GRanges objects for genes and filtered exons
genes.gr = GRanges(seqnames=genes[,1],
                   ranges=IRanges(start=as.numeric(genes[,2]),
                                  end=as.numeric(genes[,3])),
                   strand=genes[,4])

filterC_count.gr = GRanges(seqnames=filterC_count[,1],
                           ranges=IRanges(start=as.numeric(filterC_count[,2]),
                                          end=as.numeric(filterC_count[,3])),
                           strand=as.character(filterC_count[,4]))

# find overlaps betweens genes and exons 
overlaps = findOverlaps(genes.gr, filterC_count.gr)

#             queryHits subjectHits 
#             <integer>   <integer> 
#   1              1           1 
# 2              1           2 
# 3              1           3 
# 4              1           4 
# 5              1           5 
# ...          ...         ... 

# find duplicated hits (exons should only be in one gene)
removeExons = subjectHits(overlaps)[which(duplicated(subjectHits(overlaps)) == T)]

# find unique exons which should be removed
# length(removeExons) = 33312
removeExons = unique(removeExons)

# remove exons from filter list from before 
filterD_count = filterC_count
filterD_count = filterC_count[-removeExons,]
rownames(filterD_count) <- NULL

################
##  Filter E  ##
################

# paste together chromosome name and gene name
seqnames.l = apply(filterD_count[,c(1,5)], 1, function(x) paste(x, collapse=':'))

# filter out exons which olverap with exons inside a gene
# this ensures that we have a real consecutive triplett
# this can happen due to different versions of an exons (alternative spliced exon)
filterD_count.gr = GRanges(seqnames=seqnames.l,
                        ranges=IRanges(start=as.numeric(filterD_count[,2]),
                                       end=as.numeric(filterD_count[,3])),
                        strand=as.character(filterD_count[,4]))

# reduce filterD_count.gr (overlapping exons within a gene are merged together
# to one big region)
reduced_filterD_count.gr = reduce(filterD_count.gr)

# find overlaps between the reduced gr and the whole gr
overlapsD = findOverlaps(reduced_filterD_count.gr, filterD_count.gr)

# get subject hits
hits = subjectHits(overlapsD)

# get read IDs of the subjectHits
allreads = filterD_count[hits,9]

# which regions are not duplicated (first entry)
notDuplicated = which(duplicated(queryHits(overlapsD)) == F )

# which regions are duplicated (second, third ... entry)
Duplicated = which(duplicated(queryHits(overlapsD)) == T )

# create new matrix out of reduced_filterD_count.gr
filterE_count = as.data.frame(reduced_filterD_count.gr)

# change chromosome name and add Gene ID column
chrGeneID = str_split_fixed(as.character(filterE_count[,1]), ":", 2)

filterE_count[,1] = chrGeneID[,1]
filterE_count = cbind(filterE_count, chrGeneID[,2])
colnames(filterE_count)[6] = 'GeneID'

# add read IDs
filterE_count = cbind(filterE_count, as.character(allreads[notDuplicated]))
filterE_count[,7] = as.character(filterE_count[,7])

# add reads of the duplicated entries
rows = queryHits(overlapsD)[Duplicated]
unique_rows = unique(rows)

# create array which holds reads which need to be added 
addedreads = allreads[Duplicated]

print(length(unique_rows))
mix.reads.l = as.character(c(numeric(length(unique_rows))))

for ( i in 1:length(unique_rows) ){
  reads.row = filterE_count[unique_rows[i],7]
  reads.row = unlist(strsplit(as.character(reads.row), split=';'))
  
  # reads of duplicated row which I want to add
  add.entries = which(rows == unique_rows[i])
  
  # read I want to add to entry 
  addedreads.entry = addedreads[add.entries]
  
  # split theses reads 
  added.reads.l = unlist(strsplit(as.character(addedreads.entry), split=';'))
  
  # mix them together with original entry
  mix.reads = c(reads.row,added.reads.l)
  mix.reads = unique(mix.reads)
  
  # paste them together 
  mix.reads.l[i] = paste(mix.reads, collapse=';')
  if(i %% 1000 == 0){print(i)}
}

filterE_count[unique_rows,7] = mix.reads.l

################
##  Filter F  ##
################

# same as filter A in case of the filtering before
filterF_count = filterE_count

uniqueStart = unique(filterF_count[,c(2,5)])
genStartsCount = count(uniqueStart[,2])
genOnlyOneStart = genStartsCount[which(genStartsCount[,2] < 3),1]

if ( length(which(filterF_count[,5] %in% genOnlyOneStart)) != 0 ) {
  filterF_count = filterF_count[-which(filterF_count[,5] %in% genOnlyOneStart),]
}

uniqueEnd = unique(filterF_count[,c(3,5)])
genEndsCount = count(uniqueEnd[,2])
genOnlyOneEnd = genEndsCount[which(genEndsCount[,2] < 3),1]

if ( length(which(filterF_count[,5] %in% genOnlyOneEnd)) != 0 ) { 
  filterF_count = filterF_count[-which(filterF_count[,5] %in% genOnlyOneEnd),]
}

#############
##  Count  ##
#############

# Get read count for junction of an inclusion and exclusion of middle exon of a 
# exon triplet C1-A -C2. Where C1-A + C2-A = number of reads for an inclusion of A. 
# Hence C1 - C2 = number of reads for an exclusion of A. 
counts = filterF_count
counts = counts[-4]
counts = counts[-6]
for ( i in 1:2 ){ counts = cbind(counts,numeric(nrow(filterF_count))) }
colnames(counts) = c('chr','start','end','strand','genID','Ninc','Nexc') 

# split the readIDs for each exon 
readIDs = strsplit(as.character(filterF_count[,7]), split=';')

# sliding window function that takes 3 exons and evaluates Ninc and Nexc
slideCountFunction = function(x) {
  c1 = x[[1]]
  a = x[[2]]
  c2 = x[[3]]

  Ninc = length(intersect(c1, a)) + length(intersect(a, c2))
  
  # intersect(c1,c2,a) do not count reads which spans over all three exons
  Nexc = length(intersect(c1, c2)) - length( intersect(intersect(c1,c2),a) )

  return( c(Ninc , Nexc) )
}

# fast version for a sliding window over a sequence 
wapply <- function(x, width, by = NULL, FUN = NULL, ...)
{
  FUN <- match.fun(FUN)
  if (is.null(by)) by <- width
  
  lenX <- length(x)
  SEQ1 <- seq(1, lenX - width + 1, by = by)
  SEQ2 <- lapply(SEQ1, function(x) x:(x + width - 1))
  
  OUT <- lapply(SEQ2, function(a) FUN(x[a], ...))
  OUT <- base:::simplify2array(OUT, higher = TRUE)
  return(OUT)
}

# start sliding window 
start = Sys.time()
countarray = wapply(readIDs, width=3, by=1, FUN=slideCountFunction)
end = Sys.time()
print(end-start)

# fill Ninc and Nexc into matrix
counts[2:(nrow(counts)-1),6] = countarray[1,]
counts[2:(nrow(counts)-1),7] = countarray[2,]

######################
##  total coverage  ##
######################

# Total coverage of each junction C1-A-C2 is (Ninc/2) + Nexc
counts2 = counts
counts2 = cbind(counts,numeric(nrow(counts)))
counts2[,8] = (counts2[,6] / 2) + counts2[,7]
colnames(counts2)[8] = 'total'

#############################
##  First Estimate of PSI  ##
#############################

# First estimate of PSI is related to paper from Xiong et. al. 2014 'The human splicing 
# code reveals new insights into the genetic determinant of disease'
counts3 = counts2
counts3 = cbind(counts3,numeric(nrow(counts3)))
counts3[,9] = (counts2[,6] / 2) / counts3[,8]
colnames(counts3)[9] = "apriori_PSI" 

##############################
##  Second Estimate of PSI  ##
##############################

# Get prior probability for an inclusion from Beta distribution:
# assuming inclusion and exclusion of an exon are equally likely (50:50).
# So building up a beta distribution with Beta(10,10) and include the
# read counts for inclusion and exclusion.
# If you set the constants (10,10) even higher you force the prior to be more like 0.5.
offset = 10
a = ((counts2[,6] / 2) + offset)
b = (counts3[,7] + offset)

# With beta distribution and likelihood of a bernoulli distribution you get an 
# estimate for the maximum a posterior probability for the number of reads for 
# an inclusion. The short formula below can be derived from a bayesian approach
# p(u,X) = p(X,u) * p(u)
# Where p(u) the pior probability distribution is our beta distribution.
# Where p(X,u) is our likelihood from the bernoulli distribution: 
# p(X,u) = product [u^(Ninc/2) * (1-u)^Nexc].
# Taking log: log(p(u,X)) = log(p(X,u)) + log(p(u)) 
# and maximize by setting first derivative to 0 you will end up in the formula below.
map = ( (counts2[,6]/2) + a - 1 )  / ( counts3[,8] + b + a - 2 ) 

# bind map to matrix and write file
counts4 = counts3
counts4 = cbind(counts4,numeric(nrow(counts4)))
counts4[,10] = map

colnames(counts4)[10] = 'max_aposterior_PSI' 

checkpoint1 = paste0('data/PSI.10.',celltype,'.tsv')

if ( file.exists(checkpoint1) )
  file.remove(checkpoint1)

write.table(counts4, file = checkpoint1, sep='\t')

print('[FINISH]')




