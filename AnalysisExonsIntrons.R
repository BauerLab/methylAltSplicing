
library(plyr)
library(GenomicRanges)
library(ggplot2)
library(Rmisc)
library(BSgenome.Hsapiens.UCSC.hg19)
library(bsseq)
library(prodlim)      # for row matching
library(mixtools)     # to get mixture gaussian approximation
library(stringr)      # to count substring in a string
library(data.table)   # for quick import of data

#cs = F
#shortening = F
#shift = 0.0
#inlength = 500
pipe = function(cs, shortening, shift, inlength, zero, celltype, plotcellnames) {
  
  # write to log-file
  log.file = paste0("plots/log/",cs,'.',shortening,'.txt')
  file.create(log.file)
  sink(file=log.file)
  
  ########################
  ## Get cufflink data  ##
  ########################
  
  # list of replica
  patients.l = c('1','2')

  # big matrix
  cufflink.gene.df = matrix(ncol=7)
  colnames(cufflink.gene.df) = c("chr","start","end","strand","geneID","expression","meanexpr")

  # check if big matrix file already exists if not then create it
  checkpoint0 = paste0('data/methPSI.cuff.bigtable.',celltype,'.tsv')
  
  # for multiplott of ggplot
  plots <- list()
  
  postscript(file = paste0('plots/log/histogramCov',celltype,'.eps'), width=1400, height=800,  family="Times")
  par(cex = 2.0)
  
  if ( file.exists(checkpoint0) == F ) {

    # go over all patients
    for (i in 1:length(patients.l) ) {

      print(patients.l[i])

      # read in cufflink data for patient
      cufflink.exons.df <- read.delim(paste0('data/',celltype,'_',patients.l[i],'/cuff.exons.',celltype,'_',patients.l[i],'.tsv'), header=F)
      colnames(cufflink.exons.df) = c("chr","start","end","score","strand","geneID","trID","FPKM","cov")

      # make ggplot to see the distribution of the coverage
      quantiles_cov = quantile(cufflink.exons.df[,9])
      p <- ggplot(cufflink.exons.df, aes(x=cov)) + 
              theme(text=element_text(family="Times")) + 
              geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                             binwidth=.01,
                             colour="black", fill="grey") +
              xlab('Read Coverage') +
              ylab('Density')  +
              scale_x_continuous(trans='log10', breaks=c(quantiles_cov[2], quantiles_cov[4]),
                                 labels = c(round(quantiles_cov[2],digits=1), round(quantiles_cov[4],digits=1))) +
              # ggtitle(paste0('Distribution of Read Coverage for each Exon of Sample ',i)) +
              ggtitle(paste0(plotcellnames,' Replica ',i)) +
              geom_vline(xintercept = quantiles_cov[2], color = 'red', size = 1) +
              geom_vline(xintercept = quantiles_cov[4], color = 'red', size = 1) +
              theme(text = element_text(size=25), axis.title.y=element_text(margin=margin(0,20,0,0)),
                    axis.title.x=element_text(margin=margin(20,0,0,0)), plot.margin = unit(c(1,1,1,1),'cm'))
      
      plots[[i]] <- p 

      ############################
      ## Prepare cufflink data  ##
      ############################

      # order exons based on geneID and start position
      ordered_cufflink.exons.df = cufflink.exons.df[order(cufflink.exons.df[,6], as.numeric(cufflink.exons.df[,2])),]

      # unique genes from cufflink data
      uniqueGenes =  unique(ordered_cufflink.exons.df[,6])

      # prepare gene matrix
      patient.gene.df = matrix(ncol=7, nrow=length(uniqueGenes))
      colnames(patient.gene.df) = c("chr","start","end","strand","geneID","expression","meanexpr")

      # collect first and last exons
      firstexons = which(duplicated(ordered_cufflink.exons.df[,6]) == F )
      lastexons = which(duplicated(ordered_cufflink.exons.df[,6], fromLast=T) == F )

      # fill matrix
      patient.gene.df[,1] = ordered_cufflink.exons.df[firstexons,1]
      patient.gene.df[,2] = ordered_cufflink.exons.df[firstexons,2]
      patient.gene.df[,3] = ordered_cufflink.exons.df[lastexons,3]
      patient.gene.df[,4] = as.character(ordered_cufflink.exons.df[lastexons,5])
      patient.gene.df[,5] = as.character(ordered_cufflink.exons.df[lastexons,6])
      
      row.matches = match(ordered_cufflink.exons.df[,6], uniqueGenes)

      # chr    start       end         strand geneID expression
      # [1,] "chr1" "760728"    "762886"    "1"    "1"    "1"
      # [2,] "chr1" "901877"    "911245"    "3"    "2"    "0"
      # [3,] "chr1" "6472478"   "6484730"   "1"    "3"    "0"
      # [4,] "chr1" "79353619"  "79355169"  "2"    "4"    "1"
      # [5,] "chr6" "150070579" "150132044" "3"    "5"    "1"
      # [6,] "chr6" "150139934" "150186199" "1"    "6"    "1"

      # find out if gene is expressed or silenced
      for (j in 1:length(uniqueGenes) ){

        exons.l = which(row.matches == j)

        # get score FPKM and coverage for the exon
        score.l = ordered_cufflink.exons.df[exons.l,4]
        FPKM.l = ordered_cufflink.exons.df[exons.l,8]
        cov.l = ordered_cufflink.exons.df[exons.l,9]

        # look if one of the exons has significant scores which
        # would identify the gene as expressed
        if ( max(score.l) > 1 && max(FPKM.l) > 1.0  && max(cov.l) > 20.0 ){
          patient.gene.df[j,6] = 1
          patient.gene.df[j,7] = mean(cov.l)
        } else {
          patient.gene.df[j,6] = 0
        }

        if(j%%1000 == 0){print(j)}
      }

      # bind the genes observed for each patient to the big matrix
      cufflink.gene.df = rbind(cufflink.gene.df, patient.gene.df)

    }

    # change chromosome names from numbers to chr1 etc.
    cufflink.gene.df[,1] = paste0('chr',cufflink.gene.df[,1])
    
    # write table if the file not already exist
    write.table(cufflink.gene.df, file = checkpoint0, sep='\t', append=F)
  } else {
    # make only plots
    for (i in 1:length(patients.l) ) {
      
      print(patients.l[i])
      
      # read in cufflink data for patient
      cufflink.exons.df <- read.delim(paste0('data/',celltype,'_',patients.l[i],'/cuff.exons.',celltype,'_',patients.l[i],'.tsv'), header=F)
      colnames(cufflink.exons.df) = c("chr","start","end","score","strand","geneID","trID","FPKM","cov")
      
      # make ggplot to see the distribution of the coverage
      quantiles_cov = quantile(cufflink.exons.df[,9])
      
      p <- ggplot(cufflink.exons.df, aes(x=cov)) +
              theme(text=element_text(family="Times")) + 
              geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                             binwidth=.01,
                             colour="black", fill="grey") +
              xlab('Read Coverage') +
              ylab('Density')  +
              scale_x_continuous(trans='log10', breaks=c(quantiles_cov[2], quantiles_cov[4]),
                                 labels = c(round(quantiles_cov[2],digits=1), round(quantiles_cov[4],digits=1))) +
              ggtitle(paste0(plotcellnames, ' Replica ',i)) +
              geom_vline(xintercept = quantiles_cov[2], color = 'red', size = 1) +
              geom_vline(xintercept = quantiles_cov[4], color = 'red', size = 1) +
              theme(text = element_text(size=25), axis.title.y=element_text(margin=margin(0,20,0,0)),
                    axis.title.x=element_text(margin=margin(20,0,0,0)), plot.margin = unit(c(1,1,1,1),'cm'))

      plots[[i]] <- p 
    }
  }
  multiplot(plotlist = plots, cols = 2)
  dev.off()
  
  # read table
  cufflink.gene.df = read.delim(checkpoint0, header=T)

  # remove first row (it has NAs in it)
  cufflink.gene.df = cufflink.gene.df[-1,]

  # change strand of genes '.' to '*' because of GenomeRanges
  cufflink.gene.df[,4] = as.character(cufflink.gene.df[,4])
  cufflink.gene.df[which(!cufflink.gene.df[,4] %in% c('+','-')),4] = '*'
  
  ###############
  ## Filter 1  ##
  ###############
  
  # remove first and last exon of each gene
  filter1_psi.df = psi.df
  firstexons = which(duplicated(filter1_psi.df[,5]) == F )
  lastexons = which(duplicated(filter1_psi.df[,5], fromLast=T) == F )
  
  rowNumExons = c(firstexons, lastexons)
  
  # get only unique rows
  rowNumExons = unique(rowNumExons)
  
  # remove exons
  filter1_psi.df = filter1_psi.df[-rowNumExons,]
  rownames(filter1_psi.df) <- NULL
  
  #################
  ## Group Exons ##
  #################
  
  if( length(which(is.na(filter1_psi.df[,10]))) != 0 ) {
    print('[INFO] some MAP are NA')
    filter1_psi.df[which(is.na(filter1_psi.df[,10])),10] = 0
  }
  
  if ( celltype != 'IMR90' && celltype != 'Gm12878' && celltype != 'H1hesc' ) {
  
    # set seed for the EM algorithm
    set.seed(123)
    
    # model mixture of two gaussian distribution on the whole population 
    mixmdl = normalmixEM(as.numeric(filter1_psi.df[,10]))
    
    # get the mean 
    mu = mixmdl$mu
    
    # get the standard deviation 
    sigma = mixmdl$sigma
    
    # get the raw ordered data 
    ordered_x = unique(mixmdl$x[order(mixmdl$x, decreasing=F)])
      
    # get first distribution 
    firstDnorm.df = as.data.frame(list(ordered_x, dnorm(ordered_x, mean=mu[1], sd=sigma[1])))
    firstDnorm.df[,2] = log(firstDnorm.df[,2] + 1)
    colnames(firstDnorm.df) = c('x','density')
    
    # get second distribution 
    secondDnorm.df = as.data.frame(list(ordered_x, dnorm(ordered_x, mean=mu[2], sd=sigma[2])))
    secondDnorm.df[,2] = log(secondDnorm.df[,2] + 1)
    colnames(secondDnorm.df) = c('x','density')
    
    # make ggplot to see the distribution of MAP_PSI + mixture model 
    postscript(file = 'plots/log/histogramMAP_PSI.eps', width = 600, height = 600,  family="Times")
    print(ggplot(filter1_psi.df, aes(x=max_aposterior_PSI)) + 
      theme(text=element_text(family="Times")) +  
      geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                     binwidth=.05,
                     colour="white", fill="grey") +
                     xlab('Maximum A Posterior Probability of Inclusion') +
                     ylab('Log Density') +
      scale_y_continuous(trans='log1p', breaks=c(0,1,3,4,5,9,10)) +
      geom_line(data=firstDnorm.df, aes(x=x, y=density), size=1.0, color='red') +
      geom_line(data=secondDnorm.df, aes(x=x, y=density), size=1.0, color='black', linetype="dotted") +
      geom_vline(xintercept = (mu[1] + sigma[1]), color = 'red', size = 1) + 
      geom_vline(xintercept = (mu[2] - sigma[2]), color = 'black', size = 1, linetype="dotted") +
      annotate("text", x = c((mu[1] + sigma[1]) - 0.04, (mu[2] - sigma[2]) + 0.04), 
               y = c(-0.05,-0.05), 
               label = c(round((mu[1] + sigma[1]),3), round((mu[2] - sigma[2]),3))) +
      ggtitle(plotcellnames) +
        theme(text = element_text(size=25), axis.title.y=element_text(margin=margin(0,20,0,0)),
              axis.title.x=element_text(margin=margin(20,0,0,0)), plot.margin = unit(c(1,1,1,1),'cm'))
    )
    dev.off()
    
    if (cs) {
      thres = mu[2] - sigma[2]
    } else {
      thres = mu[1] + sigma[1]
    }
  } else {

    postscript(file = paste0('plots/log/histogramMAP_PSI_',celltype,'.eps'), width = 600, height = 600,  family="Times")
    print(ggplot(filter1_psi.df, aes(x=max_aposterior_PSI)) +
            theme(text=element_text(family="Times")) + 
            geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                           binwidth=.05,
                           colour="white", fill="grey") +
            xlab('Maximum A Posterior Probability of Inclusion') +
            ylab('Log Density') +
            scale_y_continuous(trans='log1p', breaks=c(0,1,3,4,5,9,10)) +
            scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
            geom_vline(xintercept = 0.6, color = 'red', size = 1) + 
            geom_vline(xintercept = 0.8, color = 'red', size = 1) +
            ggtitle(plotcellnames) +
            theme(text = element_text(size=25), axis.title.y=element_text(margin=margin(0,20,0,0)),
                  axis.title.x=element_text(margin=margin(20,0,0,0)), plot.margin = unit(c(1,1,1,1),'cm'))
    )
    dev.off()

    apirori.df <- as.data.frame(filter1_psi.df[,9])
    colnames(apirori.df) <- 'Priori'
    apirori.df[which(is.na(apirori.df[,1])),1] = 0.0
    
    postscript(file = paste0('plots/log/histogramPSI_',celltype,'.eps'), width = 600, height = 600,  family="Times")
    print(ggplot(apirori.df, aes(x=Priori)) + 
            theme(text=element_text(family="Times")) + 
            geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                           binwidth=.05,
                           colour="white", fill="grey") +
            xlab('Probability of Inclusion') +
            ylab('Log Density') +
            scale_y_continuous(trans='log1p', breaks=c(0,1,3,4,5,9,10)) +
            #scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
            #geom_vline(xintercept = 0.6, color = 'red', size = 1) + 
            #geom_vline(xintercept = 0.8, color = 'red', size = 1) +
            ggtitle(plotcellnames) +
            theme(text = element_text(size=25), axis.title.y=element_text(margin=margin(0,20,0,0)),
                  axis.title.x=element_text(margin=margin(20,0,0,0)), plot.margin = unit(c(1,1,1,1),'cm'))
    )
    dev.off()
     
    if (cs) {
      thres = 0.8
    } else {
      thres = 0.4
    }
  
  }

  inclusions = filter1_psi.df[which(as.numeric(filter1_psi.df[,10]) >= thres),]
  exclusions = filter1_psi.df[which(as.numeric(filter1_psi.df[,10]) <= thres),]
  
  ################################################################
  ##  Remove Exons belonging to chromosome X, Y and M  Filter 2 ##
  ################################################################
  
  # find rows
  removeXYM.inc.l = c(which(inclusions[,1] == 'chrX'), which(inclusions[,1] == 'chrY'), which(inclusions[,1] == 'chrM'))
  removeXYM.exc.l = c(which(exclusions[,1] == 'chrX'), which(exclusions[,1] == 'chrY'), which(exclusions[,1] == 'chrM'))
  
  # remove rows
  removeXYM_inclusions = inclusions[-removeXYM.inc.l,]
  removeXYM_exclusions = exclusions[-removeXYM.exc.l,]
  
  chosenExons = c()
  if(cs) {
    chosenExons = removeXYM_inclusions
  } else {
    chosenExons = removeXYM_exclusions
  }
  
  # remove exons where Ninc = 0 and Nexc = 0 
  remove_exons.l = intersect(which(chosenExons[,6] == 0 ), which(chosenExons[,7] == 0 ))
  
  if ( length(remove_exons.l) != 0 ) {
    chosenExons = chosenExons[-remove_exons.l,]
  }
  
  rownames(chosenExons) = c(1:nrow(chosenExons))
  
  print(paste('[NOTE] before cufflink filtering number of exons =',nrow(chosenExons)))
  
  ####################################
  ##  Check with Cufflink Filter 3  ##
  ####################################

  chosenExons.gr = GRanges(seqnames=as.character(chosenExons[,1]),
                            ranges=IRanges(start=as.numeric(chosenExons[,2]),
                                           end=as.numeric(chosenExons[,3])),
                            strand=as.character(chosenExons[,4]))


  cufflink.gr = GRanges(seqnames=as.character(cufflink.gene.df[,1]),
                        ranges=IRanges(start=as.numeric(cufflink.gene.df[,2]),
                                       end=as.numeric(cufflink.gene.df[,3])),
                        strand=as.character(cufflink.gene.df[,4]))

  # included/excluded exons have to be inside of a region the cufflink data
  overlaps = findOverlaps(chosenExons.gr, cufflink.gr, type='within')

  # exon hits
  exonhits = queryHits(overlaps)

  # cufflink region hits
  cufflinkhits = subjectHits(overlaps)

  # unique exons
  uniqe_exonhits = unique(exonhits)

  # list holds the bitwise vector from cufflink
  # 0 = silenced, 1 = expressed
  expression.l = numeric(length(uniqe_exonhits))

  # probability of expression (number of 1 divided by total num observations)
  prob.expression.l = numeric(length(uniqe_exonhits))
  
  # mean expression
  # mean.expression.l = numeric(length(uniqe_exonhits))

  # go over all unqie exons
  for ( i in 1:length(uniqe_exonhits) ){
    # pick entries for each individual exons in the overlaps
    entries = which(exonhits == uniqe_exonhits[i] )

    # get cufflink data and look if exon is identified as expressed or not
    expression = cufflink.gene.df[cufflinkhits[entries],6]

    # collapse expression values into a string
    expression.l[i] = paste0(expression, collapse='')

    # get probability of expression
    prob.expression.l[i] = mean(expression)
    if(i%%1000 == 0){print(i)}
    
    # get mean expression
    # mean.expression.l[i] = cufflink.gene.df[cufflinkhits[entries],7]
    
  }

  # bind to matrix cufflink expression values
  # offset = -1 for exons which werent coffered by cufflink
  chosenExons = cbind(chosenExons, -1)
  colnames(chosenExons)[11] = 'CuffExp'

  # bind to matrix cufflink expression values
  # offset = -1 for exons which werent coffered by cufflink
  chosenExons = cbind(chosenExons, -1)
  colnames(chosenExons)[12] = 'MeanCuffExp'

  chosenExons[uniqe_exonhits, 11] = expression.l
  chosenExons[uniqe_exonhits, 12] = prob.expression.l

  # find genes which could be silenced
  maybe_silenced.genes.l = chosenExons[intersect(which(chosenExons[,12] <= 0.5), which(chosenExons[,12] != -1)),5]
  maybe_silenced.genes.l = unique(maybe_silenced.genes.l)

  # get max read coverage for the gene (exon with the maximum read coverage)
  # offset = -1 for gene which are not silenced and should not be removed
  chosenExons = cbind(chosenExons, -1)
  colnames(chosenExons)[13] = 'MaxReads'

  print(length(maybe_silenced.genes.l))
  
  # go over all genes which could be possibly be silenced
  for ( i in 1:length(maybe_silenced.genes.l) ){
    # get rows for individual gene
    rows = which(chosenExons[,5] %in% maybe_silenced.genes.l[i] )
    
    # check if all exons could be silenced or if some exons
    # did not overlap with cufflinks data and therefore have a -1 as offset
    if ( length(which(chosenExons[rows,12] == 1)) != 0 ) {
      chosenExons[rows,13] = -1
    } else {
      # get exons total number of reads and look for the maximum value
      max.reads = max(chosenExons[rows,8])
      chosenExons[rows,13] = max.reads
    }
    
    if(i%%1000 == 0){print(i)}
  }
    
  # remove whole gene if the threshold is not -1 (offset)
  remove.l = which(chosenExons[,13] != -1)
  
  # check if everything went right
  if ( length(which(chosenExons[remove.l,12] == 1.0)) != 0 ) {
    print('[ERROR] something went wrong with the cufflinks filtering')
  }

  print(paste('[NOTE] removing',length(remove.l),'samples because gene could be silenced'))

  if ( length(remove.l) != 0 ) {
    chosenExons = chosenExons[-remove.l,]
  }

  #######################################
  ## Get Expression Rate for each Exon ##
  #######################################
  
  expression.df <- as.data.frame(expression.df)
  
  m <- row.match(chosenExons[,1:4], expression.df[,1:4])
  chosenExons <- cbind(chosenExons, expression.df[m,ncol(expression.df)])
  
  #########################################
  ##  Remove bad sampled Exons Filter 4  ##
  #########################################
  
  # make ggplot to see the distribution of the average read depth for Ninc and Nexc
  t = 5
  postscript(file = paste0('plots/log/histogramMaxReads',celltype,'_',cs,'.eps'), width = 800, height = 800,  family="Times")
  print(ggplot(chosenExons, aes(x=total)) +
          theme(text=element_text(family="Times")) + 
          geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                         binwidth=.01,
                         colour="black", fill="grey") +
          xlab('Number of splitted Reads') +
          ylab('Density')  +
          scale_x_continuous(trans='log10', breaks = c(t, 100, 1000, 10000)) +
          ggtitle(plotcellnames) +
          geom_vline(xintercept = t, color = 'red', size = 1) +
          theme(text = element_text(size=25), axis.title.y=element_text(margin=margin(0,20,0,0)))
  )
  dev.off()
  
  # look how many bad samples you have
  badSamples = which(chosenExons[,8] < t)
  print(paste('[NOTE] remove total number of reads <',t,'=', length(badSamples)))
  
  if ( length(badSamples) != 0 ) {
    chosenExons = chosenExons[-badSamples,]
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
  
  # checkpoint 
  checkpoint1 = paste0('data/unfiltered/PSI.unfiltered.chosenExons.',cs,'.',celltype,'.tsv')
  if ( file.exists(checkpoint1) ) {
    file.remove(checkpoint1)
  }
  write.table(chosenExons, file = checkpoint1, sep='\t', append=F)
  
  # continue pipeline 
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
  
  
  #########################################################
  ## Remove introns which overlaps with exons (Filter 5) ##
  #########################################################
  
  overlaps1 = findOverlaps(all.exons.gr, Upintrons.gr)
  removeIntrons1 = unique(subjectHits(overlaps1))
  
  overlaps2 = findOverlaps(all.exons.gr, Dointrons.gr)
  
  removeIntrons2 = unique(c(removeIntrons1,unique(subjectHits(overlaps2))))
  
  # testing how many introns overlap with another version of the exon 
  chosenExons.test.gr = GRanges(seqnames=chosenExons[,1],
                                  ranges=IRanges(start=as.numeric(chosenExons[,2]),
                                                 end=as.numeric(chosenExons[,3])),
                                strand=chosenExons[,4])
  
  unique_queries = unique(c(queryHits(overlaps1), queryHits(overlaps2)))
  unique_queries.gr = all.exons.gr[unique_queries]
  print(paste('[NOTE] Number of samples I uneseccary lose, due to different versions of the exon =',
    length(which(countOverlaps(chosenExons.test.gr, unique_queries.gr) > 1)) ) )
  
  # remove samples which are not good
  filtered_Upintrons.df = Upintrons.df[-removeIntrons2,]
  filtered_Dointrons.df = Dointrons.df[-removeIntrons2,]
  chosenExons = chosenExons[-removeIntrons2,]
  
  print(paste0('[NOTE] number of upstream introns for testing = ',nrow(filtered_Upintrons.df)))
  print(paste0('[NOTE] number of exons for testing = ',nrow(chosenExons)))
  print(paste0('[NOTE] number of downstream introns for testing = ',nrow(filtered_Dointrons.df)))
  
  # checkpoint
  checkpoint2 = paste0('data/filtered/PSI.filtered.chosenExons.',cs,'.',shift,'.',inlength,'.',celltype,'.tsv')
  if ( file.exists(checkpoint2) ) {
    file.remove(checkpoint2)
  }
  write.table(chosenExons, file = checkpoint2, sep='\t', append=F)
  
  #######################################
  ## Methylation content and C content ##
  #######################################

  reference = BSgenome.Hsapiens.UCSC.hg19
  
  #############
  ### Exons ### 
  #############
  
  # create GRanges object out of the exons
  grExons <- GRanges(seqnames = chosenExons[,1], 
                     ranges = IRanges(start = as.numeric(chosenExons[,2]), 
                                      end = as.numeric(chosenExons[,3])), 
                     strand = chosenExons[,4])
  
  # count methylation loci per exon
  overlapcountsExons <- countOverlaps(grExons,mCpGs.gr)
  
  # create datatable
  Ematrix <- matrix(ncol=7, nrow=length(grExons))
  colnames(Ematrix) <- c("ID", "#Methylathions", "seqlength", 
                         'lengthNorm', 'CpG', 'mCpG/CpG', 'FPKM')
  
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
  
  Ematrix[,7] <- chosenExons[,ncol(chosenExons)]
  
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
  
  ##################
  ### DOWNSTREAM ### 
  ##################
  
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
  
  print('[NOTE] found introns')
  
  #########################
  ##  Methylation Plots  ##
  #########################
  
  # lists holding data
  first.df = c()
  second.df = c()
  
  # make plots for methyaltion content and mCpG/CpG ratio 
  methC = F
  for ( i in 1:2 ) {
    
    path_for_plots = "plots/"
    
    if(shortening) {
      path_for_plots = paste0(path_for_plots,'short/')
    } else {
      path_for_plots = paste0(path_for_plots,'unshort/')
    }
    
    print(paste('[NOTE] path of image = ', path_for_plots))
    
    if (cs == F) {
      sgn = '<='
      ge = 'exclusion'
    } else {
      sgn = '>='
      ge = 'inclusion'
    }
    
    title = 'Methylation Rate'
    values.up.l = as.numeric(UpImatrix[,4])
    values.ex.l = as.numeric(Ematrix[,4])
    values.do.l = as.numeric(DoImatrix[,4])
    boxplotlimits = c(0.0, 0.4)
    barplotlimits = c(0.0, 0.05)
    
    if ( methC == T ) {
      title = 'mCpG/CpG Ratio'
      values.up.l = as.numeric(UpImatrix[,6])
      values.ex.l = as.numeric(Ematrix[,6])
      values.do.l = as.numeric(DoImatrix[,6])
      boxplotlimits = c(0.0, 2.0)
      barplotlimits = c(0.0, 2.0)
    }
    
    # need datafram for ggplot
    data.df = data.frame(c(rep('Up', nrow(UpImatrix)), rep('Ex', nrow(Ematrix)), rep('Do', nrow(DoImatrix))), 
                         c(values.up.l, values.ex.l, values.do.l))
    
    colnames(data.df) = c('group', 'values')
    # stop alphabetic ordering of ggplot
    data.df$group <- factor(data.df$group, levels=unique(data.df$group))
    
    if ( methC == F ){
      first.df = data.df
    } else {
      second.df = data.df
    }
    
    # change to mCpG/CpG-ratio
    methC = T   
  }
    
  return( list(first.df, second.df, Ematrix[,3], Ematrix[,7]) )
  sink()
}

# evaluates the significances of the tri-region (up,ex,do)
# i.e. how often the exon is more methylated than the up/do-introns
signMeth = function(up,ex,do){
  
  k = 0
  for(i in 1:nrow(ex) ){
    if ( ex[i,4] > up[i,4] & ex[i,4] > do[i,4] ) {
      k = k + 1
    }
  }
  
  # total number of different methylations
  n = nrow(ex)
  
  # sample portion
  p_hat = k/n
  
  return(p_hat)
}

###############
## Get data  ##
###############

args = commandArgs(trailingOnly = TRUE)
celltype <- args[1]
if(is.na(celltype)){
  warning("[WARNING] define correct celltype: IMR90, Gm12878 and H1hesc")
}

plotcellnames = ''

if ( celltype == 'IMR90' ) { 
  plotcellnames <- 'IMR-90'
  # load methylation data 
  # WGBS for methylated CpGs
  mCpGs.df = fread('data/mCpG_IMR90.tsv', sep='\t', header=TRUE)
  
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
  CpGs.df <- fread('data/CpG_IMR90.tsv', sep='\t', header=TRUE)
  
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
  plotcellnames <- 'GM12878'
  # load methylation data 
  # WGBS for methylated CpGs
  mCpGs.df = fread('data/mCpG_Gm12878.tsv', sep='\t', header=TRUE)
  
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
  CpGs.df <- fread('data/CpG_Gm12878.tsv', sep='\t', header=TRUE)
  
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
  plotcellnames <- 'H1-hESC'
  # load methylation data 
  # WGBS for methylated CpGs
  mCpGs.df = fread('data/mCpG_H1hesc.tsv', sep='\t', header=TRUE)
  
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
  CpGs.df <- fread('data/CpG_H1hesc.tsv', sep='\t', header=TRUE)
  
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
  psi.df = read.delim("/data/PSI.10.H1hesc.tsv", header=T)
  row.names(psi.df) = NULL
}

# Load Expression data
expression.df = fread('data/new_counts.tsv', sep='\t', header=FALSE)

keep = which(!as.character(expression.df[[1]]) %in% c('chrM','chrX','chrY'))
expression.df = expression.df[keep,]

width <- abs(expression.df[[3]] - expression.df[[2]])
count <- expression.df[[8]]
N <- sum(as.numeric(count))

expression.df =  cbind(expression.df,
                        mapply(function(C,L) ((10e9*C)/(N*L)), 
                               as.numeric(count), 
                               as.numeric(width)))

#####################
## Get other data  ##
#####################

# get exons for hg19 
all.exons.df = fread('data/sorted_exons.bed', sep='\t', header=FALSE)

######################
## Creathe folders  ##
######################

dirs.l = c('plots/log',
           '/plots/whole')

# delete folders if they are exist
for ( i in 1:length(dirs.l) ) {
  if ( file.exists(dirs.l[i]) ) { 
    unlink(dirs.l[i], recursive = T)
    dir.create(dirs.l[i])
  } else {
    dir.create(dirs.l[i])
  }
}

###################
## Run Pipeline  ##
###################

f = F
s = T
zero = F

# function(cs, shortening, shift, inlength, zero)
first = pipe(f, F, 0.0, 500, zero, celltype, plotcellnames)
closeAllConnections()
second = pipe(s, F, 0.0, 500, zero, celltype, plotcellnames)
closeAllConnections()

# delete files in plots folder 
if ( length(list.files(path='plots/',pattern = "\\.eps$")) != 0 ) { 
  unlink('plots/*.eps')
}

# make plot about length distribution of the exons
length.df = c(rep(as.character(as.numeric(f)),length(first[[3]])), 
              rep(as.character(as.numeric(s)),length(second[[3]])))
length.df = cbind(data.frame(length.df), c(first[[3]], second[[3]]))
colnames(length.df) = c('group','values')

# t.test for the length 
length.t <- t.test(first[[3]], second[[3]])$p.value

postscript(file = paste0('plots/boxplot_length_compare','.eps'), family="Times")      
print(ggplot(length.df, aes(x=group, y=values))  + 
        theme(text=element_text(family="Times")) + 
        stat_boxplot(geom ='errorbar') + 
        geom_boxplot(fill=c('orange1','lightblue')) +
        guides(fill=FALSE) + 
        xlab("") +
        ylab('Length of Exon') +
        scale_y_continuous(trans='log10') +
        ggtitle(plotcellnames) +
        annotate("text", x = 1.5, y = max(c(first[[3]], second[[3]]))
                 , label = format(length.t), size = 5) +
        theme(text = element_text(size=20), axis.title.y=element_text(margin=margin(0,20,0,0)))
)
dev.off()

# make plot about expression rate 
mCpG_CpG_Ratio_table <- rbind(first[[2]], second[[2]])
mCpG_Cpg_Ratios <- mCpG_CpG_Ratio_table[which(mCpG_CpG_Ratio_table[,1] == 'Ex'),2]

expressionslist <- c(first[[4]], second[[4]])

hundredperzent_Exons <- expressionslist[which(mCpG_Cpg_Ratios >= 1.0)]
zeroperzent_Exons <- expressionslist[which(mCpG_Cpg_Ratios <= 0.0)]
inbetween_Exons <- expressionslist[intersect(which(mCpG_Cpg_Ratios > 0.0), which(mCpG_Cpg_Ratios < 1.0)) ]

# t.test for the length 
t1 <- t.test(inbetween_Exons, hundredperzent_Exons)$p.value
t2 <- t.test(inbetween_Exons, zeroperzent_Exons)$p.value
t3 <- t.test(hundredperzent_Exons, zeroperzent_Exons)$p.value

expr.df = c(rep('Average',length(inbetween_Exons)), 
            rep('Hyper',length(hundredperzent_Exons)),
            rep('Hypo',length(zeroperzent_Exons)))
expr.df = cbind(data.frame(expr.df), c(inbetween_Exons, hundredperzent_Exons, zeroperzent_Exons))
colnames(expr.df) = c('group','values')

postscript(file = paste0('plots/boxplot_epression_compare','.eps'),  family="Times")      
print(ggplot(expr.df, aes(x=group, y=values))  + 
        theme(text=element_text(family="Times")) + 
        stat_boxplot(geom ='errorbar') + 
        geom_boxplot(fill=c('orange1','lightblue','green')) +
        guides(fill=FALSE) + 
        xlab("") +
        ylab('FPKM') +
        scale_y_continuous(trans='log10') +
        ggtitle(plotcellnames) +
        annotate("text", x = c(1.5,2.5), 
                 y = rep(max(c(inbetween_Exons, hundredperzent_Exons, zeroperzent_Exons), na.rm=TRUE),2)
                 , label = c(format(t1), format(t3)), size = 5) +
        theme(text = element_text(size=20), axis.title.y=element_text(margin=margin(0,20,0,0)))
)
dev.off()

type = 'whole'
methC = F
for ( i in 2:2 ) {

  if ( i == 2 ){
    methC = T
  }

  f = F
  s = T
  
  if ( methC == F  ) {
    first.df = first[[1]]
    second.df = second[[1]]
  } else {
    first.df = first[[2]]
    second.df = second[[2]]
  }
  
  f.up = first.df[which(first.df[,1] == 'Up'),2]
  f.ex = first.df[which(first.df[,1] == 'Ex'),2]
  f.do = first.df[which(first.df[,1] == 'Do'),2]
  
  s.up = second.df[which(second.df[,1] == 'Up'),2]
  s.ex = second.df[which(second.df[,1] == 'Ex'),2]
  s.do = second.df[which(second.df[,1] == 'Do'),2]
  
  # if i = 3 then look only at exon which can not be methylated (no CpG sites)
  if ( i == 3 ) {
  
    methC = T
    
    type = 'only_notpossible'
      
    # look for NAs in exons for the mCpg/Cpg ratio
    f.nas.ex = which(f.ex == -1)
    s.nas.ex = which(s.ex == -1)
  
    # remove NAs from data
    f.up = f.up[f.nas.ex]
    f.ex = f.ex[f.nas.ex]
    f.do = f.do[f.nas.ex]
  
    s.up = s.up[s.nas.ex]
    s.ex = s.ex[s.nas.ex]
    s.do = s.do[s.nas.ex]
  }
  
  # if i = 4 then look at exon which can be methylated (without exon with no CpG sites)
  if ( i == 4 ) {
    
    methC = T
    
    type = 'without_notpossible'
    
    # look for NAs in exons for the mCpg/Cpg ratio
    f.nas.ex = which(f.ex == -1)
    s.nas.ex = which(s.ex == -1)
    
    # remove NAs from data
    f.up = f.up[-f.nas.ex]
    f.ex = f.ex[-f.nas.ex]
    f.do = f.do[-f.nas.ex]
    
    s.up = s.up[-s.nas.ex]
    s.ex = s.ex[-s.nas.ex]
    s.do = s.do[-s.nas.ex]
  }
  

  # set the introns with NAs to zero
  f.up[which(f.up == -1)] = 0
  f.do[which(f.do == -1)] = 0
  f.ex[which(f.ex == -1)] = 0

  s.up[which(s.up == -1)] = 0
  s.do[which(s.do == -1)] = 0
  s.ex[which(s.ex == -1)] = 0
  
  t_up = t.test(f.up, s.up)
  t_ex = t.test(f.ex, s.ex)
  t_do = t.test(f.do, s.do)
  
  print(t_up)
  print(t_ex)
  print(t_do)
  
  # change dataframes
  first.df = data.frame( c( rep('Up', length(f.up)), rep('Ex', length(f.ex)),
                            rep('Do', length(f.do)) ),
                         c(f.up, f.ex, f.do) )
  colnames(first.df) = c('group', 'values')

  second.df = data.frame( c( rep('Up', length(s.up)), rep('Ex', length(s.ex)),
                             rep('Do', length(s.do)) ),
                          c(s.up, s.ex, s.do) )
  colnames(second.df) = c('group', 'values')
  
  # create plot
  first.gg = summarySE(first.df, measurevar='values', groupvars='group')
  second.gg = summarySE(second.df, measurevar='values', groupvars='group')
  
  third.gg = first.gg
  third.gg[,2:5] = (log2(first.gg[,2:5]))/log2(second.gg[,2:5])
  
  title = 'Methylation Rate'
  barplotlimits = c(0.0, 2.0)
  textlimits = c( (2.0 - 0.1), (2.0 - 0.1), (2.0 - 0.1) )
  if ( methC == T ) {
    title = 'mCpG/CpG Ratio'
    barplotlimits = c(0.0, 2.5)
    textlimits = c( (2.5 - 0.1), (2.5 - 0.1), (2.5 - 0.1) )
    if( zero == T ){
      barplotlimits = c(0.0, 2.5)
      textlimits = c( (2.5 - 0.1), (2.5 - 0.1), (2.5 - 0.1) )
    }
  }
  
  postscript(file = paste0('plots/',type,'/barplot_',methC,'_foldChange','.eps'),  family="Times")  
  print(ggplot(third.gg, aes(x=group, y=values, fill=group)) + 
          theme(text=element_text(family="Times")) + 
          geom_bar(position=position_dodge(), stat="identity") +
          guides(fill=FALSE) + 
          xlab("") +
          ylab("Fold-Change") +
          ggtitle(paste('Fold-Change of the',title,'of Exclusions to Inclusions')) + 
          annotate("text", x = c(1,2,3), 
                 y = textlimits,
                 label = c(format(third.gg[1,3]),format(third.gg[2,3]),format(third.gg[3,3]))
          ) + 
          scale_y_continuous(limits = barplotlimits)
  )
  dev.off()
  
  if ( f ) {
    f.char = as.character(1)
  } else {
    f.char = as.character(0)
  }  
  
  if ( s ) {
    s.char = as.character(1)
  } else {
    s.char = as.character(0)
  }  
  
  
  first.df[,1] = paste0(first.df[,1],f.char)
  second.df[,1] = paste0(second.df[,1],s.char)
  
  third.df <- first.df[which(first.df[,1] == 'Up0'),]
  third.df <- rbind(third.df, second.df[which(second.df[,1] == 'Up1'),])
  third.df <- rbind(third.df, first.df[which(first.df[,1] == 'Ex0'),])
  third.df <- rbind(third.df, second.df[which(second.df[,1] == 'Ex1'),])
  third.df <- rbind(third.df, first.df[which(first.df[,1] == 'Do0'),])
  third.df <- rbind(third.df, second.df[which(second.df[,1] == 'Do1'),])
  
  colnames(third.df) = c('group', 'values')
  # stop alphabetic ordering of ggplot
  third.df$group <- factor(third.df$group, levels=unique(third.df$group))
  
  third.gg = summarySE(third.df, measurevar='values', groupvars='group')
 
  print("\n")
  print("summary")
  print(third.gg)
 
  title = 'Methylation Rate'
  barplotlimits = c(0.0, 0.05)
  if ( methC == T ) {
    title = 'mCpG/CpG Ratio'
    barplotlimits = c(0.0, 1.0)
  }
  
  postscript(file = paste0('plots/',type,'/barplot_',methC,'_compare','.eps'),  family="Times")
  print(
    ggplot(third.gg, aes(x=group, y=values, fill=c('red','red','blue','blue','yellow','yellow'))) + 
          theme(text=element_text(family="Times")) + 
          geom_bar(position=position_dodge(), stat="identity") +
          geom_errorbar(aes(ymin=values-ci, ymax=values+ci),
                        width=.2,                    
                        position=position_dodge(.9)) +
          guides(fill=FALSE) + 
          xlab("") +
          ylab(title) +
          annotate("text", x = c(1.5,3.5,5.5), 
                    y = c(max(barplotlimits) - 0.005, max(barplotlimits) - 0.005, max(barplotlimits) - 0.005), 
                    label = c(format(t_up[[3]]),format(t_ex[[3]]),format(t_do[[3]]))
                   ) + 
          scale_y_continuous(limits = barplotlimits) + 
          ggtitle(paste(title,'of Inclusions (1) to Exclusions (0)'))
  )
  dev.off()
  
  boxplotlimits = c(0.0, 0.3)
  if ( methC == T ) {
    boxplotlimits = c(0.0, 1.5)
  }
  
  # build up boxplot
  postscript(file = paste0('plots/',type,'/boxplot_',methC,'_compare','.eps'),  family="Times")      
  print(ggplot(third.df, aes(x=factor(group), y=values))  + 
          theme(text=element_text(family="Times")) + 
          stat_boxplot(geom ='errorbar') + 
          geom_boxplot(fill=c('#ff4040','#6ca6cd','#ff4040','#6ca6cd','#ff4040','#6ca6cd')) +
          guides(fill=FALSE) + 
          xlab("") +
          ylab(title) +
          annotate("text", x = c(1.5,3.5,5.5), size = 5,
                   y = c(max(boxplotlimits) - 0.005, max(boxplotlimits) - 0.005, max(boxplotlimits) - 0.005), 
                   label = c(format(t_up[[3]]),format(t_ex[[3]]),format(t_do[[3]]))
          ) +
          scale_y_continuous(limits = boxplotlimits) + 
          ggtitle(plotcellnames) +
          theme(text = element_text(size=20), axis.title.y=element_text(margin=margin(0,20,0,0)))
  )
  dev.off()
  
  
  ############################################
  ## Plot differences in the subpopulations ##
  ############################################
  
  if ( i != 3 ) {
    # write to log-file
    log.file = paste0('plots/',type,'/plotDifferences_',methC,'.txt')
    file.create(log.file)
    sink(file=log.file)
    
    # look at 'high' or 'low' methylated population
    population.l = c('low','mid','high')
    
    # concentrate on upstream intron, exon or downstream intron 
    sub.list = c('ex')
    
    # set barplot limits based on methylation rate or mCpG/CpG ratio 
    barplotlimits = c(0.0,0.2)
    
    # set boxplot limits based on methylation rate or mCpG/CpG ratio 
    boxplotlimits = c(0.0,0.3)
    
    look.f.up = f.up
    look.f.ex = f.ex
    look.f.do = f.do
    
    look.s.up = s.up
    look.s.ex = s.ex
    look.s.do = s.do
    
    if (  methC == T ) {
      barplotlimits = c(0.0, 1.2)
      boxplotlimits = c(0.0, 1.2)
      
      population.l = c('low','mid','high', 'Hundred Percent','Zero Percent')
        
      # take the extremes out
      keep.f = intersect(which(look.f.ex > 0.0), which(look.f.ex < 1.0))
      keep.s = intersect(which(look.s.ex > 0.0), which(look.s.ex < 1.0))
      
      look.f.up = look.f.up[keep.f]
      look.f.ex = look.f.ex[keep.f]
      look.f.do = look.f.do[keep.f]
      
      look.s.up = look.s.up[keep.s]
      look.s.ex = look.s.ex[keep.s]
      look.s.do = look.s.do[keep.s]
    }
    
    for ( j in 1:length(population.l) ){
      # create for each onject in sub.list a barplot 
      for ( i in 1:length(sub.list) ){
      
        # get onject of sub.list
        subpoupulation = sub.list[i]
        population = population.l[j]
        
        #
        # For low methylation take objects which have values less than the 25th quantile 
        #
        if ( population == 'low' ){      
          
          # list for exons 
          if ( subpoupulation == 'ex' ) {
            f.sub.l = which(look.f.ex <= quantile(look.f.ex)[2])
            s.sub.l = which(look.s.ex <= quantile(look.s.ex)[2])
            sub = 'Exons'
          }
        }
        
        # For high methylation take objects which have values bigger than the 75th quantile 
        if ( population == 'high' ){
          
          # list for exons 
          if ( subpoupulation == 'ex' ) {
            f.sub.l = which(look.f.ex >= quantile(look.f.ex)[4])
            s.sub.l = which(look.s.ex >= quantile(look.s.ex)[4])
            sub = 'Exons'
          }
        } 
        
        if ( population == 'mid' ){
          
          # list for exons 
          if ( subpoupulation == 'ex' ) {
            f.sub.l = intersect(which(look.f.ex < quantile(look.f.ex)[4]), 
                                which(look.f.ex > quantile(look.f.ex)[2]))
            s.sub.l = intersect(which(look.s.ex < quantile(look.s.ex)[4]), 
                                which(look.s.ex > quantile(look.s.ex)[2]))
            sub = 'Exons'
          }
        }
        
        if ( population == 'Zero Percent'  ) {
          # look which exons are 0% methylated 
          look.f.up = f.up
          look.f.ex = f.ex
          look.f.do = f.do
          
          look.s.up = s.up
          look.s.ex = s.ex
          look.s.do = s.do
          
          f.sub.l = which(look.f.ex == 0.0)
          s.sub.l = which(look.s.ex == 0.0)
        }
        
        if ( population == 'Hundred Percent' ) {
          # look which exons are 100% methylated 
          look.f.up = f.up
          look.f.ex = f.ex
          look.f.do = f.do
          
          look.s.up = s.up
          look.s.ex = s.ex
          look.s.do = s.do
          
          f.sub.l = which(look.f.ex == 1.0)
          s.sub.l = which(look.s.ex == 1.0)
        }
        
        
        # collect all values
        f.upi.l = look.f.up[f.sub.l]
        s.upi.l = look.s.up[s.sub.l]
        
        f.exo.l = look.f.ex[f.sub.l]
        s.exo.l = look.s.ex[s.sub.l]
          
        f.doi.l = look.f.do[f.sub.l]
        s.doi.l = look.s.do[s.sub.l]
           
        print(population)
        print(sub)
        print(length(f.upi.l))
        print(length(s.upi.l))
        print('FoldChanges: Exons / Introns')
        print(log2(mean(c(f.exo.l,s.exo.l)))/log2(mean(c(f.upi.l,s.upi.l))))
        print(log2(mean(c(f.exo.l,s.exo.l)))/log2(mean(c(f.doi.l,s.doi.l))))
        print('FoldChanges: ExcExon / IncExc')
        print(log2(mean(f.exo.l))/log2(mean(s.exo.l)))
        print('')
        
        # get t.test
        t_up = t.test(f.upi.l, s.upi.l)[[3]]
        t_do = t.test(f.doi.l, s.doi.l)[[3]]
    
        # test for zero means
        if ( mean(f.exo.l) == 0.0 && mean(s.exo.l) == 0.0 ){
          t_ex = 'zero means'
        }
        if ( mean(f.upi.l) == 0.0 && mean(s.upi.l) == 0.0 ){
          t_up = 'zero means'
        }
        if ( mean(f.doi.l) == 0.0 && mean(s.doi.l) == 0.0 ){
          t_do = 'zero means'
        }
        
        if( population == 'Hundred Percent' || population == 'Zero Percent' ){
          t_ex = 1.0
        } else {
          t_ex = t.test(f.exo.l, s.exo.l)[[3]]
        }
        
        # need datafram for ggplot
        data.df = data.frame(c(rep(paste0('Up',f.char), length(f.upi.l)), rep(paste0('Up',s.char), length(s.upi.l)), 
                               rep(paste0('Ex',f.char), length(f.exo.l)), rep(paste0('Ex',s.char), length(s.exo.l)),
                               rep(paste0('Do',f.char), length(f.doi.l)), rep(paste0('Do',s.char), length(s.doi.l))), 
                             c(f.upi.l, s.upi.l, f.exo.l, s.exo.l, f.doi.l, s.doi.l))
        
        colnames(data.df) = c('group', 'values')
        # stop alphabetic ordering of ggplot
        data.df$group <- factor(data.df$group, levels=unique(data.df$group))
        # create summary for the barplot
        data.gg = summarySE(data.df, measurevar='values', groupvars='group')
        
        postscript(file = paste0('plots/',type,'/barplot_',methC,'_',population,'_',subpoupulation,'_compare','.eps'),  family="Times")
        print(
          ggplot(data.gg, aes(x=group, y=values, fill=c('Q1','Q3','Q1','Q3','Q1','Q3'))) + 
                theme(text=element_text(family="Times")) + 
                geom_bar(position=position_dodge(), stat="identity") +
                geom_errorbar(aes(ymin=values-ci, ymax=values+ci),
                              width=.2,                    
                              position=position_dodge(.9)) +
                guides(fill=FALSE) + 
                xlab("") +
                ylab(title) +
                scale_y_continuous(limits = barplotlimits) + 
                annotate("text", x = c(1.5,3.5,5.5), size = 5,
                         y = c(barplotlimits[2] - 0.05, barplotlimits[2] - 0.05, barplotlimits[2] - 0.05), 
                         label = c(format(t_up),format(t_ex),format(t_do))
                ) + 
                ggtitle(plotcellnames) +
            theme(text = element_text(size=20), axis.title.y=element_text(margin=margin(0,20,0,0)))
              )
        dev.off()
        
        
        # build up boxplot
        postscript(file = paste0('plots/',type,'/boxplot_',methC,'_',population,'_',subpoupulation,'_compare','.eps'),  family="Times")      
        print(ggplot(data.df, aes(x=group, y=values))  + 
                theme(text=element_text(family="Times")) + 
                stat_boxplot(geom ='errorbar') + 
                geom_boxplot(fill=c('orange1','orange1','lightblue','lightblue','gold','gold')) +
                guides(fill=FALSE) + 
                xlab("") +
                ylab(title) +
                annotate("text", x = c(1.5,3.5,5.5), size = 5,
                         y = c(boxplotlimits[2] - 0.05, boxplotlimits[2] - 0.05, boxplotlimits[2] - 0.05), 
                         label = c(format(t_up),format(t_ex),format(t_do))
                ) +
                scale_y_continuous(limits = boxplotlimits) + 
                ggtitle(plotcellnames) 
        )
        dev.off()
      }
    }
    sink()
    closeAllConnections()
  }
}




