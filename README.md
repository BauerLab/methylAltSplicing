

PreProcessing_RNAseqData.sbatch:
The batch job was used to map the RNA sequencing data to the available exon annotations from GENCODE using bedops, samtools and bedtools.
- input:
	1) output folder
	2) first bam file (reads)
	3) second bam file (reads)
- output:
	1) list of exons with mapped reads


ProbabilityInclusion.R:
Was used to calculate the maximum a posterior probability for the inclusion of the investigated exons.
- input:
	1) data from PreProcessing_RNAseqData.sbatch
	2) cell type (IMR90, Gm12878 or H1hesc) 
- output:
	1) A table located in the directory at the end of the script containing the maximum a posterior probabilities for the inclusion of the exons.


PreProcessingWGBS.R:
The script is used to declare the given CpGs as mCpGs, when the methylation rate: x => t1 && x <= t2. The script can also be modified to use a smoothed distribution with the bsseq package. 
- input:
  	1) bed files (WGBS of the cell)
	2) cell type (IMR90, Gm12878 or H1hesc)
	3) lower threshold for a CpG to be defined as a mCpG (t1)
	4) upper threshold for a CpG to be defined as a mCpG (t2)
- output:
	1) A list of identified mCpGs for each cell type.

runCufflinks.sbatch:
The batch job was used to apply Cufflinks to the transcript data of the available reads.
- input:
	1) output folder
	2) bam file
	3) cell type (IMR90, Gm12878 or H1hesc)
- output:
	1) Cufflinks standard output


AnalysisExonsIntrons.R:
The script was used to look at some characteristics of the exons and introns (like the length comparison or the Cpg/mCpG ratio comparison) and to filter our some of the exons. 
- input:
	1) data from ProabilityInclusion.R
	2) data from PreProcessingWGBS.R.R
	3) data of runCufflinks.sbatch
	4) cell type (IMR90, Gm12878 or H1hesc)
- output:
	1) various plots


CreateFeatures.R:
The Script was used to define the feature matrix for the ANN and GBM. The cell type has to be changed in the beginning of the script.
- input:
	1) data from AnalysisExonsIntrons.R
- output:
	1) list of features for the cell type

WriteANNMatrix.R:
The scripts writes the input matrix for the deep learning algorithm. Cell type can be changed in the beginning of the script.
- input:
	1) data from CreateFeatures.R
- output:
	1) Matrix for the training of an ANN.

AnalyseProfilesPlots:
The script was used to plot the methylation profiles before the data was fed to the ANN.
- input:
	1) data from CreateFeatures.R
- output
	1) various plots

GBM_Metropolis-Hastings:
The was used to run a Metropolis Hastings algorithm to fine tune the parameter set of a GBM. The cell type can be changed in the beginning of the script. There is also a vector (in the beginning) which selects the features for the optimisation.
- input:
	1) data from CreateFeatures.R
  	2) TRUE or FALSE 
		- TRUE = The script uses a default parameter set for the fine tuning found by grid search.
		- FALSE = The script uses a random parameter set and tries to fine tune the set.
- output:
	1) A log file that contains the last line of MH algorithm (best model).
	2) (optional) If the script is run as a batch job, then the log file should contain all steps of the MH.

ANN_Metropolis-Hastings.R:
T script runs a Metroplis Hastings algorithm to fine tune the parameter setting for an ANN. The architecture stays constant but all other parameters changes. The cell type can be changed in the beginning of the script. The features selected for the optimisation can be changed in the section "Load Data".
- input:
	1) data from WriteANNMatrix.R
  	2) units (number of units in the hidden layer)
  	3) batch size (start value of the batch size)
  	4) rate (start value of the learning rate)
  	5) dropout (start value of the dropout) 
  	6) regularisation (start value of the l2 regularisation) 
  	7) momentum (start value of the momentum parameter)
- output:
	1) A log file that contains the last line of MH algorithm (best model).
	2) (optional) If the script is run as a batch job, then the log file should contain all steps of the MH.

ANN.R:
The script used to analyse the prediction of the best ANN.
- input:
	1) data from WriteANNMatrix.R
- output:
	1) various plots

GBM.R:
The script used to analyse the prediction of the best GBM.
- input:
	1) data from WriteANNMatrix.R
- output:
	1) various plots

