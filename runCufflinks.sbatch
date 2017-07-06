#!/bin/sh -e

#SBATCH --time=20:00:00
#SBATCH --mem=80gb
#SBATCH --ntasks-per-node=10

module load bedtools
module load samtools
module load cufflinks

# $1 = output folder
# $2 = bam file
# $3 = celltype and file number

echo '[INIT] sort reads'

samtools sort $2 -o $1'/sorted_'$3'.bam'

echo '[FINISH]'

echo '[INIT] cufflinks'

cufflinks $1'/sorted_'$3'.bam' -p 10

echo '[FINISH]'



