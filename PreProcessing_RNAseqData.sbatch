#!/bin/sh -e

#SBATCH --time=10:00:00
#SBATCH --mem=80gb

module load bedtools
module load samtools

# $1 = output folder
# $2 = first bam file
# $3 = seond bam file

Raw_MergedFile=$1'/raw_merged'

echo '[INIT] samtools merge'

samtools merge -f $Raw_MergedFile'.bam' $2 $3

echo "[INIT] bedtools bamtobed"

bedtools bamtobed -bed12 -i $Raw_MergedFile'.bam' > $Raw_MergedFile'.bed'

echo "[FINISH]"

echo '[INIT] get split reads from reads.bed'

awk '$10 != 1 {print $0}' $Raw_MergedFile'.bed' > $1'/split_reads.bed'

rm $Raw_MergedFile'.bed'
rm $Raw_MergedFile'.bam'

echo '[FINISH]'

echo '[INIT] sort reads'

sort-bed --max-mem 2G $1'/split_reads.bed' > $1'/sorted_split_reads.bed'
rm $1'/split_reads.bed'

echo '[FINISH]'

echo '[INIT] change read IDs'

# change IDs to be memory efficient and change space to tab delim else error from bedmap
awk '{$4 = "r"count++; print}' $1'/sorted_split_reads.bed' | sed 's/ /\t/g' > $1'/chIDs_split_reads.bed'
rm $1'/sorted_split_reads.bed'

echo '[FINISH]'

echo '[INIT] convert from bed12 to bed6 (from blocked bed to single block bed)'

bedtools bed12tobed6 -i $1'/chIDs_split_reads.bed' > $1'/chIDs_split_reads_bed6.bed' 

echo '[FINISH]'

echo '[INIT] bedops bedmap'
# count overlaping splitted reads with exons with at least 8 bp overlap
bedmap --ec --skip-unmapped --echo-map-id --echo --count --bp-ovr 8 $1'/sorted_exons.bed' $1'/chIDs_split_reads_bed6.bed' > $1'/WGintersect.bed'
 
echo '[FINISH]'

echo '[INIT] get id-list'

# get from the bed file the chromosome name, start and end position, strand and ids
# (transcript, gene, exon from ensmbl)
cut -f 1,2,3,6,10 $1'/WGintersect.bed' | tr -s '\|' '\t' | cut -f 2,3,4,5,6,7 | cut -d ';' -f 1,2,10 | tr -s '\"\;' '\t' | awk -v OFS='\t' '{print $1, $2, $3, $4, $6, $8, $10}' > $1'/new_id-list.tsv'

echo '[FINISH]'

echo '[INIT] get scores'

# get from the bed file the score
cut -f 10 $1'/WGintersect.bed' | cut -d"|" -f 2 | tr -s  '\"\;\|' '\t' > $1'/new_score-list.tsv'

echo '[FINISH]'

echo '[INIT] get IDs'

# get from the bed file the read IDs
cut -f 1 $1'/WGintersect.bed' | cut -d"|" -f 1 > $1'/readID-list.tsv'

echo '[FINISH]'

echo '[INIT] combine IDs, scores and read IDs'

# combine ids, scores and read IDs
paste $1'/new_id-list.tsv' $1'/new_score-list.tsv' $1'/readID-list.tsv' > $1'/new_exon-score-list.tsv'

rm $1'/readID-list.tsv'
rm $1'/new_id-list.tsv'
rm $1'/new_score-list.tsv'
rm $1'/WGintersect.bed'

echo '[FINISH]'

echo '[INIT] create results'

# get unique lines (this should automatically prevent that you maybe have different read counts
# for the same exon > you would have two or more lines with same exon in the file
# which you have to check further)
cat $1'/new_exon-score-list.tsv' | cut -f 1,2,3,4,5,6,7,8,9 | sort | uniq | sed -e 's/^[ \t]*//' > $1'/new_counts.tsv'

rm $1'/new_exon-score-list.tsv'

echo '[FINISH]'

