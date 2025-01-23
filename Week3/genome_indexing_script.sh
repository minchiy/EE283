#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --job-name=genome_index

module load bwa/0.7.17
module load samtools/1.14
module load picard-tools/2.27.1
module load hisat2/2.2.1

cd /pub/minchiy/ee283_work
mkdir -p genome_indices

ref="genome_indices/dm6.fa"
gtf="genome_indices/dm6.ncbiRefSeq.gtf.gz"

bwa index $ref
samtools faidx $ref
java -jar /opt/apps/picard-tools/2.27.1/picard.jar CreateSequenceDictionary R=$ref O=genome_indices/dm6.dict

python hisat2_extract_splice_sites.py $gtf > genome_indices/dm6.ss
python hisat2_extract_exons.py $gtf > genome_indices/dm6.exon
hisat2-build -p 8 --exon genome_indices/dm6.exon --ss genome_indices/dm6.ss $ref genome_indices/dm6_trans
