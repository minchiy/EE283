#!/bin/bash
#SBATCH --array=1-12
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --job-name=rna_align

module load hisat2/2.2.1
module load samtools/1.14

cd /pub/minchiy/ee283_work/Raw_data/RNA_Seq
ls *_R1.fq.gz | sed 's/_R1.fq.gz//' > prefixes.txt

prefix=`head -n $SLURM_ARRAY_TASK_ID prefixes.txt | tail -n 1`
hisat2_index="/pub/minchiy/ee283_work/genome_indices/dm6_trans"

hisat2 -p $SLURM_CPUS_PER_TASK -x $hisat2_index \
-1 ${prefix}_R1.fq.gz -2 ${prefix}_R2.fq.gz | \
samtools view -bS - > aligned/${prefix}.bam

samtools sort aligned/${prefix}.bam -o aligned/${prefix}.sort.bam
samtools index aligned/${prefix}.sort.bam
rm aligned/${prefix}.bam
