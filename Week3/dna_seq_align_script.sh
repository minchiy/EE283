#!/bin/bash
#SBATCH --array=1-12
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --job-name=dna_align

module load bwa/0.7.17
module load samtools/1.14

cd /pub/minchiy/ee283_work/Raw_data/DNA_Seq
mkdir -p aligned

ls *_1.fq.gz | sed 's/_1.fq.gz//' > prefixes.txt
ref="/pub/minchiy/ee283_work/genome_indices/dm6.fa"
prefix=`head -n $SLURM_ARRAY_TASK_ID prefixes.txt | tail -n 1`

bwa mem -t $SLURM_CPUS_PER_TASK -M $ref \
  ${prefix}_1.fq.gz ${prefix}_2.fq.gz | \
  samtools view -bS - > aligned/${prefix}.bam

samtools sort aligned/${prefix}.bam -o aligned/${prefix}.sort.bam
rm aligned/${prefix}.bam
