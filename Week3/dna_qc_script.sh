#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --job-name=dna_qc
#SBATCH --output=dna_qc_%j.out

module load fastqc/0.11.9
module load java/1.8.0
module load trimmomatic/0.39

cd /pub/minchiy/ee283_work/Raw_data/DNA_Seq
mkdir -p fastqc_results trimmed

# Run on ADL06_1 pair
fastqc -o fastqc_results ADL06_1_1.fq.gz ADL06_1_2.fq.gz

java -jar /opt/apps/trimmomatic/0.39/trimmomatic-0.39.jar PE \
  ADL06_1_1.fq.gz ADL06_1_2.fq.gz \
  trimmed/ADL06_1_1_trimmed.fq.gz trimmed/ADL06_1_1_unpaired.fq.gz \
  trimmed/ADL06_1_2_trimmed.fq.gz trimmed/ADL06_1_2_unpaired.fq.gz \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

fastqc -o fastqc_results trimmed/ADL06_1_1_trimmed.fq.gz trimmed/ADL06_1_2_trimmed.fq.gz
