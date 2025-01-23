#!/bin/bash
#SBATCH --array=1-12
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --job-name=atac_align

module load bwa/0.7.17
module load samtools/1.14
module load java/1.8.0
module load trimmomatic/0.39

cd /pub/minchiy/ee283_work/Raw_data/ATAC_Seq
mkdir -p aligned_trimmed aligned_untrimmed trimmed

ls *_R1.fq.gz | sed 's/_R1.fq.gz//' > prefixes.txt
ref="/pub/minchiy/ee283_work/genome_indices/dm6.fa"
prefix=`head -n $SLURM_ARRAY_TASK_ID prefixes.txt | tail -n 1`

# Trimmomatic
java -jar /opt/apps/trimmomatic/0.39/trimmomatic-0.39.jar PE \
  ${prefix}_R1.fq.gz ${prefix}_R2.fq.gz \
  trimmed/${prefix}_R1_trimmed.fq.gz trimmed/${prefix}_R1_unpaired.fq.gz \
  trimmed/${prefix}_R2_trimmed.fq.gz trimmed/${prefix}_R2_unpaired.fq.gz \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Align trimmed reads
bwa mem -t $SLURM_CPUS_PER_TASK -M $ref \
  trimmed/${prefix}_R1_trimmed.fq.gz trimmed/${prefix}_R2_trimmed.fq.gz | \
  samtools view -bS - > aligned_trimmed/${prefix}.bam

# Align untrimmed reads
bwa mem -t $SLURM_CPUS_PER_TASK -M $ref \
  ${prefix}_R1.fq.gz ${prefix}_R2.fq.gz | \
  samtools view -bS - > aligned_untrimmed/${prefix}.bam

# Sort both BAMs
samtools sort aligned_trimmed/${prefix}.bam -o aligned_trimmed/${prefix}.sort.bam
samtools sort aligned_untrimmed/${prefix}.bam -o aligned_untrimmed/${prefix}.sort.bam

# Clean up
rm aligned_trimmed/${prefix}.bam aligned_untrimmed/${prefix}.bam
