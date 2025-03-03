#!/bin/bash
#SBATCH --job-name=problem2
#SBATCH --partition=standard
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=6:00:00


set -euo pipefail

cd /pub/minchiy/EE283/Week4
exec > >(tee problem2.log) 2>&1

salloc -A ecoevo283 --ntasks=2 srun --pty /bin/bash -i

module load samtools/1.15.1
module load bedtools2/2.30.0
module load blast/2.12.0

WORKDIR="/pub/minchiy/EE283/Week4"
BLAST_DB="$WORKDIR/problem2/dfam_db"

A4="$WORKDIR/A4_chrX_region.bam"
A5="$WORKDIR/A5_chrX_region.bam"

mkdir -p "$WORKDIR/problem2"
cd "$WORKDIR/problem2"

echo "Getting read IDs mapping to the interval..."
samtools view "$A4" chrX:1880000-2000000 | cut -f1 > A4.IDs.txt
samtools view "$A5" chrX:1880000-2000000 | cut -f1 > A5.IDs.txt

echo "Extracting reads mapping to different chromosomes..."
samtools view "$A4" | grep -f A4.IDs.txt | \
    awk '{if($3 != "X"){printf(">%s\n%s\n",$1,$10)}}' > A4.fa

samtools view "$A5" | grep -f A5.IDs.txt | \
    awk '{if($3 != "X"){printf(">%s\n%s\n",$1,$10)}}' > A5.fa

echo "Calculating statistics..."
echo "Number of reads in A4 FASTA:"
cat A4.fa | grep ">" | wc -l
echo "Number of reads in A5 FASTA:"
cat A5.fa | grep ">" | wc -l
echo "Total IDs in A4:"
wc -l A4.IDs.txt
echo "Total IDs in A5:"
wc -l A5.IDs.txt

echo "Setting up conda environment..."
source ~/miniforge3/etc/profile.d/conda.sh
source ~/miniforge3/etc/profile.d/mamba.sh

if ! conda env list | grep -q "spades_env"; then
    echo "Creating SPAdes environment..."
    mamba create -y -n spades_env spades
fi

conda activate spades_env

echo "Running SPAdes assembly..."
mkdir -p assembly/spades-default
spades.py -o assembly/spades-default/ -s A4.fa --isolate > A4.messages.txt

echo "Assembly results saved in assembly/spades-default/contigs.fasta"

conda deactivate

echo "Preparing BLAST database..."
mkdir -p "$BLAST_DB"
cd "$BLAST_DB"

if [ ! -f "dfam.fasta" ]; then
    echo "Downloading Dfam database..."
    wget https://www.dfam.org/releases/current/families/Dfam.fasta.gz
    gunzip Dfam.fasta.gz
fi

if [ ! -f "dfam.nhr" ]; then
    echo "Making BLAST database..."
    makeblastdb -in dfam.fasta -dbtype nucl -out dfam
fi

cd "$WORKDIR/problem2"

echo "Running BLAST analysis..."
blastn -query assembly/spades-default/contigs.fasta \
       -db "$BLAST_DB/dfam" \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
       -evalue 1e-10 \
       -num_threads 16 \
       -out blast_results.txt

echo "Top BLAST hits:"
sort -k12,12nr blast_results.txt | head -n 10

echo "Pipeline completed"
