#!/bin/bash
#SBATCH --job-name=dna_seq_problem1
#SBATCH --partition=standard
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=6:00:00
set -euo pipefail

cd /pub/minchiy/EE283/Week4
exec > >(tee problem1.log) 2>&1

DATADIR="/pub/minchiy/ee283_work/Raw_data/DNA_Seq"
WORKDIR="/pub/minchiy/EE283/Week4"
REF="$WORKDIR/dmel-all-chromosome-r6.13.fasta"
TMPDIR="$WORKDIR/tmp"

mkdir -p "$TMPDIR"

rm -f $WORKDIR/*.tmp.* $WORKDIR/*.temp.*

module load bwa/0.7.17
module load samtools/1.15.1

declare -A samples=(
    ["A4"]="ADL06"
    ["A5"]="ADL09"
)

for strain in "${!samples[@]}"; do
    sample="${samples[$strain]}"
    echo "Processing ${strain} (${sample})..."
    
    if [ -f "${WORKDIR}/${strain}_chrX_region.bam" ]; then
        echo "Final BAM for ${strain} already exists, skipping..."
        continue
    fi
    
    for rep in {1..3}; do
        echo "Processing replicate ${rep}..."
        R1="${DATADIR}/${sample}_${rep}_1.fq.gz"
        R2="${DATADIR}/${sample}_${rep}_2.fq.gz"
        temp_bam="${TMPDIR}/${sample}_${rep}.temp.bam"
        sorted_bam="${TMPDIR}/${sample}_${rep}.sorted.bam"
        
        echo "Aligning replicate ${rep}..."
        bwa mem -t 16 "$REF" "$R1" "$R2" | \
            samtools view -b -h -o "$temp_bam" -
        
        echo "Sorting BAM file..."
        samtools sort -T "$TMPDIR/${sample}_${rep}" \
                     -o "$sorted_bam" "$temp_bam"
        
        region_bam="${TMPDIR}/${strain}_${rep}_chrX_region.bam"
        echo "Extracting region..."
        samtools view -b -h -q 30 -o "$region_bam" \
            "$sorted_bam" chrX:1880000-2000000
        
        # 清理臨時文件
        rm -f "$temp_bam" "$sorted_bam"
    done
    
    # 合併重複
    echo "Merging replicates for ${strain}..."
    samtools merge -f "${WORKDIR}/${strain}_chrX_region.bam" \
        "${TMPDIR}/${strain}"_*_chrX_region.bam
    samtools index "${WORKDIR}/${strain}_chrX_region.bam"
    
    # 清理中間文件
    rm -f "${TMPDIR}/${strain}"_*_chrX_region.bam
done

# 清理臨時目錄
rm -rf "$TMPDIR"

echo "Pipeline completed successfully"
echo "Final outputs:"
ls -lh ${WORKDIR}/*_chrX_region.bam*
