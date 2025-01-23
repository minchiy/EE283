#!/bin/bash

SourceDir="/data/class/ecoevo283/public/Bioinformatics_Course/RNAseq/RNAseq384plex_flowcell01"
DestDir="/pub/minchiy/ee283_work/Raw_data/RNA_Seq"
File="/data/class/ecoevo283/public/Bioinformatics_Course/RNAseq/RNAseq384_SampleCoding.txt"

mkdir -p "$DestDir"

if [[ ! -f "$File" ]]; then
    echo "Error: Sample coding file $File not found!"
    exit 1
fi

tail -n +2 "$File" | while IFS=$'\t' read -r SampleNumber Multiplexi5index Lane i7index PlateName PlateRow PlateColumn PlateWell RILcode TissueCode Replicate FullSampleName; do
    custom_filename="${SampleNumber}_${Multiplexi5index}_${Lane}_${i7index}_${PlateName}_${PlateRow}${PlateColumn}_${PlateWell}_${RILcode}_${TissueCode}_${Replicate}_${FullSampleName}"
    
    for project in "$SourceDir"/Project_*; do
        READ1=$(find "$project" -type f -name "${SampleNumber}_${i7index}_${Lane}_R1_*.fastq.gz")
        READ2=$(find "$project" -type f -name "${SampleNumber}_${i7index}_${Lane}_R2_*.fastq.gz")
        
        if [[ -n $READ1 && -n $READ2 ]]; then
            ln -s "$READ1" "${DestDir}/${custom_filename}_R1.fq.gz"
            ln -s "$READ2" "${DestDir}/${custom_filename}_R2.fq.gz"
            echo "Linked: $custom_filename"
            break
        fi
    done
    
    if [[ -z $READ1 || -z $READ2 ]]; then
        echo "Warning: FASTQ files not found for $custom_filename"
    fi
done
