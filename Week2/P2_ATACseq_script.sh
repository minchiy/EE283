#!/bin/bash

SourceDir="/data/class/ecoevo283/public/Bioinformatics_Course/ATACseq"
DestDir="/pub/minchiy/ee283_work/Raw_data/ATAC_Seq"
File="${SourceDir}/README.ATACseq.txt"


mkdir -p "$DestDir"


if [[ ! -f "$File" ]]; then
    echo "Error: File $File not found!"
    exit 1
fi


tail -n +2 "$File" | grep -v -E '^$|^ED|^WD' | while IFS=$'\t' read -r barcode genotype tissue bioRep; do
    
    if [[ -z "$barcode" || -z "$genotype" || -z "$tissue" || -z "$bioRep" ]]; then
        continue
    fi

    
    READ1=$(find "$SourceDir" -type f -iname "*_${barcode}_R1.fq.gz")
    READ2=$(find "$SourceDir" -type f -iname "*_${barcode}_R2.fq.gz")

    
    if [[ -n $READ1 ]]; then
        ln -s "$READ1" "${DestDir}/${genotype}_${tissue}_${bioRep}_R1.fq.gz"
    else
        echo "Warning: READ1 not found for barcode $barcode"
    fi
    if [[ -n $READ2 ]]; then
        ln -s "$READ2" "${DestDir}/${genotype}_${tissue}_${bioRep}_R2.fq.gz"
    else
        echo "Warning: READ2 not found for barcode $barcode"
    fi

   
    echo "Processed: Barcode=$barcode, Genotype=$genotype, Tissue=$tissue, Replicate=$bioRep"
done


