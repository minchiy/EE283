#!/bin/bash

SourceDir_DNAseq="/data/class/ecoevo283/public/Bioinformatics_Course/DNAseq"
DestDir_DNAseq="/pub/minchiy/ee283_work/Raw_data/DNA_Seq"

mkdir -p "$DestDir_DNAseq"

FILES_DNAseq="$SourceDir_DNAseq/*"
for f in $FILES_DNAseq
do
   ff=$(basename "$f")
   echo "Processing DNAseq $ff file..."
   ln -s "$SourceDir_DNAseq/$ff" "$DestDir_DNAseq/$ff"
done
