#!/bin/bash
files=(chip.bed control.bed)
for file in ${files[@]}; do
    outprefix=${file%%.*}
    echo "Doing $outprefix now."
    bedtools getfasta -s -name -fi ~/Data/Genomes/yeast-genome.fa -bed $file -fo ${outprefix}.fa
    bwa aln ~/Data/Genomes/yeast-genome.fa ${outprefix}.fa > ${outprefix}.sai
    bwa samse ~/Data/Genomes/yeast-genome.fa ${outprefix}.sai ${outprefix}.fa > ${outprefix}.sam
    samtools view -b ${outprefix}.sam > ${outprefix}.bam
    samtools sort ${outprefix}.bam > ${outprefix}.sorted.bam
    samtools index ${outprefix}.sorted.bam
done
