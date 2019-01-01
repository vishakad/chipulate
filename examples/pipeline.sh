#!/bin/bash
files=(fasta-basicExample2.tsv.chip_reads.fastq fasta-basicExample2.tsv.control_reads.fastq)
for file in ${files[@]}; do
    outprefix=${file%%.fastq}
    echo "Doing $outprefix now."
    #bedtools getfasta -s -name -fi ~/Data/Genomes/yeast-genome.fa -bed $file -fo ${outprefix}.fa
    bwa aln ~/Data/Genomes/yeast-genome.fa ${outprefix}.fastq > ${outprefix}.sai
    bwa samse ~/Data/Genomes/yeast-genome.fa ${outprefix}.sai $file > ${outprefix}.sam
    samtools view -b ${outprefix}.sam > ${outprefix}.bam
    samtools sort ${outprefix}.bam > ${outprefix}.sorted.bam
    samtools index ${outprefix}.sorted.bam
done
macs2 callpeak -t fasta-basicExample2.tsv.chip_reads.sorted.bam -c fasta-basicExample2.tsv.control_reads.sorted.bam --nomodel --extsize 200 -n blah -g 1.21e7
igv fasta-basicExample2.tsv.chip_reads.sorted.bam fasta-basicExample2.tsv.control_reads.sorted.bam
