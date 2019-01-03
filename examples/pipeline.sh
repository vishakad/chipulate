#!/bin/bash
files=(fastq-test.chip_reads.fastq fastq-test.control_reads.fastq)
for file in ${files[@]}; do
    outprefix=${file%%.fastq}
    echo "Doing $outprefix now."
    bwa index yeast-genome.fa
    bwa aln yeast-genome.fa ${outprefix}.fastq > ${outprefix}.sai
    bwa samse yeast-genome.fa ${outprefix}.sai $file > ${outprefix}.sam
    samtools view -b ${outprefix}.sam > ${outprefix}.bam
    samtools sort ${outprefix}.bam > ${outprefix}.sorted.bam
    samtools index ${outprefix}.sorted.bam
done
igv fastq-test.chip_reads.sorted.bam fastq-test.control_reads.sorted.bam
