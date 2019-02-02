#!/bin/bash
prefixes=(fastq-test.chip_reads fastq-test.control_reads)
for outprefix in ${prefixes[@]}; do
    echo "Doing $outprefix now."
    bwa index yeast-genome.fa
    #Uncomment the following two lines for single-end test run.
    #bwa aln yeast-genome.fa ${outprefix}.fastq > ${outprefix}.sai
    #bwa samse yeast-genome.fa ${outprefix}.sai $file > ${outprefix}.sam
    #Uncomment the following line for paired-end test run.
    bwa mem yeast-genome.fa ${outprefix}_R1.fastq ${outprefix}_R2.fastq > ${outprefix}.sam
    samtools view -b ${outprefix}.sam > ${outprefix}.bam
    samtools sort ${outprefix}.bam > ${outprefix}.sorted.bam
    samtools index ${outprefix}.sorted.bam
done
igv fastq-test.chip_reads.sorted.bam fastq-test.control_reads.sorted.bam
