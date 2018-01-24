#!/bin/bash

ref=$1
sample=$2
seq1=$3
seq2=$4

bwa.0.7 mem ${ref} ${seq1} ${seq2} -R '@RG\tID:Pool\tSM:${sample}\tPL:illumina\tLB:pool\tPU:${seq}'|samtools view -Sb - >${sample}.bam
java -jar AddionalTools/picard.jar  SortSam  I=${sample}.bam O=${sample}.sorted.bam  SORT_ORDER=coordinate
java -jar AddionalTools/picard.jar MarkDuplicates I=${sample}.sorted.bam O=${sample}.soted.mkdup.bam REMOVE_DUPLICATES=T  M=${sample}.marked_dup.txt
