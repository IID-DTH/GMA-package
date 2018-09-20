
#!/bin/bash

ref=$1
sample=$2
seq1=$3
seq2=$4

$q1=$5
$q2=$6

bwa.0.7 mem ${ref} ${seq1} ${seq2} -R '@RG\tID:Pool\tSM:${sample}\tPL:illumina\tLB:pool\tPU:${seq}'|samtools view -Sb - >${sample}.bam
samtools view -F 4 -q $q1 -h  ${sample}.bam | perl -ne 'print if(/^@/||(/AS:i:(\d+)\tXS:i:(\d+)/&&($1-$2)> $q2))'|grep -v "XA:Z\|SA:Z" |samtools view -Sb - >${sample}.uniq.bam
java -jar AddionalTools/picard.jar  SortSam  I=${sample}.uniq.bam O=${sample}.sorted.bam  SORT_ORDER=coordinate
java -jar AddionalTools/picard.jar MarkDuplicates I=${sample}.sorted.bam O=${sample}.soted.mkdup.bam REMOVE_DUPLICATES=T  M=${sample}.marked_dup.txt
