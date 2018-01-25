#!/bin/bash
bam=$1;
sample=$2
gt=$3;
ploidy=$4;
ref=$5;

echo "./SV_corrected_SNP_calling.sh <bam file> <output_sample> <gt file> <ploidy> <reference fasta>\n";


#print "bam:$bam\nSV_gt:$SV_gt\nploidy";

#${bai}=$bam.".bai";
if [ ! -f "${bam}.bai" ]; 
then 
	samtools index ${bam};
fi


if [ ! -f "${sample}.raw.table" ];
then
java -jar AddionalTools/GenomeAnalysisTK.jar -T UnifiedGenotyper -I ${bam} -R ${ref} -ploidy ${ploidy} -o ${sample}.raw.vcf -dcov 100000;
java -jar AddionalTools/GenomeAnalysisTK.jar -T VariantsToTable -V ${sample}.raw.vcf -R ${ref} -o ${sample}.raw.table -F CHROM -F POS -F REF -F ALT -F QUAL -F AC -F AN -F DP
fi

if [ -f "${sample}.corr.table" ]; 
then 
rm ${sample}.corr.table;
fi

pl=`expr ${ploidy} - 1`

for i in $( seq 1 $pl );do
	
	if [ ! -f "${sample}.gt.$i.vcf" ];
	then
		awk -v p=$i '{if($5==p){print $1":"$2+1"-"$3}}' ${gt} > ${sample}.gt.$i.intervals;
		java -jar AddionalTools/GenomeAnalysisTK.jar -T UnifiedGenotyper -I ${bam} -L ${sample}.gt.$i.intervals -R ${ref} -ploidy $i -o ${sample}.gt.$i.vcf -dcov 100000;
	fi	
		java -jar AddionalTools/GenomeAnalysisTK.jar -T VariantsToTable -V ${sample}.gt.$i.vcf -R ${ref} -o ${sample}.gt.$i.table -F CHROM -F POS -F REF -F ALT -F QUAL -F AC -F AN -F DP ;
		cat ${sample}.gt.$i.table >> ${sample}.corr.table	
done	

awk -v p=${ploidy} '{if($5==p){print $1":"$2+1"-"$3}}' ${gt} > ${sample}.gt.${ploidy}.intervals;
java -jar AddionalTools/GenomeAnalysisTK.jar  -T SelectVariants  -R ${ref} -V ${sample}.raw.vcf -L ${sample}.gt.${ploidy}.intervals  -o ${sample}.gt.${ploidy}.vcf 
java -jar AddionalTools/GenomeAnalysisTK.jar -T VariantsToTable -V ${sample}.gt.${ploidy}.vcf -R ${ref} -o ${sample}.gt.${ploidy}.table -F CHROM -F POS -F REF -F ALT -F QUAL -F AC -F AN -F DP ;

cat ${sample}.gt.${ploidy}.table >> ${sample}.corr.table
sed '/CHROM/d' ${sample}.corr.table |sort -k 1,1 -k 2,2n|sed '1iCHROM\tPOS\tREF\tALT\tQUAL\tAC\tAN\tDP' >tmp
mv tmp ${sample}.corr.table

rm ${sample}.gt.*
