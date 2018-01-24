#!/bin/bash
bam=$1;
#gt=$2;
ploidy=$2;
ref=$3;

#print "bam:$bam\nSV_gt:$SV_gt\nploidy";

#${bai}=$bam.".bai";
if [ ! -f "${bam}.bam.bai" ]; 
then 
	samtools index ${bam}.bam;
fi


if [ ! -f "${bam}.raw.table" ];
then
java -jar GenomeAnalysisTK.jar -T UnifiedGenotyper -I ${bam}.bam -R ${ref} -ploidy ${ploidy} -o ${bam}.raw.vcf -dcov 100000;
java -jar GenomeAnalysisTK.jar -T VariantsToTable -V ${bam}.raw.vcf -R ${ref} -o ${bam}.raw.table -F CHROM -F POS -F REF -F ALT -F QUAL -F AC -F AN -F DP
fi

if [ -f "${bam}.corr.table" ]; 
then 
rm ${bam}.corr.table;
fi

pl=`expr ${ploidy} - 1`

for i in $( seq 1 $pl );do
	
	if [ ! -f "${bam}.gt.$i.vcf" ];
	then
		awk -v p=$i '{if($5==p){print $1":"$2+1"-"$3}}' ${bam}.gt > ${bam}.gt.$i.intervals;
		java -jar GenomeAnalysisTK.jar -T UnifiedGenotyper -I ${bam}.bam -L ${bam}.gt.$i.intervals -R ${ref} -ploidy $i -o ${bam}.gt.$i.vcf -dcov 100000;
	fi	
		java -jar GenomeAnalysisTK.jar -T VariantsToTable -V ${bam}.gt.$i.vcf -R ${ref} -o ${bam}.gt.$i.table -F CHROM -F POS -F REF -F ALT -F QUAL -F AC -F AN -F DP ;
		cat ${bam}.gt.$i.table >> ${bam}.corr.table	
done	

awk -v p=${ploidy} '{if($5==p){print $1":"$2+1"-"$3}}' ${bam}.gt > ${bam}.gt.${ploidy}.intervals;
java -jar GenomeAnalysisTK.jar  -T SelectVariants  -R ${ref} -V ${bam}.raw.vcf -L ${bam}.gt.${ploidy}.intervals  -o ${bam}.gt.${ploidy}.vcf 
java -jar GenomeAnalysisTK.jar -T VariantsToTable -V ${bam}.gt.${ploidy}.vcf -R ${ref} -o ${bam}.gt.${ploidy}.table -F CHROM -F POS -F REF -F ALT -F QUAL -F AC -F AN -F DP ;
cat ${bam}.gt.${ploidy}.table >> ${bam}.corr.table
