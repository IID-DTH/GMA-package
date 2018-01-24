# GMA-package
Genetics Modal Assignments

1.Prepare the fasta and its index file

   BWA index building
   
2.GATK index building

```{sh}
samtools faidx test_data/ref/NC_016845.fa
```

2. Fastq mapping and bam processing

```{sh}
./bam_Processing.sh test_data/ref/NC_016845.fa test_data/bam_result/test test_data/raw_fastq/test.1.fq test_data/raw_fastq/test.2.fq
```

3. CNV and SNV detection

```{sh}
./SV_identify_pool_gp.pl -bam test_data/bam_result/test.sorted.mkdup -genome test_data/ref/NC_016845.genome -type dep -ploidy 10 >test_data/bam_result/test.sorted.mkdup.gt_gq
./SV_corrected_SNP_calling.sh test_data/bam_result/test.sorted.mkdup 10 test_data/ref/NC_016845.fa
```
