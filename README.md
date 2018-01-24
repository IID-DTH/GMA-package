# GMA-package

Genetics Modal Assignments

1. Prepare the reference fasta and its index file

   The first step is to prepare the reference file. Here we use BWA to map the reads and GATK to call SNV variatns. The index file should be prepared with BWA and picard programme. 
   Moreover, the GMA package need samtools, bedtools, GATK, Picard and Annovar tools to be pre-installed.
   

+ BWA index building

```{sh}
   bwa index test_data/ref/NC_016845.fa
```

+ GATK index building

   The [GATK index building](https://gatkforums.broadinstitute.org/gatk/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference) is refer to the GATK forums.
   
```{sh}
java -jar picard.jar CreateSequenceDictionary R=test_data/ref/NC_016845.fa O=test_data/ref/NC_016845.dict 
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

4. CNV and SNV annotation


5. CNV and SNV comparation
