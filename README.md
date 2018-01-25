# GMA-package

Genetics Modal Assignments

1. Prepare the reference fasta and its index file

   The first step is to prepare the reference file. Here we use BWA to map the reads and GATK to call SNV variatns. The index file should be prepared with BWA and picard programme. 
   Moreover, the GMA package need samtools, bedtools, GATK, Picard and Annovar tools to be pre-installed. These tools were include in the AddionalTools dictionary. 
   

+ BWA index building

```{sh}
   bwa index test_data/ref/NC_016845.fa
```

+ GATK index building

   The GATK index building is refer to the GATK forums (https://gatkforums.broadinstitute.org/gatk/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference). 
 
```{sh}
java -jar AddionalTools/picard.jar CreateSequenceDictionary R=test_data/ref/NC_016845.fa O=test_data/ref/NC_016845.dict 
samtools faidx test_data/ref/NC_016845.fa
```

+ Genome file for bedtools

   The genome file is as input for bedtools to calculate the depth in every position.

```{sh}
cut -f 1-2 test_data/ref/NC_016845.fa.fai >test_data/ref/NC_016845.genome
```

+ Annotation file for Annovar

   The SNV annotation is based on Annovar programme. The Annovar programe provide tools for personal gene definition databases ( http://annovar.openbioinformatics.org/en/latest/user-guide/gene/#create-your-own-gene-definition-databases-for-non-human-species).

```{sh}
mkdir AddionalTools/annovar/kpndb/
cd AddionalTools/annovar/kpndb/
../retrieve_seq_from_fasta.pl -format refGene -seqdir ./ NC_016845_refGene.txt
mv NC_016845_refGene.txt.fa NC_016845_refGeneMrna.fa
cd ../../../
```

2. Fastq mapping and bam processing

   The input fastq was mapped to the reference genome and use Picard to sort bam and remove duplicate reads.

```{sh}
./bam_Processing.sh test_data/ref/NC_016845.fa test_data/bam_result/test test_data/raw_fastq/test.1.fq test_data/raw_fastq/test.2.fq
```

3. CNV and SNV detection

   The CNV detection is based on the GMA bayesian model. There are two mode in the package. Th "dep" mode output is the ploidy result in every position. 

```{sh}
./SV_identify_pool_gp.pl -bam test_data/bam_result/test.sorted.mkdup -genome test_data/ref/NC_016845.genome -type bga -ploidy 10 >test_data/bam_result/test.sorted.mkdup.gt
```
The "bga" mode result report ploidy result in BedGraph format and as the input for the SNV detection.

```{sh}
./SV_identify_pool_gp.pl -bam test_data/bam_result/test.sorted.mkdup -genome test_data/ref/NC_016845.genome -type dep -ploidy 10 >test_data/bam_result/test.sorted.mkdup.gt_gq
```

The SNV detection is 

```{sh}
./SV_corrected_SNP_calling.sh test_data/bam_result/test.sorted.mkdup 10 test_data/ref/NC_016845.fa
```



4. SNV annotation



```{sh}
./AddionalTools/annovar/annotate_variation.pl test.avinput ../AddionalTools/annovar/kpndb/ -buildver NC_016845
```

5. CNV and SNV comparation

We compare the CNV/SNV files from two population with Fisher Exact test with R programme. 
The compare could also take the CNV dep result as input to calculate the undetect positions. 

```{sh}
./compare_2_table.pl -f1 <(grep -v "," test1.table) -f2 <(grep -v "," test2.table) -p1 10 -p2 10 -g1 test1.gt_gq  -g2 test2.gt_gq >test1_2.compare 
Rscript fisher_result.r test1_2.compare
```
