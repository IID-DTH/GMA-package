

# GPA-package

Genetics Polymorphisms Assignments

1. Prepare the reference fasta and its index file

   The first step is to prepare the reference file. Here we use BWA to map the reads and GATK to call SNV variatns. The index file should be prepared with BWA and Picard programme. 
   Moreover, the GPA package need SAMtools, bedtools, GATK, Picard and Annovar tools to be pre-installed. These tools could be add to  AddionalTools dictionary. 

Download URLs

[SAMtools](https://github.com/samtools/samtools)

[bedtools](https://github.com/arq5x/bedtools2/blob/master/docs/index.rst)

[GATK](https://software.broadinstitute.org/gatk/download/)

[Picard](http://broadinstitute.github.io/picard/)

[ANNOVAR](http://annovar.openbioinformatics.org/en/latest/user-guide/download/)


+ BWA index building

```{sh}
   bwa index NC_016845.fa
```

+ GATK index building

   The GATK index building is refer to the GATK forums (https://gatkforums.broadinstitute.org/gatk/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference). 
 
```{sh}
java -jar YOUR_PICARD_PATHWAY/picard.jar CreateSequenceDictionary R=NC_016845.fa O=NC_016845.dict 
samtools faidx NC_016845.fa
```

+ Genome file for bedtools

   The genome file is as input for bedtools to calculate the depth in every position.
```{sh}
cut -f 1-2 NC_016845.fa.fai >NC_016845.genome
```

+ Annotation file for Annovar

   The SNV annotation is based on Annovar programme. The Annovar programe provide tools for personal gene definition databases ( http://annovar.openbioinformatics.org/en/latest/user-guide/gene/#create-your-own-gene-definition-databases-for-non-human-species).

```{sh}
mkdir YOUR_ANNOVAR_PATHWAY/annovar/kpndb/
cd YOUR_ANNOVAR_PATHWAY/annovar/kpndb/
../retrieve_seq_from_fasta.pl -format refGene -seqdir ./ NC_016845_refGene.txt
mv NC_016845_refGene.txt.fa NC_016845_refGeneMrna.facd ../../../
```

2. Fastq mapping and bam processing

   The input fastq was mapped to the reference genome using BWA mem programme with default parameters. We filtered the unmapped reads, low quality reads and multiple mapping reads. The low quality reads were defined as reads that mapping quality lower than 30. The multiple mapping reads were defined as reads that the difference of Phred score between best alignment and the secondary alignment more than 10. Then we use Picard to sort bam and remove duplicate reads. This process could be modified in your particularly application.

```{sh}
bwa.0.7 mem NC_016845.fa test.1.fq.gz test.2.fq.gz -R '@RG\tID:Pool\tSM:test\tPL:illumina\tLB:pool\tPU:test'|samtools view -Sb - >bam_result/test.bam
samtools view -F 4 -q 30 -h  bam_result/test.bam | perl -ne 'print if(/^@/||(/AS:i:(\d+)\tXS:i:(\d+)/&&($1-$2)>10))'|grep -v "XA:Z\|SA:Z" |samtools view -Sb - >bam_result/test.uniq.bam
java -jar AddionalTools/picard.jar  SortSam  I=bam_result/test.uniq.bam O=bam_result/test.sorted.bam  SORT_ORDER=coordinate
java -jar AddionalTools/picard.jar MarkDuplicates I=bam_result/test.sorted.bam O=bam_result/test.soted.mkdup.bam REMOVE_DUPLICATES=T  M=bam_result/test.marked_dup.txt
```

For simply using our default parameters, you could use the bam_Processing.sh programme.

```{sh}
./bam_Processing.sh NC_016845.fa bam_result/test test.1.fq.gz test.2.fq.gz
```

3. CNV and SNV detection

   The CNV detection is based on the GMA bayesian model. There are two mode in the package. Th "dep" mode output is the ploidy result in every position. The result maybe used in SNV compare mode.

```{sh}
./SV_identify_pool_gp.pl -bam bam_result/test.sorted.mkdup -genome NC_016845.genome -type dep -ploidy 3 >bam_result/test.sorted.mkdup.gt_gq
```

The "bga" mode result report ploidy result in BedGraph format and as the input for the SNV detection.

```{sh}
./SV_identify_pool_gp.pl -bam bam_result/test.sorted.mkdup -genome NC_016845.genome -type bga -ploidy 3 >bam_result/test.sorted.mkdup.gt
```

The SNV detection is based on GATK UnifiedGenotper as SNV caller, while taking the BedGraph ploidy result into consideration.  

```{sh}
./SV_corrected_SNP_calling.sh bam_result/test.sorted.mkdup.bam bam_result/test test_data/bam_result/test.sorted.mkdup.gt 3 NC_016845.fa
```

4. SNV annotation

The SNV annotation use ANNOVAR programme with personal made annotataion files. 

```{sh}
 awk '{print $1"\t"$2"\t"$2"\t"$3"\t"$4}'  test.corr.table |sed '1d' >test.corr.avinput
./AddionalTools/annovar/annotate_variation.pl test.corr.avinput ./YOUR_ANNOVAR_PATHWAY/annovar/kpndb/ -buildver NC_016845
```

5. CNV and SNV comparation

We compare the CNV/SNV files from two population with Fisher Exact test with R programme. 

```{sh}
./compare_2_table.pl -f1  test1.table -f2  test2.table -p1 10 -p2 10 >test1_2.compare 
Rscript fisher_result.r test_data/geno_result/test1_2.compare
```

The compare could also take the CNV dep result as input to calculate the undetect positions.

```{sh}
./compare_2_table.pl -f1  test1.table -f2  test2.table -p1 10 -p2 10 -g1 test1.gt_gq  -g2 test2.gt_gq >test1_2.compare 
Rscript fisher_result.r test_data/geno_result/test1_2.compare
```

