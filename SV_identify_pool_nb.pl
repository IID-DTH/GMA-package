#!/usr/bin/perl
use warnings;
use Getopt::Long;
my %prior=();
my($bam,$genome,$type,$ploidy);
GetOptions(
	'bam=s'		=> \$bam,
	'genome=s'	=> \$genome,
	'type=s'	=> \$type,
	'ploidy=s'	=> \$ploidy,

); 




my $cov=$bam.".cov";
unless(-e $cov and -s $cov){
`genomeCoverageBed -ibam $bam.bam -g $genome |awk '\$1!="genome"' >$cov`;
}

$post=$bam.".post";
unless(-e $post and -s $post){
`Rscript nb_fit.r $cov $ploidy  > $post`;
}
open(F2,"$post");
while(<F2>){
#print;
	if($_ ne "" ){
		chomp;
		@items=split;
		$gt{$items[0]}{$items[1]}=$items[5];
		$gq{$items[0]}{$items[1]}=$items[6];
	}else{
		die "#Error in $bam: Error in estimate genotyping. \n";
	}
}
close F2;


if($type eq "bga"){
$dep=$bam.".bga";
#print "$dep\n";
unless(-e $dep and -s $dep){
	`genomeCoverageBed -ibam $bam.bam -g $genome -bga >$dep`;
}
#$bga=shift;
$last_geno=0;
$last_start=0;
$last_end=0;
$last_chr=0;
$last_genoQ=0;
$last_depth=0;
open(F1,$dep);
while(<F1>){
	chomp;
	($chr,$start,$end,$depth)=split;
	unless(exists $gt{$chr}{$depth}){
		die "#Error in $bam: Error in estimate genotyping. \n";
	}
	$genotype=$gt{$chr}{$depth};
	$genoQ=$gq{$chr}{$depth};
#	print "$chr\t$start\t$end\t$depth\t$cut_0{$chr}\t$cut_1{$chr}\t$genotype\t$last_geno\n";
	if($last_chr ne $chr ){
		if($last_chr ne 0){
			$last_length=$last_end-$last_start;
			$last_genoQA=sprintf("%.3f",($last_genoQ/$last_length));
			$last_depthA=sprintf("%.3f",($last_depth/$last_length));
			print "$last_chr\t$last_start\t$last_end\t$last_length\t$last_geno\t$last_depthA\t$last_genoQA\n";
		}
		$last_chr=$chr;
		$last_geno=$genotype;
		$last_start=$start;
		$last_end=$end;
		$last_genoQ=$genoQ * ($end-$start);
		$last_depth=$depth * ($end-$start);
	}elsif($genotype ne $last_geno){
		if($last_end ne 0){
			$last_length=$last_end-$last_start;
			$last_genoQA=sprintf("%.3f",($last_genoQ/$last_length));
			$last_depthA=sprintf("%.3f",($last_depth/$last_length));
			print "$last_chr\t$last_start\t$last_end\t$last_length\t$last_geno\t$last_depthA\t$last_genoQA\n";
		}
		$last_start =$start;
		$last_end =$end;
		$last_geno=$genotype;
		$last_genoQ=$genoQ * ($end-$start);
		$last_depth=$depth * ($end-$start);
#		print "$last_genoQ\t$last_depth\n";
	}elsif($genotype eq $last_geno){
		$last_end=$end;
		$last_genoQ += $genoQ * ($end-$start);
		$last_depth += $depth * ($end-$start);
#		print "$last_genoQ\t$last_depth\n";
	}	
}
$last_length=$last_end-$last_start;
$last_genoQA=sprintf("%.3f",($last_genoQ/$last_length));
$last_depthA=sprintf("%.3f",($last_depth/$last_length));
print "$last_chr\t$last_start\t$last_end\t$last_length\t$last_geno\t$last_depthA\t$last_genoQA\n";
close F1;
	
}else{
$dep=$bam.".dep";
#print "$dep\n";
unless(-e $dep and -s $dep){
	`genomeCoverageBed -ibam $bam.bam -g $genome -d >$dep`;
}
open(F3,$dep);
while(<F3>){
	chomp;
	($chr,$pos,$depth)=split;
	unless(exists $gt{$chr}{$depth}){
                die "#Error in $bam: Error in  estimate genotyping. \n";
        }
	$genotype=$gt{$chr}{$depth};
	$genoQ=$gq{$chr}{$depth};
	print "$chr\t$pos\t$depth\t$genotype\t$genoQ\n";
}
close F3;
}

