#!/usr/bin/perl
#===============================================================================
#
#         FILE: GATK_SNPCall.pl
#
#        USAGE: This program is used for SNP calling by GATK
#
#  DESCRIPTION:
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Yuan Chen
#      COMPANY: MMRL, Duke University Medical Center
#      VERSION:
#      CREATED:
#     REVISION:
#===============================================================================
use strict;
use Getopt::Long;
use Benchmark;

my %opts;

GetOptions(\%opts,"a=s","b=s","r=s","o=s","n=s","u:s","d:i","q:i","m:i");
if((!defined($opts{a})) ||(!defined($opts{b}))||(!defined($opts{o})) ||(!defined($opts{r}))|| !defined($opts{n})){
    &Usage();
}

my $fastqA=$opts{a};
my $fastqB=$opts{b};
my $filein=$opts{o};
my $reffile=$opts{r};
my $RGID=$opts{n};
my $RGPU=$opts{u};
my $call_conf=$opts{q};
my $emit_conf=$opts{m};
my $dcov=$opts{d};

$RGPU||="CAGATC";
$call_conf||=30.0;
$emit_conf||=10.0;
$dcov||=100;

die "Can't find input file: $fastqA\n" if(!(-e $fastqA));
die "Can't find input file: $fastqB\n" if(!(-e $fastqB));
die "Can't find reference file: $reffile\n" if (!(-e $reffile));
my $time=localtime();
print "Job started on $time\n";
my $timestamp1=Benchmark->new;
my $command;
$command="bwa index $reffile";
print stderr "$command\n\n";
`$command`;
$command="bwa aln $reffile $fastqA > $fastqA.sai";
print stderr "$command\n\n";
`$command`;
$command="bwa aln $reffile $fastqB > $fastqB.sai";
print stderr "$command\n\n";
`$command`;
$command="bwa sampe $reffile $fastqA.sai $fastqB.sai $fastqA $fastqB > $filein";
print stderr "$command\n\n";
`$command`;
$command="samtools view -bS $filein > $filein.bam";
print stderr "$command\n\n";
`$command`;
$command="samtools sort $filein.bam $filein.sort";
print stderr "$command\n\n";
`$command`;
$command="samtools rmdup $filein.sort.bam $filein.sort.rmdup.bam";
print stderr "$command\n\n";
`$command`;
$command="java -jar ~/bin/picard/MarkDuplicates.jar INPUT=$filein.sort.rmdup.bam OUTPUT=$filein.sort.rmdup.picard.bam METRICS_FILE=picard_info.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT";
print stderr "$command\n\n";
`$command`;
$command="samtools calmd -Abr $filein.sort.rmdup.picard.bam $reffile > $filein.sort.rmdup.picard.baq.bam";
print stderr "$command\n\n";
`$command`;
$command="java -jar ~/bin/picard/AddOrReplaceReadGroups.jar INPUT=$filein.sort.rmdup.picard.baq.bam OUTPUT=$filein.sort.rmdup.picard.baq.group.bam RGID=$RGID RGLB=$RGID RGPL=illumina RGPU=$RGPU RGSM=$RGID VALIDATION_STRINGENCY=SILENT";
print stderr "$command\n\n";
`$command`;
$command="samtools index $filein.sort.rmdup.picard.baq.group.bam";
print stderr "$command\n\n";
`$command`;
$command="java -jar ~/bin/GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reffile -I $filein.sort.rmdup.picard.baq.group.bam -o GATK.mapping.intervals";
print stderr "$command\n\n";
`$command`;
$command="java -jar ~/bin/GATK/GenomeAnalysisTK.jar -T IndelRealigner -R $reffile -I $filein.sort.rmdup.picard.baq.group.bam -targetIntervals GATK.mapping.intervals -o $filein.sort.rmdup.picard.baq.group.realign.bam";
print stderr "$command\n\n";
`$command`;
$command="java -jar ~/bin/GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $reffile -I $filein.sort.rmdup.picard.baq.group.realign.bam -o $filein.sort.rmdup.picard.baq.group.realign.bam.snp.vcf -dcov $dcov -stand_call_conf $call_conf -stand_emit_conf $emit_conf";
print stderr "$command\n\n";
`$command`;
$command="java -jar ~/bin/GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $reffile -I $filein.sort.rmdup.picard.baq.group.realign.bam --genotype_likelihoods_model INDEL -o $filein.sort.rmdup.picard.baq.group.realign.bam.indel.vcf -dcov $dcov -stand_call_conf $call_conf -stand_emit_conf $emit_conf";
print stderr "$command\n\n";
`$command`;
print stderr "All Jobs finished.\n";
my $timestamp2=Benchmark->new;
my $timerun=timediff($timestamp2,$timestamp1);
print stderr "Total time spend $timerun\n";
exit(1);

sub Usage #help subprogram
{
    print << "    Usage";

	Usage: $0 [options] -a <in1.fq> -b <in2.fq> -o <output.sam> -r <refseq.fa> -n <SN>

		-u            Barcode of the sample

		-d            Maximum of coverage, default 100

                -q            Minimal qual score for SNP call, default 30

                -m            Minimal cutoff score for SNP call, default 10

    Usage

	exit(0);
};
