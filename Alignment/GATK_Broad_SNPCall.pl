#!/usr/bin/perl
#===============================================================================
#
#         FILE: GATK_SNPCall.pl
#
#        USAGE: GATK_SNPCall.pl -a <in1.fq> -b <in2.fq> -o <output.sam> -r <refseq.fa> -p <threads_num> -n <SN>
#
#  DESCRIPTION: This pipeline is according to Broad
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Yuan Chen
#      COMPANY: MMRL, Duke University Medical Center
#      VERSION: 1.0
#      CREATED: 07/11/2012
#     REVISION:
#===============================================================================
use strict;
use Getopt::Long;
use Benchmark;

my %opts;

GetOptions(\%opts,"a=s","b=s","r=s","o=s","n=s","p:n");
if((!defined($opts{a})) ||(!defined($opts{b}))||(!defined($opts{o})) ||(!defined($opts{r}))|| !defined($opts{n})){
    &Usage();
}

my $fastqA=$opts{a};
my $fastqB=$opts{b};
my $filein=$opts{o};
my $reffile=$opts{r};
my $RGID=$opts{n};
my $RGPU=$opts{u};
my $threads=$opts{p};


$RGPU||="CAGATC";
$threads||=1;

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
$command="bwa aln -t $threads $reffile $fastqA > $fastqA.sai";
print stderr "$command\n\n";
`$command`;
$command="bwa aln -t $threads $reffile $fastqB > $fastqB.sai";
print stderr "$command\n\n";
`$command`;
$command="bwa sampe $reffile $fastqA.sai $fastqB.sai $fastqA $fastqB > $filein";
print stderr "$command\n\n";
`$command`;
$command="samtools view -bS $filein > $filein.bam";
print stderr "$command\n\n";
`$command`;
$command="java -jar ~/bin/picard/AddOrReplaceReadGroups.jar INPUT=$filein.bam OUTPUT=$filein.sort.bam SORT_ORDER=coordinate RGID=$RGID RGLB=$RGID RGPL=illumina RGPU=$RGPU RGSM=$RGID VALIDATION_STRINGENCY=SILENT";
print stderr "$command\n\n";
`$command`;
$command="samtools index $filein.sort.bam";
print stderr "$command\n\n";
`$command`;
$command="java -jar ~/bin/GenomeAnalysisTK-1.4-14-g2e47336/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reffile -I $filein.sort.bam -o GATK.mapping.intervals";
print stderr "$command\n\n";
`$command`;
$command="java -jar ~/bin/GenomeAnalysisTK-1.4-14-g2e47336/GenomeAnalysisTK.jar -T IndelRealigner -R $reffile -I $filein.sort.bam -targetIntervals GATK.mapping.intervals -o $filein.sort.realign.bam";
print stderr "$command\n\n";
`$command`;
$command="samtools index $filein.sort.realign.bam";
print stderr "$command\n\n";
`$command`;
$command="java -jar ~/bin/GenomeAnalysisTK-1.4-14-g2e47336/GenomeAnalysisTK.jar -T UnifiedGenotyper -nt $threads -R $reffile -I $filein.sort.realign.bam -o $filein.sort.realign.bam.raw.vcf -dt NONE -glm BOTH -G Standard --output_mode EMIT_ALL_SITES -A AlleleBalance --computeSLOD";
print stderr "$command\n\n";
`$command`;
$command="java -jar ~/bin/GenomeAnalysisTK-1.4-14-g2e47336/GenomeAnalysisTK.jar -T VariantFiltration -R $reffile -V $filein.sort.realign.bam.raw.vcf -o $filein.sort.realign.bam.filtered.vcf --logging_level ERROR --filterExpression \"AB > 0.2 || ((DP - MQ0) < 5) || ((MQ0 / (1.0 * DP)) >= 0.5)\" --filterName LowConfidence";
print stderr "$command\n\n";
`$command`;
$command="grep \'^#\\|AC=[1-9]\' $filein.sort.realign.bam.filtered.vcf > $filein.sort.realign.bam.filtered.passed.vcf";
print stderr "$command\n\n";
`$command`;
$command="grep \'^#\' $filein.sort.realign.bam.filtered.passed.vcf > $filein.sort.realign.bam.filtered.final.vcf";
print stderr "$command\n\n";
`$command`;
$command="grep \'PASS\' $filein.sort.realign.bam.filtered.passed.vcf \| awk \'\$4==\"A\"\|\|\$4==\"C\"\|\|\$4==\"G\"\|\|\$4==\"T\"\' \| awk \'\$5==\"A\"\|\|\$5==\"C\"||\$5==\"G\"||\$5==\"T\"\'  \>\> $filein.sort.realign.bam.filtered.final.vcf";
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

	Usage: $0 [options] -a <in1.fq> -b <in2.fq> -o <output.sam> -r <refseq.fa> -n <SN> -p <threads_num>

    Usage

	exit(0);
};
