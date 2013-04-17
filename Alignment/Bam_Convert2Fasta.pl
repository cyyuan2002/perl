#!/usr/bin/perl
#===============================================================================
#
#         FILE:Bam_Convert2Fasta.pl
#
#        USAGE:Bam_Convert2Fasta.pl <Bam_File> <Sites ex:chr1:1000-2000> <Output_File>
#
#  DESCRIPTION:This program is used to get a map region from Bam file and convert it to fasta and qual file
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Yuan Chen
#      COMPANY: Division of Infectious Disease, DUMC
#      VERSION: 1.0
#      CREATED: 10/13/2012
#     REVISION:
#===============================================================================

use strict;
use File::Temp;

my ($parBamfile,$parSites,$parOutfile,$parQualScore)=@ARGV;
if(@ARGV<3){
    die "Usage:$0 <Bam_File> <Sites ex:chr1:1000-2000> <Output_File>";
}

my %refchroms;
my @tempfiles;

open (my $fh_BAM, "samtools view -h $parBamfile |");
while(<$fh_BAM>){
    chomp();
    if(/^\@SQ/){
        my ($chrom)=($_=~/SN\:(\S+)/);
        my ($len)=($_=~/LN\:(\d+)/);
        $refchroms{$chrom}=$len;
    }
    last if(/^\@RG/);
}
close $fh_BAM;

my @sites=split(" ",$parSites);
my @newsites;
foreach my $site(@sites){
    $site=~/(\S+)\:(\d+)-(\d+)/;
    my $chrom=$1;
    my $sitestart=$2;
    my $siteend=$3;
    die "Input Error: start_site must less than end_site\n" if($sitestart > $siteend);
    die "Error: Can't find \"$chrom\" in $parBamfile\n" if(!exists($refchroms{$chrom}));
    if($siteend>$refchroms{$chrom}){
        print stderr "Warning\: End_site ($siteend) is larger than $chrom length ($refchroms{$chrom})\n";
        print stderr "Using $refchroms{$chrom} instead of $siteend\n";
        $siteend=$refchroms{$chrom};
    }
    my $newsite="$chrom\:$sitestart\-$siteend";
    push(@newsites,$newsite);
}

my $editedsites=join(" ",@newsites);

$parQualScore ||=33;

my $Tempfile=File::Temp::tempnam(".","Temp");
my $TempBam=$Tempfile.".bam";

`samtools index $parBamfile` if(! -e ("$parBamfile.bai"));

`samtools view -bh $parBamfile $editedsites > $TempBam`;
`bam2fastq -o $Tempfile\#.fastq $TempBam`;

push(@tempfiles,$TempBam);

my $OutFasfile=$parOutfile.".fas";
my $OutQualfile=$OutFasfile.".qual";

open(my $fh_fasfile,">$OutFasfile");
open(my $fh_qualfile,">$OutQualfile");

my @fastqs=("$Tempfile\_1.fastq","$Tempfile\_2.fastq");
push (@tempfiles,@fastqs);

foreach my $fastq(@fastqs){
    open(my $fh_fastq,"$fastq") || die "Can't open temp fastq file: $fastq\n";
    my $i=0;
    while(<$fh_fastq>){
        $i++;
        if($i==1){
            $_=~s/@//;
            print $fh_fasfile ">$_";
            print $fh_qualfile ">$_";
        }
        elsif($i==2){
            print $fh_fasfile "$_";
        }
        elsif($i==4){
            chomp();
            my @quals=split("",$_);
            my @phradscore;
            foreach my $qual(@quals){
                my $q = ord ($qual) - $parQualScore;
                push(@phradscore,$q);
            }
            print $fh_qualfile join(" ",@phradscore),"\n";
            $i=0;
        }
    }
    close $fh_fastq;
}

close $fh_fasfile;
close $fh_qualfile;

foreach my $tempfile(@tempfiles){
    unlink $tempfile;
}

exit(1);
