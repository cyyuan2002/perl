#!/usr/bin/perl -w

#===============================================================================
#
#         FILE: VCF_CoverageCheck.pl
#
#        USAGE: The script is used to check the mapping coverage based on the
#               SNP information in the VCF file
#
#  DESCRIPTION:
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Yuan Chen
#      COMPANY: Division of Infectious Disease, DUMC
#      VERSION: 1.0
#      CREATED: 
#     REVISION:
#===============================================================================

use strict;

my ($vcf_file,$covr_file)=@ARGV;
die "Usage:$0 <VCF_File> <Coverage_File>\n" if (@ARGV<2);

my %SNPs;
my @snpsites;
my $lastchrom="";
my @chroms;
my @snpinfors;

open(my $fh_vcffile,"$vcf_file") || die "Can't open $vcf_file files";
while(<$fh_vcffile>){
    chomp();
    if(/^#/){
        next;
    }
    else{
        my @lines=split(/\t/,$_);
        if($lastchrom ne $lines[0]){
            if($lastchrom ne ""){
                my @tmpsites=@snpsites;
                $SNPs{$lastchrom}->{'sites'}=\@tmpsites;
                my @tmpinfors=@snpinfors;
                $SNPs{$lastchrom}->{'infors'}=\@tmpinfors;
            }
            $lastchrom=$lines[0];
            @snpsites=();
            @snpinfors=();
            push(@chroms,$lines[0]);
        }
        push (@snpsites,$lines[1]);
        push (@snpinfors,$_);
    }
}
my @sites=@snpsites;
$SNPs{$lastchrom}->{'sites'}=\@sites;
my @tmpinfors=@snpinfors;
$SNPs{$lastchrom}->{'infors'}=\@tmpinfors;

open(my $fh_covr_file,"$covr_file");
@snpsites=();
@snpinfors=();
my $ispass=0;
my $sitescount=0;
while(<$fh_covr_file>){
    my @lines=split(/\t/,$_);
    if($lastchrom ne $lines[0]){
        @snpsites=@{$SNPs{$lines[0]}->{'sites'}};
        @snpinfors=@{$SNPs{$lines[0]}->{'infors'}};
        $lastchrom=$lines[0];
        $sitescount=0;
        $ispass=0;
    }
    next if($ispass==1);
    if($lines[1]==$snpsites[$sitescount]){
        if($lines[3]<1){
            print "$snpinfors[$sitescount]\n";  
        }
        $sitescount++;
    }
    elsif($lines[1]>$snpsites[$sitescount]){
        print "$snpinfors[$sitescount]\n";
        $sitescount++;
    }
    $ispass=1 if($sitescount==@snpsites);
}
close $fh_covr_file;

exit(1);