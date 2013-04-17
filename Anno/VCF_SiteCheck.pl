#!/usr/bin/perl -w

#===============================================================================
#
#         FILE: VCF_SiteCheck.pl
#
#        USAGE: 
#
#  DESCRIPTION: This script is used to filter out the snps that not
#               covered by all the genomes.
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

my ($file_list,$vcf_output)=@ARGV;
die "Usage:$0 <File_lst> <VCF_output>\n" if (@ARGV<2);

open(my $fh_filelist, $file_list) || die "Can't open file: $file_list\n";
#-----sample of file_list--------
#C.neo_A1.pass.vcf  A1.sam.sort.realign.bam.covr
#C.neo_A5.pass.vcf  A5ev.sam.sort.realign.bam.covr
#--------------------------------
my %vcffiles;

while(<$fh_filelist>){
    chomp();
    my @lines=split(/\t/,$_);
    $vcffiles{$lines[0]}=$lines[1];
}
close $fh_filelist;

my $ismissingfile=0;
foreach my $key(keys %vcffiles){
    if(!(-e $key)){
        print "Can't find file: $key\n";
        $ismissingfile=1;
    }
    if (!(-e $vcffiles{$key})){
        print "Can't find file: $vcffiles{$key}";
        $ismissingfile=1;
    }
}

die if($ismissingfile==1);

my @zipfiles;

foreach my $key(keys %vcffiles){
    `bgzip -c $key > $key.gz`;
    `tabix -p vcf $key.gz`;
    push(@zipfiles,"$key.gz");
}

my $zipfileN=join(" ",@zipfiles);

`vcf-merge $zipfileN > $vcf_output.tmp`;

my @vcfheads;
my %SNPs;
my @snpsites;
my $lastchrom="";
my @chroms;
my @snpinfors;

open(my $fh_vcfmerge,"$vcf_output.tmp") || die "Can't open $vcf_output.tmp files";
while(<$fh_vcfmerge>){
    chomp();
    if(/^#/){
        push (@vcfheads,$_);
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

$lastchrom="";
my $sitescount=0;
my %missedsites;
foreach my $key(keys %vcffiles){
    my $bam_covr_file=$vcffiles{$key};
    print "File:$bam_covr_file\n";
    open(my $fh_covr_file,"$bam_covr_file");
    my @snpsites;
    my $ispass=0;
    while(<$fh_covr_file>){
        my @lines=split(/\t/,$_);
        if($lastchrom ne $lines[0]){
            @snpsites=@{$SNPs{$lines[0]}->{'sites'}};
            $lastchrom=$lines[0];
            $sitescount=0;
            $ispass=0;
        }
        next if($ispass==1);
        if($lines[1]==$snpsites[$sitescount]){
            if($lines[3]<1){
                $missedsites{$lines[0]}->{$sitescount}=1;
            }
            $sitescount++;
        }
        elsif($lines[1]>$snpsites[$sitescount]){
            $missedsites{$lines[0]}->{$sitescount}=1;
            $sitescount++;
        }
        $ispass=1 if($sitescount==@snpsites);
    }
    close $fh_covr_file;
}

open(my $fh_vcfout,">$vcf_output");
foreach my $key(@vcfheads){
    print $fh_vcfout "$key\n";
}
foreach my $key(@chroms){
    my @snpsites=@{$SNPs{$key}->{'sites'}};
    my @infors=@{$SNPs{$key}->{'infors'}};
    for(my $i=0;$i<@snpsites;$i++){
        print $fh_vcfout "$infors[$i]\n" if(!exists($missedsites{$key}->{$i}));
    }
}
close $fh_vcfout;

exit(0);
