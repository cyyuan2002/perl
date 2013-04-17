#!/usr/bin/perl
use strict;
use Getopt::Long;

##This program is used to find the ortholog genes from ortholog file
##ortholog file
#coer|100025     conc|15914      2.396
#coer|100045     conc|77812      0.11
#coer|100164     conc|85752      0.182
#coer|100209     conc|17764      1.335
#coer|100280     conc|80186      0.563
#coer|10040      conc|166923     2.396
#coer|10052      conc|8143       0.917
#coer|100537     conc|78372      2.396

##input file
#AMD1    YML035C
#BUB2    YMR055C
#CCT2    YIL142W
#COX15   YER141W

my %opts;
GetOptions(\%opts,"i=s","o=s","s=s","help");
if((!defined $opts{i} || !defined $opts{o} || !defined $opts{s})){
    &Usage();
}

my $glfile=$opts{i};
my $orthfile=$opts{o};
my $species=$opts{s};

my %genelist;
my %orthogenes;

open(my $fh_glfile,$glfile) || die "Can't open gene list file: $glfile\n";
while(<$fh_glfile>){
    chomp();
    my @lines=split(/\t/,$_);
    $genelist{$lines[1]}=$lines[0];
}
close $fh_glfile;

open(my $fh_orthfile,$orthfile) || die "Can't open orthologous file: $orthfile\n";
while(<$fh_orthfile>){
    chomp();
    my @lines=split(/\t/,$_);
    my @inforA=split(/\|/,$lines[0]);
    my @inforB=split(/\|/,$lines[1]);
    if($inforA[0] eq $species){
        $orthogenes{$inforA[1]}=(exists($orthogenes{$inforA[1]})) ? "$orthogenes{$inforA[1]},$lines[1]" : $lines[1];
    }
    if($inforB[0] eq $species){
        $orthogenes{$inforB[1]}=(exists($orthogenes{$inforB[1]})) ? "$orthogenes{$inforB[1]},$lines[0]" : $lines[0];
    }
}
close $fh_orthfile;

foreach my $key (keys %genelist){
    if(exists($orthogenes{$key})){
        print "$genelist{$key}\|$key\t$orthogenes{$key}\n";
    }
}

exit(1);


sub Usage(){
  print << "    Usage";

	Usage:  $0 

	<options>
		-i     Gene list
		-o     Orthologous file which is built by orthomcl
		-s     Species of the gene list (same as the species in orthologous file)
                
    Usage

	exit(0);
};