#!/usr/bin/perl
use strict;

my ($parAnnoFile,$parFilein)=@ARGV;

open(my $fh_annofile,$parAnnoFile) || die "Can't open file: $parAnnoFile\n";
my %annogene;
my %genes;
while(<$fh_annofile>){
    chomp();
    my @lines=split(/\t/,$_);
    my @genes=split(/,/,$lines[7]);
    for(my $i=0;$i<@genes;$i++){
	$annogene{$genes[$i]}=$lines[8];
    }
    $genes{$lines[8]}=$lines[7];
}
close $fh_annofile;

open (my $fh_filein,$parFilein) || die "Can't open file: $parFilein\n";
while(<$fh_filein>){
    chomp();
    my @lines=split(/\t/,$_);
    my $geneid;
    if(exists($annogene{$lines[1]})){
	$geneid=$annogene{$lines[1]};
	print "$lines[0]\t$lines[2]\t$genes{$geneid}\t$geneid\n";
    }
    else{
	$geneid=$lines[1];
	print "$lines[0]\t$lines[2]\t$geneid\n";
    }

}
close $fh_filein;
exit(0);
