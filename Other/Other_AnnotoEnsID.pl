#!/usr/bin/perl
use strict;

my ($parEnsIDfile,$parAnnofile)=@ARGV;
my %geneName;
my %geneDes;
my %geneEnsPro;


open(my $ensfile,$parEnsIDfile) || die "Can't open file: $parEnsIDfile\n";
while(<$ensfile>){
	chomp();
	my @lines=split(/\t/,$_);
	$geneName{$lines[0]}=$lines[1];
	$geneDes{$lines[0]}=$lines[2];
	$geneEnsPro{$lines[0]}=$lines[3];
}
close $ensfile;

open(my $annofile,$parAnnofile)|| die "Can't open file: $parAnnofile\n";
while(<$annofile>){
	chomp();
	my @lines=split(/\t/,$_);
	my @ids=split(/,/,$lines[7]);
	my $genename;
	my $genedes;
	my $geneproid;
	my $geneid;
	for(my $i=0;$i<@ids;$i++){
		if($genename eq ""){
			if(exists($geneName{$ids[$i]})){
				$genename=$geneName{$ids[$i]};
				$genedes=$geneDes{$ids[$i]};
				$geneproid=$geneEnsPro{$ids[$i]};
				$geneid=$ids[$i];
				#print "$lines[8]\t$genename\t$genedes\t$geneproid\n";
			}
		}
		else{
			next if(!exists($geneName{$ids[$i]}));
			if($genename ne $geneName{$ids[$i]}){
				#print "Error:\t$geneid\t$ids[$i]\n";
				$genename=$genename."/".$geneName{$ids[$i]};
				$genedes=$genedes."/".$geneDes{$ids[$i]};
				$geneproid=$geneproid."/".$geneEnsPro{$ids[$i]};
			}
		}
	}
	if($genename ne ""){
		print "$lines[8]\t$genename\t$genedes\t$geneproid\n";
	}
}
close $annofile;
exit(0); 