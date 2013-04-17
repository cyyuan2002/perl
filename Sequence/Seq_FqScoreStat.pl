#!/usr/bin/perl
use strict;

use Compress::Zlib;

my ($parFileFq)=shift;
my $parFileMode;


my @FileNA=split(/\./,$parFileFq);
my $Filetype=$FileNA[scalar(@FileNA)-1];

if($Filetype eq "gz"){
    $parFileMode=1;
}
else{
    $parFileMode=0;
}

if($parFileMode==1){
	my $fhfqfile=gzopen($parFileFq,"rb");
	my $linecount;
	while($fhfqfile->gzreadline(my $line)){
		$linecount++;
		if($linecount==4){
			my @scores=split("",$line);
			for(my $i=0;$i<@scores-1;$i++){
				my $score=ord($scores[$i])-64;
				print "$score\t";
			}
			$linecount=0;
			print "\n";
		}	
	}
	close $fhfqfile;
}
else{
	open(my $fhfqfile,$parFileFq);
	my $linecount;
	while(<$fhfqfile>){
		$linecount++;
		if($linecount==4){
			my @scores=split("",$_);
			for(my $i=0;$i<@scores-1;$i++){
				my $score=ord($scores[$i])-64;
				print "$score\t";
			}
			print "\n";
			$linecount=0;
		}
	}
	close $fhfqfile;
}

