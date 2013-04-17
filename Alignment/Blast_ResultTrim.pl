#!/usr/bin/perl
use strict;

my ($parBlastFile,$parHitsNum,$parOutFile,$parOutNohits)=@ARGV;
if(@ARGV<2){
    print "Usage:$0 Blast_file Hit_num [Output_file] [Output_nohits_file]\n";
    exit(1);
}

open(my $blastfile,"$parBlastFile") || die "Can't open file: $parBlastFile\n";
if($parOutFile eq ""){
    $parOutFile = $parBlastFile.".trim";
}
if($parOutNohits eq ""){
    $parOutNohits=$parBlastFile.".nohit";
}

open(my $outfile,">$parOutFile");
open(my $outnohits,">$parOutNohits");

my $blastinfo;
my $isrecord=0;
my $resultcount=0;
my $nohits=0;
my $queryName;
while(<$blastfile>){
    my $line=$_;
    if(/^\w?BLAST\w?/){
	if ($blastinfo ne ""){
	    print $outfile $blastinfo;
	}
	$blastinfo="";
	$blastinfo=$line;
	$isrecord=1;
	$nohits=0;
	$resultcount=0;
    }
    elsif(/Query=\s(\S+)/){
	$queryName=$1;
	$blastinfo=$blastinfo.$line;
    }
    elsif(/^Sequences producing significant alignments/){
	$blastinfo=$blastinfo.$line;
	my $linetemp=<$blastfile>;
	$blastinfo=$blastinfo.$linetemp;
	$blastinfo=$blastinfo."\n";
	for(my $i=0;$i<$parHitsNum;$i++){
	    $linetemp=<$blastfile>;
	    if($linetemp=~/\S+\s+\d+\s+\S+/){
		$blastinfo=$blastinfo.$linetemp;
	    }
	    else{
		last;
	    }
	}
	$isrecord=0;
	$blastinfo=$blastinfo."\n";
    }
    elsif(/^>/){
	$resultcount++;
	if($resultcount<=$parHitsNum){
	    $blastinfo=$blastinfo.$line;
	    $isrecord=1;
	}
	else{
	    $isrecord=0;
	}
    }
    elsif(/No hits found/){
	$nohits=1;
	print $outnohits "$queryName\n";
	$blastinfo="";
    }
    elsif(/Database/){
	if($nohits==0){
	    $isrecord=1;
	    $blastinfo=$blastinfo.$line;
	}
	else{
	    $isrecord=0;
	}
    }

    else{
	if($isrecord==1){
	    $blastinfo=$blastinfo.$line;
	}
    }
}

if ($blastinfo ne ""){
    print $outfile $blastinfo;
}
close $blastfile;
exit(0);
