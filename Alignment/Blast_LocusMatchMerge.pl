#!/usr/bin/perl
use strict;

my ($file1,$file2,$dirfile)=@ARGV;
open (my $fh_file1,$file1) || die "Can't open file: $file1\n";
my %geneids;
my %locus;
my %locusdir;
my %dirs;
my %gene;
while(<$fh_file1>){
    chomp();
    my @lines=split(/\t/,$_);
    if($lines[2]==2){
	if(!exists($geneids{$lines[5]})){
	    $geneids{$lines[5]}=$lines[0];
	    $dirs{$lines[5]}=$lines[4];
	}
	else{
	    $geneids{$lines[5]}.=",$lines[0]";
	    $dirs{$lines[5]}.=",$lines[4]";
	}
    }
    $gene{$lines[5]}=$lines[3];
}
close $fh_file1;

open(my $fh_dirfile,$dirfile) || die "Can't open file: $dirfile\n";
while(<$fh_dirfile>){
    chomp();
    my @lines=split(/\t/,$_);
    $locus{$lines[0]}=$lines[1];
    $locusdir{$lines[0]}=$lines[2];
}
close $fh_dirfile;

open(my $fh_file2,$file2) || die "Can't open file: $file2\n";
while(<$fh_file2>){
    chomp();
    my @lines=split(/\t/,$_);
    my @ids=split(/,/,$lines[0]);
    my @Dirs=split(/,/,$lines[1]);
    my $idNs;
    my $dirs;
    for(my $i=0;$i<@ids;$i++){
	my $id=$ids[$i];
	#if($id eq "Locus_5732_0"){
	#    sleep 1;
	#}
	my $dir=$locusdir{$id};
	if($Dirs[$i] eq "-"){
	    $dir=~tr/[\+\-]/[\-\+]/;
	}
	if($dirs eq ""){
	    $dirs=$dir;
	}
	else{
	    $dirs.=",$dir";
	}
	if($idNs eq ""){
	    $idNs=$locus{$id};
	}
	else{
	    $idNs.=",$locus{$id}";
	}
    }

    if($lines[3] ne ""){
	if(exists($geneids{$lines[3]})){
	    $geneids{$lines[3]}.=",$idNs";
	    $dirs{$lines[3]}.=",$dirs";
	}
	else{
	    $geneids{$lines[3]}=$idNs;
	    $dirs{$lines[3]}=$dirs;
	}
    }
    else{
	if(exists($geneids{$lines[2]})){
	    $geneids{$lines[2]}.=",$idNs";
	    $dirs{$lines[2]}.=",$dirs";
	}
	else{
	    $geneids{$lines[2]}=$idNs;
	    $dirs{$lines[2]}=$dirs;
	}
    }
}
close $fh_file2;

foreach my $key(keys %geneids){
    print "$geneids{$key}\t$dirs{$key}\t$key\n";
}

exit(0);
