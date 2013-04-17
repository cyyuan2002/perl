#!/usr/bin/perl
##This program is used to filter the short sequences from Fasta sequences;

use strict;

my ($filein,$length,$fileout)=@ARGV;
if(@ARGV<2){
    print "Usage:program fastafile length outfile\n";
    exit;
}
if($fileout eq ""){
    $fileout="$filein.lft";
}

my @fileN=split(/\./,$filein);

if(!(-e $filein)){
    print stderr "Can't find file $filein\n";
    exit;
}
open (fileout,">$fileout");

my $filetype=$fileN[scalar(@fileN)-1];

if($filetype eq "gz"){
    use Compress::Zlib;
    my $fh=gzopen($filein,"rb");
    my $seq;
    my $seqN;
    while($fh->gzreadline(my $line)){
	chomp($line);
	if($line=~/^>\S+/){
	    if(length($seq)>=$length){
		print fileout "$seqN\n$seq\n";
	    }
	    $seqN=$line;
	    $seq="";
	}
	else{
	    $seq=$seq.$line;
	}
    }
    if(length($seq)>=$length){
	print fileout "$seqN\n$seq\n";
    }
    $fh->gzclose();
}
else{
    open(filein,$filein);
    my $seq;
    my $seqN;
    while(<filein>){
	chomp();
	if(/^>\S+/){
	    if(length($seq)>=$length){
		print fileout "$seqN\n$seq\n";
	    }
	    $seqN=$_;
	    $seq="";
	}
	else{
	    $seq=$seq.$_;
	}
    }
    if(length($seq)>=$length){
	print fileout "$seqN\n$seq\n";
    }
    close filein;
}

close fileout;
