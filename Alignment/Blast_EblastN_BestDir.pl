#!/usr/bin/perl
use strict;

my ($filein,$fileout,$mode)=@ARGV;
if(@ARGV<2){
    print stderr "Usage:$0 input_file output_file\n";
    exit(0);
}

my ($fhin,$fhout);
open($fhin,$filein) || die "Can't open file $filein\n";
open($fhout,">$fileout");

my $QueryN;
my $QueryS;

while(<$fhin>){
    chomp();
    my @lines=split(/\t/,$_);
    if($QueryN ne $lines[0]){
	my $querys=$lines[2];
	my $querye=$lines[3];
	my $subs=$lines[4];
	my $sube=$lines[5];
	my ($dirq,$dirs);
	if($querys<$querye){
	    $dirq=1;
	}
	else{
	    $dirq=-1;
	}
	if($subs<$sube){
	    $dirs=1;
	}
	else{
	    $dirs=-1;
	}
	my $dir=$dirq*$dirs;
	if($dir>0){
	    print $fhout "$lines[0]\t$lines[11]\t+\n";
	}
	else{
	    print $fhout "$lines[0]\t$lines[11]\t-\n";
	}
	#print $fhout "$_";
	$QueryN=$lines[0];
    }
}
close $fhin;
close $fhout;
