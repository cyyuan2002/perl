#!/usr/bin/perl
use strict;

my ($filein,$fileout)=@ARGV;
if(@ARGV<2){
    print stderr "Usage: $0 Input_File(fastq) Output_File(fasta)\n";
    exit(0);
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
    my $i;
    while($fh->gzreadline(my $line)){
	$i++;
	if($i==1){
	    $line=~s/@//;
	    print fileout ">$line";
	}
	elsif($i==2){
	    print fileout "$line";
	}
	elsif($i==4){
	    $i=0;
	}
    }
    $fh->gzclose();
}
#elsif($filetype eq "bz2"){
#    use Compress::Bzip2;
#    my $fh=bzopen($filein,"r");
#    my $i;
#    while($fh->bzreadline(my $line)){
#	$i++;
#	$line=~s/@//;
#	if($i==1){
#		print fileout ">$line";
#	}
#	elsif($i==2){
#		print fileout "$line";
#	}
#	elsif($i==4){
#		$i=0;
#	}
#    }
#    $fh->gzclose();
#}
else{
    open(file,$filein);
    my $i;
    while(<file>){
	$i++;
	$_=~s/@//;
	if($i==1){
		print fileout ">$_";
	}
	elsif($i==2){
		print fileout "$_";
	}
	elsif($i==4){
		$i=0;
	}
    }
    close file;
}

close fileout;

