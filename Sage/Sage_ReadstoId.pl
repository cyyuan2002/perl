#!/usr/bin/perl
use strict;

my ($listfile,$readsfile)=@ARGV;

if(@ARGV<2){
	    print "Usage:$0 Listfile Readsfile\n";
	        exit;
}

open (listfile, $listfile) || die "Can't open file: $listfile\n";
print stderr "Reading Listfile...";
my %seqid;
while(<listfile>){
	    chomp();
	        my @lines=split(/\t/,$_);
		    $seqid{$lines[1]}=$lines[0];
}
close listfile;
print stderr "Reading finished!\n";

open(readsfile,$readsfile) || die "Can't open file: $readsfile\n";
while(<readsfile>){
	    chomp();
	        my @lines=split(/\t/,$_);
		    if(!exists($seqid{$lines[0]})){
			            print stderr "$can't find $lines[0]\n";
				        }
		        else{
				        print "$seqid{$lines[0]}\t$lines[1]\n";
					    }
}
close readsfile;
