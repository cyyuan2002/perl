#!/usr/bin/perl
#This file is used to extract scaffold sequences from the assembly file of SOAPdenovo

use strict;

my $file=shift;

die "Usage:$0 <scafseq_file>" if(@ARGV<1);

my $isoutput=0;

open(my $fh_file,$file) || die "Can't open file $file\n";
while(<$fh_file>){
    if(/^>/){
        if(/scaffold/){
            print "$_";
            $isoutput=1;
        }
        else{
            $isoutput=0;
        }
    }
    else{
        print "$_" if($isoutput);
    }
}
close $fh_file;
exit(1);

