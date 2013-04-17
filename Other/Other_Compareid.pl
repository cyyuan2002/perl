#!/usr/bin/perl
use strict;

my ($file1,$field1,$file2,$field2,$mode)=@ARGV;
if(@ARGV<3){
    print "Usage:$0 id_file1 field1 id_file2 field2 mode(1 different|2 same)\n";
    exit(1)
}
my %ids;
open(my $fh_file1,$file1) || die "Can't open file $file1\n";
while(<$fh_file1>){
    chomp();
    my @lines=split(/\t/,$_);
    $ids{$lines[$field1]}=1;
}
close $fh_file1;

open(my $fh_file2,$file2) || die "Can't open file $file2\n";
while(<$fh_file2>){
    chomp();
    my @lines=split(/\t/,$_);
    if($mode==1){
		if(!exists($ids{$lines[$field2]})){
		    print "$_\n";
		}
    }
    else{
		if(exists($ids{$lines[$field2]})){
		    print "$_\n";
		}
    }
}
close $fh_file2;
exit(0);
