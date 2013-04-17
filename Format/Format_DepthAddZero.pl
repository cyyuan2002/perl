#!/usr/bin/perl
#This program is used to full fill coverage information from samtools depth.
#It will add zero before and after the valued coordinates.
use strict;

my ($File_Cov)=@ARGV;
if (@ARGV<1) {
    die "Usage: <Coverage_File>\n";
}

my @cov_info;
my %seq_len;

my $lastpos=0;
my $lastchrom="";
open(my $fh_filein, "<","$File_Cov") or die "Can't open file: $File_Cov\n";
while (<$fh_filein>) {
    chomp();
    my @lines=split(/\t/,$_);
    if ($lastchrom ne $lines[0]) {
        if ($lastchrom ne "") {
            my $pos=$lastpos+1;
            print "$lastchrom\t$pos\t0\n";
        }
        $lastpos=0;
        $lastchrom=$lines[0];
    }
    
    if ($lastpos < $lines[1]-1) {
        if ($lastpos != 0 && $lastpos < $lines[1]-2) {
            my $L_pos=$lastpos+1;
            print "$lines[0]\t$L_pos\t0\n";
        }
        my $pos=$lines[1]-1;
        print "$lines[0]\t$pos\t0\n";
    }
    print "$_\n";
    $lastpos=$lines[1];
}
my $pos=$lastpos+1;
print "$lastchrom\t$pos\t0\n";
close $fh_filein;
