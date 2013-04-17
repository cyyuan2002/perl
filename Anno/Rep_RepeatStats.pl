#!/usr/bin/perl -w

#===============================================================================
#
#         FILE: Rep_RepeatStats.pl
#
#        USAGE: 
#
#  DESCRIPTION:
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Yuan Chen
#      COMPANY: Division of Infectious Disease, DUMC
#      VERSION: 1.0
#      CREATED: 
#     REVISION:
#===============================================================================

use strict;
my $RepFile=shift;

my %repeatcounts;
my %repeatslength;
my %repeatstype;

open(my $fh_repfile,"$RepFile") || die "Can't open file: $RepFile\n";
while(<$fh_repfile>){
    $_=&trim($_);
    next if(!/^\d/);
    my @lines=split(/\s+/,$_);
    next if($lines[10] eq "Simple_repeat" || $lines[10] eq "Low_complexity");
    my $replen=$lines[6]-$lines[5]+1;
    if(exists($repeatcounts{$lines[9]})){
        $repeatcounts{$lines[9]}++;
        $repeatslength{$lines[9]}+=$replen;
    }
    else{
        $repeatcounts{$lines[9]}=1;
        $repeatslength{$lines[9]}=$replen;
        $repeatstype{$lines[9]}=$lines[10];
    }
}
close $fh_repfile;

print "Name\tType\tCount\tTotal_length\n";
foreach my $repID (sort {$repeatcounts{$b} <=> $repeatcounts{$a}} keys %repeatcounts){
    print "$repID\t$repeatstype{$repID}\t$repeatcounts{$repID}\t$repeatslength{$repID}\n";
}

exit(0);

sub trim($)
{
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}