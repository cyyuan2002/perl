#!/usr/bin/perl

use strict;
use warnings;

my $seq;
my $seen;# for tracking homopolymers
my $i;# running length of homopolymer
my $max=0;# longest homopolymer found
my @hompol;# hompol = homopolymer counts, e.g. hompol["A"][3] is the count of AAA
my %bases; 
$bases{"A"} = 1;$bases{"C"} = 2;$bases{"G"} = 3;$bases{"T"} = 4;

# set the record separator to the '&gt;' symbol
$/=">";
<>;# remove the empty first 'sequence'
# loop over the sequences
while (<>){
    $seen = "";
    $i=1;
    # remove the trailing '>' symbol
    chomp;
    # split the entry into individual lines based on the newline character
    my @lines = split(/\n/,$_);
    # the header is the first line (now without the '>' symbol)
    my $header = shift @lines;
    # the sequence is the rest
    my $seq = join "", @lines;
    
    # loop over the bases
    for (split (//, uc($seq))) {
        # skip N's
        next if /N/;
        # when the current base is the same as the previous one
        if ($seen eq $_){$i++}
        # otherwise, the homopolymer run has ended
        else {
            $hompol[$bases{$seen}][$i]++ if $seen ne "";
            $max=$i if $i>$max;
            $i=1;
        }
        $seen=$_;
    }
    # do not forget last homopolymer
    $hompol[$bases{$seen}][$i]++ if $seen ne "";
    $max=$i if $i>$max;
}
$/="\n"; # reset the record separator

# output
print "\tA\tC\tG\tT\n";
# homopolymer length loop
for $i (1..$max){
    print $i;
    # bases loop
    for my $j (1..4){
        print "\t";
        print $hompol[$j][$i] || 0;
    }

    print "\n";
}
