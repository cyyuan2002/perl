#!/usr/bin/perl
#This script is extract the longest transcript of each gene from Trinity output
use strict;

my $TrinityOut=shift;
open(my $fh_fas,$TrinityOut) or die "Can't open file: $TrinityOut";

my %Seqs;
my $SeqName;
my $GeneName;
my %LongestTrans;
my $seq;
while (<$fh_fas>) {
    if (/^>(\S+)/) {
        if ($SeqName ne "") {
            my $seqlength=length($seq);
            $Seqs{$SeqName}=$seq;
            if (!exists($LongestTrans{$GeneName})) {
                $LongestTrans{$GeneName}->{'id'}=$SeqName;
                $LongestTrans{$GeneName}->{'len'}=$seqlength;
            }
            else{
                if ($LongestTrans{$GeneName}->{'len'} < $seqlength) {
                    $LongestTrans{$GeneName}->{'id'}=$SeqName;
                    $LongestTrans{$GeneName}->{'len'}=$seqlength;
                }
                
            }
        }
        $SeqName=$1;
        my @seqn=split(/\_/,$SeqName);
        $GeneName=$seqn[0];
        $seq="";
    }
    else{
        $seq.=$_;
    }
}
{
    my $seqlength=length($seq);
    $Seqs{$SeqName}=$seq;
    if (!exists($LongestTrans{$GeneName})) {
        $LongestTrans{$GeneName}->{'id'}=$SeqName;
        $LongestTrans{$GeneName}->{'len'}=$seqlength;
    }
    else{
        if ($LongestTrans{$GeneName}->{'len'} < $seqlength) {
            $LongestTrans{$GeneName}->{'id'}=$SeqName;
            $LongestTrans{$GeneName}->{'len'}=$seqlength;
        }
        
    }
}
close $fh_fas;

foreach my $gene (keys %LongestTrans){
    my $seqID=$LongestTrans{$gene}->{'id'};
    print ">$gene\n$Seqs{$seqID}";
}