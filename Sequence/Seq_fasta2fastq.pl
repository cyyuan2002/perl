# This program is used to merge fasta file and quality file into fastq file
#!/usr/bin/perl
use strict;

use Getopt::Long;

#----------------user option--------------

my %opts;

GetOptions(\%opts,"f:s","q:s","o:s","help");

if(!defined($opts{f}) && !defined($opts{q}) && !defined($opts{o}))
{
    Usage();
    exit();
}

my $fasfile = $opts{f};
my $qualfile = $opts{q};
my $output = $opts{o};


my %seqs;
my $seq;
my $seqname;
open (fasfile,$fasfile) || die "Can't open fasta file: $fasfile\n";
while(<fasfile>){
    chomp();
    if(/^>(\S+)\s\S+/){
        if($seqname ne ""){
            $seqs{$seqname}=$seq;
            $seq="";
        }
        $seqname=$1;
    }
    else{
        $seq=$seq.$_;
    }
}
$seqs{$seqname}=$seq;
close fasfile;

open (qualfile,$qualfile) || die "Can't open quality file: $qualfile\n";
open (output,">$output");
my $score;
$seqname="";
while (<qualfile>){
    chomp();
    if(/^>(\S+)\s\S+/){
        if($seqname ne ""){
            print output "\@$seqname\n$seqs{$seqname}\n+\n";
            my @scores=split(/\s/,$score);
            my $scorechar;
            foreach my $key (@scores){
                if ($key ne ""){
                    my $char=chr(64+$key);
                    $scorechar=$scorechar."$char";
                }
            }
            print output "$scorechar\n";
            $score="";
        }
        $seqname=$1;
    }
    else{
        $score=$score.$_;
    }
}
close qualfile;
print output "\@$seqname\n$seqs{$seqname}\n+\n";
my @scores=split(/\s/,$score);
my $scorechar;
foreach my $key (@scores){
    if ($key ne ""){
        my $char=chr(64+$key);
        $scorechar=$scorechar."$char";
    }
}
print output "$scorechar\n";

sub Usage{
    print "Usage:Program -f fasta -q quality -o output\n";
}