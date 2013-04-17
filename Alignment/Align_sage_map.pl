#!/usr/bin/perl
use strict;

my ($sagefile,$file)=@ARGV;
my $outfile=$file.".map";
my @tags;
my @revtags;

if(@ARGV<2){
    print "Usage:program sagefile reffile\n";
    exit;
}

open (sagefile,"$sagefile") || die "Can't open file $sagefile\n";
while(<sagefile>){
    chomp();
    if(/^>/){
        next;
    }
    push (@tags,$_);
    my $rev_seq=reverse($_);
    $rev_seq=~ tr/[ATGC]/[TACG]/;
    push (@revtags,$rev_seq);
}
close sagefile;

print stderr "All sequences readed!\n";


open (file,$file) || die "Can't open file\n";
open (outfile,">$outfile");
my $seqname;
my $seq;
my $count;
my %seqs;
my $seqcount;
my $count;

while(<file>){
    chomp();
    if(/^>(\S+)/){
        if ($seq ne ""){
            #$seqs{$seqname}=$seq;
            for(my $i=0;$i<@tags;$i++){
                if($seq=~/$tags[$i]/){
                    my $pos=$-[0];
                    print outfile "$tags[$i]\t+\t$seqname\t$pos\n";
                }    
                elsif($seq=~/$revtags[$i]/){
                    my $pos=$-[0];
                    print outfile "$tags[$i]\t-\t$seqname\t$pos\n";
                }
            }
        }
        $seq="";
        $seqname=$1;
        next;
    }
    $seq=$seq.$_;
    $seqcount++;
    $count++;
    if($count==20){
        print "$seqcount sequences searched\n";
        $count=0;
    }
}
$seqs{$seqname}=$seq;
close file;
