#!/usr/bin/perl
#===============================================================================
#
#         FILE:
#
#        USAGE:
#
#  DESCRIPTION: This program is used to calculate GC coverage of genome, which is used for circos
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Yuan Chen
#      COMPANY: MMRL, Duke University Medical Center
#      VERSION: 1.0
#      CREATED:
#     REVISION:
#===============================================================================

use strict;
my ($infile,$windowsize)=@ARGV;
if(@ARGV<1){
    print "Usage:$0 sequence_file <window_size(default:1000bp)>";
    exit(0);
}

$windowsize||=1000;

my $seqid;
my $seq;
open(my $fh_infile,$infile) || die "$!\n";
while(<$fh_infile>){
    chomp();
    if(/^>(\S+)/){
        if($seqid ne ""){
            &calGC($seqid,$seq);
        }
        $seqid=$1;
        $seq="";
    }
    else{
        $seq.=$_;
    }
}
close $fh_infile;
&calGC($seqid,$seq);
exit(1);

sub calGC{
    my ($seqN,$seqS)=@_;
    my $seqlen=length($seqS);
    my $a=$seqlen%$windowsize;
    my $lengthcount;
    if($a==0){
        $lengthcount=int($seqlen/$windowsize);
    }
    else{
        $lengthcount=int($seqlen/$windowsize)+1;
    }
    for(my $i=0;$i<$lengthcount;$i++){
        my $startpos=$i*$windowsize;
        my $subseq=substr($seqS,$startpos,$windowsize);
        my $sublen=length($subseq);
        my $gcnumer=$subseq=~s/[gcGC]//g;
        my $gccontent=sprintf("%.2f",$gcnumer/$sublen);
        my $endpos=$startpos+$sublen-1;
        print "$seqN $startpos $endpos $gccontent\n";
    }
}
