#!/usr/bin/perl
use strict;

## This program is used to remove duplicate entry in the fasta file, it is usually used for remove isoforms from transcripts and proteins file.
my ($fasfile, $idfield, $outfile)=@ARGV;
if(@ARGV<2){
    print "Usage:$0 fasta_file id_field <output_file>\n";
    ## id_field is the unique id for each genes.
    exit(0);
}

if($outfile eq ""){
    $outfile="$fasfile.flt";
}
my %geneinfo;
my $seq;
my $gid;
my $seqN;
open (my $fh_fasfile,$fasfile) || die "Can't open file: $fasfile\n";
while(<$fh_fasfile>){
    if(/^>/){
        if($gid ne ""){
            if (!exists($geneinfo{$gid})){
                my %gi;
                $gi{len}=length($seq);
                $gi{seq}=$seq;
                $gi{gname}=$seqN;
                $geneinfo{$gid}=\%gi;
            }
            else{
                if(length($seq)>$geneinfo{$gid}->{len}){
                    my %gi;
                    $gi{len}=length($seq);
                    $gi{seq}=$seq;
                    $gi{gname}=$seqN;
                    $geneinfo{$gid}=\%gi;
                }
            }
            $seq="";
        }
        $seqN=$_;
        s/^\>\s*//;
        s/\s+/ /g;
        s/\s*\|\s*/\|/g;
        my @gn = split(/[\s\|]/);
        $gid=$gn[$idfield-1];
    }
    else{
        $seq=$seq.$_;
    }
}
if (!exists($geneinfo{$gid})){
    my %gi;
    $gi{len}=length($seq);
    $gi{seq}=$seq;
    $gi{gname}=$seqN;
    $geneinfo{$gid}=\%gi;
    print "Dup Sequences: $gid\n";
}
else{
    if(length($seq)>$geneinfo{$gid}->{len}){
        my %gi;
        $gi{len}=length($seq);
        $gi{seq}=$seq;
        $gi{gname}=$seqN;
        $geneinfo{$gid}=\%gi;
    }
}

open (my $fh_out, ">$outfile");
foreach my $key (keys %geneinfo){
    my $seqn=$geneinfo{$key}->{gname};
    my $seqs=$geneinfo{$key}->{seq};
    print $fh_out "$seqn";
    print $fh_out "$seqs";
}
close $fh_out;
exit(1);