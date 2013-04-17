#!/usr/bin/perl -w

#===============================================================================
#
#         FILE: Seq_Change_ScafName.pl
#
#        USAGE: Seq_Change_ScafName.pl <-c coords> <-a assembled.fas> <-f ref.fas> [-o out.fas]
#
#  DESCRIPTION:This program is used to changed assembly contig/scaf names
#              according to the position of the reference genome.
#              The coords file is obtained by mummer show-coords -THrcl
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
use Getopt::Long;
my %opts;

GetOptions(\%opts,"c=s","a=s","r=s","o:s","help");

if((!defined $opts{c})||(!defined $opts{a})||(!defined $opts{r})){
    die "Usage:$0 <-c coords> <-a scaf.fas> <-r ref.fas> [-o out.fas]\n";
}

my $parCoordsfile=$opts{c};
my $parAssemfile=$opts{a};
my $parReffile=$opts{r};
my $parOutfile=$opts{o};

my @refseqID;
my %seqs;
my %orders;
my %scafdirection;
my %scafbelongs;

$parOutfile ||= "$parAssemfile.new";

open(my $fh_reffile,$parReffile) || die "Can't open file: $parReffile\n";
while(<$fh_reffile>){
    chomp();
    if(/^>(\S+)/){
        push(@refseqID,$1);
    }
}
close $fh_reffile;

my $lastseqN="";
my $lastseq;
open(my $fh_assemfile,$parAssemfile) || die "Can't open file: $parAssemfile\n";
while(<$fh_assemfile>){
    chomp();
    if(/^>(\S+)/){
        if($lastseqN ne ""){
            $seqs{$lastseqN}=$lastseq;
        }
        $lastseqN=$1;
        $lastseq="";
    }
    else{
        $lastseq.=$_;
    }
}
$seqs{$lastseqN}=$lastseq;
close $fh_assemfile;

my $REFCOL=11;
my $ASSEMCOL=12;
my $lastchrom="";
my @assembleids;

open (my $fh_coordsfile,$parCoordsfile) || die "Can't open file: $parCoordsfile\n";
my %scafmaplengths;

while(<$fh_coordsfile>){
    chomp();
    ##--scaf order 
    my @lines=split(/\t/,$_);
    if($lines[$REFCOL] ne $lastchrom){
        if(@assembleids>0){
            my @arraytemp=@assembleids;
            $orders{$lastchrom}=\@arraytemp;
        }
        @assembleids=();
        $lastchrom=$lines[$REFCOL];
    }
    push(@assembleids,$lines[$ASSEMCOL]);
    
    my $length=abs($lines[3]-$lines[2]);
    ##--scaf map length---
    if(!exists($scafmaplengths{$lines[$ASSEMCOL]}->{$REFCOL})){
        $scafmaplengths{$lines[$ASSEMCOL]}->{$lines[$REFCOL]}=$length;
    }
    else{
        $scafmaplengths{$lines[$ASSEMCOL]}->{$lines[$REFCOL]}+=$length;
    }
    
    ##--scaf map direction
    if(!exists($scafdirection{$lines[$ASSEMCOL]})){
        if($lines[2]<$lines[3]){
            $scafdirection{$lines[$ASSEMCOL]}->{'+'}=$length;
            $scafdirection{$lines[$ASSEMCOL]}->{'-'}=0;
        }
        else{
            $scafdirection{$lines[$ASSEMCOL]}->{'+'}=0;
            $scafdirection{$lines[$ASSEMCOL]}->{'-'}=$length;
        }
    }
    else{
        if($lines[2]<$lines[3]){
            $scafdirection{$lines[$ASSEMCOL]}->{'+'}+=$length;
        }
        else{
             $scafdirection{$lines[$ASSEMCOL]}->{'-'}+=$length;
        }
    }
}
my @arraytemp=@assembleids;
$orders{$lastchrom}=\@arraytemp;
close $fh_coordsfile;

foreach my $scaf(keys %scafmaplengths){
    my $maxlength;
    foreach my $refID(keys %{$scafmaplengths{$scaf}}){
        if(!exists($scafbelongs{$scaf})){
            $scafbelongs{$scaf}=$refID;
            $maxlength=$scafmaplengths{$scaf}->{$refID};
        }
        else{
            if($scafmaplengths{$scaf}->{$refID} > $maxlength){
                $scafbelongs{$scaf}=$refID;
                $maxlength=$scafmaplengths{$scaf}->{$refID};
            }
        }
    }
}



my %changedIDs;
my $scafcount=1;

open(my $fh_fileout,">$parOutfile");
foreach my $chrom (@refseqID){
    my @ids=@{$orders{$chrom}};
    foreach my $scaf(@ids){
        if(!exists($changedIDs{$scaf})){
            next if($scafbelongs{$scaf} ne $chrom);
            $scaf=~/([a-zA-Z]+)\d+/;
            my $newid="$1$scafcount";
            my $seqout=$seqs{$scaf};
            if($scafdirection{$scaf}->{'+'} < $scafdirection{$scaf}->{'-'}){
                 $seqout=~tr/atgcATGC/tacgTACG/;
                 $seqout=reverse $seqout;
            }
            print $fh_fileout ">$newid\n$seqout\n";
            $changedIDs{$scaf}=1;
            $scafcount++;
        }
    }
}

foreach my $seqID (keys %seqs){
    if(!exists($changedIDs{$seqID})){
        $seqID=~/(\w+)\d+/;
        my $newid="$1$scafcount";
        print $fh_fileout ">$newid\n$seqs{$seqID}\n";
        $changedIDs{$seqID}=1;
        $scafcount++;
    }
}


close $fh_fileout;
exit(1);
