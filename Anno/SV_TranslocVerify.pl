#!/usr/bin/perl

#===============================================================================
#
#         FILE: SV_TranslocVerify.pl
#
#        USAGE:
#
#  DESCRIPTION: This program is used to verify the translocation by mapping
#               the assembled sequence to the assembled genome
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
use File::Temp;
use File::Basename;

my %opts;
GetOptions(\%opts,"i=s","d=s","f:i","g:s","o:s");
if(! defined($opts{i}) || ! defined($opts{d})){
    &Usage();
}

my $parRstFile=$opts{i};
my $parDtlFile=$opts{d};
my $parGenomeFile=$opts{g};
my $parOutFile=$opts{o};
my $parFlankLen=$opts{f};
my $parMinIdentity=98;
my $parMinLengthCov=0.95;

$parOutFile ||= "$parRstFile.out";
$parFlankLen ||= 250;

my $parMode;
if($parGenomeFile eq ""){
    $parMode = 0;
}
else{
    die "Can't open file: $parGenomeFile\n" if (! -e($parGenomeFile));
    $parMode = 1;
}

my %detailInfo;


die "Can't open file: $parRstFile\n" if(! -e ($parRstFile));
my $tempfile="$parRstFile.tmp";

`cat $parRstFile | sort -k 2,2 -k 3,3n -k 5,5 -k 6,6n > $tempfile`;

my @infolines;

open(my $fh_dtlfile,"$parDtlFile") || die "Can't open file $parDtlFile\n";
while(<$fh_dtlfile>){
    chomp();
    next if($_ eq "");
    if(/^\@(\d+)-(\d+)/){
        if(scalar(@infolines)>0){
            &dealinfo(@infolines);
        }
        @infolines=();
    }
    push(@infolines,$_);
}
&dealinfo(@infolines);

my @provedpoints;
open(my $fh_tempfile,"$tempfile");
my $last_chrom_a;
my $last_site_a;
my $last_chrom_b;
my $last_site_b;

while(<$fh_tempfile>){
    chomp();
    my @lines=split(/\t/,$_);
    next if ($lines[12] ne "Approved");
    #print "$_\n";
    my $ID=$lines[0];
    my @inforlist=@{$detailInfo{$ID}};

    if($parMode==0){
        my %infors=%{$inforlist[0]};
        my $site_a;
        my $site_b;
        my $seqlen;
        $site_a=$infors{'site_a'};
        $site_b=$infors{'site_b'};
        $seqlen=length($infors{'seq'});

        if(@inforlist>1){
            for(my $i=1;$i<@inforlist;$i++){
                %infors=%{$inforlist[$i]};
                my $site_n_a=$infors{'site_a'};
                my $site_n_b=$infors{'site_b'};
                my $seqlen_n=length($infors{'seq'});
                if($seqlen_n > $seqlen){
                    $site_a=$site_n_a;
                    $site_b=$site_n_b;
                    $seqlen=$seqlen_n;
                }
            }
        }

        if($last_chrom_a eq $lines[1] && $last_chrom_b eq $lines[4]){
            if(abs($last_site_a - $site_a) <= $parFlankLen && abs($last_site_b -$site_b) <= $parFlankLen){
                pop(@provedpoints);
                my $newpoints="$lines[1]\t$site_a\t$lines[4]\t$site_b\t$lines[7]\n";
                push(@provedpoints,$newpoints);
            }
        }
        else{
            my $newpoints="$lines[1]\t$site_a\t$lines[4]\t$site_b\t$lines[7]\n";
            push(@provedpoints,$newpoints);
        }
        $last_chrom_a=$lines[1];
        $last_chrom_b=$lines[4];
        $last_site_a=$site_a;
        $last_site_b=$site_b;
    }
    if($parMode==1){
        my ($res,$site_a,$site_b)=&verifyfas(\@inforlist);
        if($res==1){
            if($last_chrom_a eq $lines[1] && $last_chrom_b eq $lines[4]){
                if(abs($last_site_a - $site_a) <= $parFlankLen && abs($last_site_b -$site_b) <= $parFlankLen){
                    pop(@provedpoints);
                    my $newpoints="$lines[1]\t$site_a\t$lines[4]\t$site_b\t$lines[7]\n";
                    push(@provedpoints,$newpoints);
                }
            }
            else{
                my $newpoints="$lines[1]\t$site_a\t$lines[4]\t$site_b\t$lines[7]\n";
                push(@provedpoints,$newpoints);
            }
            $last_chrom_a=$lines[1];
            $last_chrom_b=$lines[4];
            $last_site_a=$site_a;
            $last_site_b=$site_b;
        }
    }
}
close $fh_tempfile;
unlink $tempfile;

open (my $fh_outfile,">$parOutFile");

foreach my $point(@provedpoints){
    print $fh_outfile $point;
}
close $fh_outfile;

exit(0);

sub verifyfas{
    my $arrayref=shift;
    my @inforlist=@{$arrayref};
    my $site_a=0;
    my $site_b=0;
    my $seqlen=0;
    my $ispassed=0;

    for(my $i=0;$i<@inforlist;$i++){
        my %infors=%{$inforlist[$i]};
        my $seq=$infors{'seq'};
        my $res=&mapfas($seq);
        if($res==1){
            my $site_n_a=$infors{'site_a'};
            my $site_n_b=$infors{'site_b'};
            my $seqlen_n=length($infors{'seq'});
            if($seqlen_n > $seqlen){
                $site_a=$site_n_a;
                $site_b=$site_n_b;
                $seqlen=$seqlen_n;
            }
            $ispassed=1;
        }
    }

    return ($ispassed,$site_a,$site_b);

}

sub mapfas{
    my $seq=shift;
    my $seqlen=length($seq);
    my $FasTempfile=File::Temp::tempnam(".","fastmp");
    #print "Seq_len:$seqlen\n";
    $FasTempfile=basename($FasTempfile);
    $FasTempfile="$FasTempfile.fas";
    open(my $fh_outfile, ">$FasTempfile");
    print $fh_outfile ">temp\n$seq\n";
    close $fh_outfile;
    my $blastouts=`blastn -query $FasTempfile -subject $parGenomeFile -outfmt 6 -evalue 1e-20`;
    #print "$blastouts\n";
    unlink $FasTempfile;
    my @blastout=split(/\n/,$blastouts);
    my $ispassed=0;

    foreach my $blastline(@blastout){
        my @lines=split(/\t/,$blastline);
        if($lines[2]>=$parMinIdentity && $lines[3]>=($seqlen*$parMinLengthCov)){
            if($ispassed==0){
                $ispassed=1;
            }
            else{
                $ispassed=0;
                last;
            }
        }
    }

    return ($ispassed);
}


sub dealinfo{
    my @lines=@_;
    my $isseq=0;
    my $seq;
    my $ID;
    my $chrom_a;
    my $site_a;
    my $chrom_b;
    my $site_b;

    foreach my $line(@lines){
        if($line=~/^\@(\d+)\-(\d+)\t(\S+)\t(\S+)/){
            $ID=$1;
            my @sites_a=split(/\,/,$3);
            my @sites_b=split(/\,/,$4);
            my $chrom_a=$sites_a[0];
            my $chrom_b=$sites_b[2];
            if($chrom_a ne $chrom_b){
                if($sites_b[0]  eq $chrom_a){
                    $site_a=$sites_b[1];
                    $site_b=$sites_b[3];
                }
                else{
                    $site_a=$sites_b[3];
                    $site_b=$sites_b[1];
                }
            }
            else{
                if($sites_a[1] < $sites_a[3]){
                    if($sites_b[1] < $sites_b[3]){
                        $site_a=$sites_b[1];
                        $site_b=$sites_b[3];
                    }
                    else{
                        $site_a=$sites_b[3];
                        $site_b=$sites_b[1];
                    }
                }
                else{
                    if($sites_b[1] < $sites_b[3]){
                        $site_a=$sites_b[3];
                        $site_b=$sites_b[1];
                    }
                    else{
                        $site_a=$sites_b[1];
                        $site_b=$sites_b[3];
                    }
                }
            }
        }
        if($line=~/^\>(\d+)-(\d+)/){
            $isseq=1;
        }
        else{
            $seq.=$line if($isseq);
        }
    }

    my @infolist;
    my %infors;

    $infors{'site_a'}=$site_a;
    $infors{'site_b'}=$site_b;
    $infors{'seq'}=$seq;
    if(exists($detailInfo{$ID})){
        @infolist=@{$detailInfo{$ID}};
    }
    push(@infolist,\%infors);
    $detailInfo{$ID}=\@infolist;
}


sub Usage #help subprogram
{
    print << "    Usage";

	Usage: $0 -i <Result_File> -d <Detail_File> -g <Genome_File>

    Usage

	exit(1);
};
