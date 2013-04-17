#!/usr/bin/perl
#===============================================================================
#
#         FILE: SV_Breakdancer_Filter.pl
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
#      VERSION:
#      CREATED: 2013-01-2
#     REVISION:
#===============================================================================
use strict;
use File::Basename;

my ($BD_File,$Sample_Name,$REF_File)=@ARGV;
die "Usage:$0 <Breakdancer_output> <Sample_name> <Reference_name>" if(@ARGV<3);

die "Can't open file $BD_File\n" if(!-e($BD_File));

my %RefSeq;

open(my $fh_ref,"$REF_File") ||  die "Can't open file $REF_File";
my $seqN;
my $seqS;
my @SeqNs;
my @SeqLens;
my $reference=basename($REF_File);
while(<$fh_ref>){
    chomp();
    if(/^>(\S+)/){
        if($seqN ne ""){
            $RefSeq{$seqN}=$seqS;
            push(@SeqLens,length($seqS));
        }
        $seqN=$1;
        push(@SeqNs,$seqN);
        $seqS="";
    }
    else{
        $seqS.=$_;
    }
}
$RefSeq{$seqN}=$seqS;
push(@SeqLens,length($seqS));
close $fh_ref;


my $BD_FILTER_SCORE=90;
my $BD_MERGE_OUT=50;
my $BD_MERGE_IN=100;
my $BD_File_Base=basename($BD_File);

my $BD_CTX_File=$BD_File_Base.".ctx";
my $BD_INV_File=$BD_File_Base.".inv";
my $BD_ITX_File=$BD_File_Base.".itx";
my $BD_DEL_File=$BD_File_Base.".del";
my $BD_INS_File=$BD_File_Base.".ins";

my $Is_CTX=0;
my $Is_INV=0;
my $Is_ITX=0;
my $Is_DEL=0;
my $Is_INS=0;

open(my $fh_ctx,">$BD_CTX_File");
open(my $fh_inv,">$BD_INV_File");
open(my $fh_itx,">$BD_ITX_File");
open(my $fh_del,">$BD_DEL_File");
open(my $fh_ins,">$BD_INS_File");

open(my $fh_BD,"$BD_File");
while(<$fh_BD>){
    next if (/^#/);
    my @lines=split(/\t/,$_);
    next if($lines[8]<$BD_FILTER_SCORE);
    if($lines[6] eq "CTX") {print $fh_ctx "$lines[0]\t$lines[1]\t$lines[3]\t$lines[4]\t$lines[6]\t$lines[7]\t$lines[8]\t$lines[9]\n";$Is_CTX=1;}
    elsif ($lines[6] eq "INV") {print $fh_inv "$lines[0]\t$lines[1]\t$lines[3]\t$lines[4]\t$lines[6]\t$lines[7]\t$lines[8]\t$lines[9]\n";$Is_INV=1;}
    elsif ($lines[6] eq "ITX") {print $fh_itx "$lines[0]\t$lines[1]\t$lines[3]\t$lines[4]\t$lines[6]\t$lines[7]\t$lines[8]\t$lines[9]\n";$Is_ITX=1;}
    elsif ($lines[6] eq "DEL") {print $fh_del "$lines[0]\t$lines[1]\t$lines[3]\t$lines[4]\t$lines[6]\t$lines[7]\t$lines[8]\t$lines[9]\n";$Is_DEL=1;}
    elsif ($lines[6] eq "INS") {print $fh_ins "$lines[0]\t$lines[1]\t$lines[3]\t$lines[4]\t$lines[6]\t$lines[7]\t$lines[8]\t$lines[9]\n";$Is_INS=1;}
}
close $fh_BD;
close $fh_ctx;
close $fh_inv;
close $fh_itx;
close $fh_del;
close $fh_ins;

if($Is_CTX){
    `cat $BD_CTX_File | sort -k 1,1 -k 2,2n -k 3,3 -k 4,4n > $BD_CTX_File.tmp`;
    my $mergefile=&MergeFile("$BD_CTX_File.tmp");
    &VCFConvert($mergefile,"$BD_CTX_File.vcf");
    unlink "$BD_CTX_File.tmp";
    unlink $mergefile;
}
unlink $BD_CTX_File;

if($Is_INV){
    `cat $BD_INV_File | sort -k 1,1 -k 2,2n -k 3,3 -k 4,4n > $BD_INV_File.tmp`;
    my $mergefile=&MergeFile("$BD_INV_File.tmp");
    &VCFConvert($mergefile,"$BD_INV_File.vcf");
    unlink "$BD_INV_File.tmp";
    unlink $mergefile;
}
unlink $BD_INV_File;

if($Is_INS){
    `cat $BD_INS_File | sort -k 1,1 -k 2,2n -k 3,3 -k 4,4n > $BD_INS_File.tmp`;
    my $mergefile=&MergeFile("$BD_INS_File.tmp");
    &VCFConvert($mergefile,"$BD_INS_File.vcf");
    unlink "$BD_INS_File.tmp";
    unlink $mergefile;
}
unlink $BD_INS_File;

if($Is_DEL){
    `cat $BD_DEL_File | sort -k 1,1 -k 2,2n -k 3,3 -k 4,4n > $BD_DEL_File.tmp`;
    my $mergefile=&MergeFile("$BD_DEL_File.tmp");
    &VCFConvert($mergefile,"$BD_DEL_File.vcf");
    unlink "$BD_DEL_File.tmp";
    unlink $mergefile;
}
unlink $BD_DEL_File;

if($Is_ITX){
    `cat $BD_ITX_File | sort -k 1,1 -k 2,2n -k 3,3 -k 4,4n > $BD_ITX_File.tmp`;
    my $mergefile=&MergeFile("$BD_ITX_File.tmp");
    &VCFConvert($mergefile,"$BD_ITX_File.vcf");
    unlink "$BD_ITX_File.tmp";
    unlink $mergefile;
}
unlink $BD_ITX_File;


sub MergeFile{
    my $Filename=shift;
    my $mergefile=$Filename.".mrg";
    open (my $fh_filein,$Filename) || die "Can't open file $Filename\n";
    open (my $fh_mergefile,">$mergefile");
    my $lastchrom;
    my $lastsite;
    my $lastmatchchrom;
    my $lastmatchsite;
    my $lastsupportreads;
    while(<$fh_filein>){
        chomp();
        my @lines=split(/\t/,$_);
        if($lines[0] eq $lastchrom){
            if($lines[2] eq $lastmatchchrom){
                if($Is_CTX == 0){
                    if($lines[1] >= $lastsite-$BD_MERGE_OUT && $lines[1] <= $lastsite+$BD_MERGE_IN && $lines[3] >= $lastmatchsite-$BD_MERGE_IN && $lines[3] <= $lastmatchsite+$BD_MERGE_OUT){
                        if($lines[3] > $lastmatchsite){
                            $lastmatchsite=$lines[3];
                        }
                        $lastsupportreads+=$lines[7];
                        next;                         
                    }
                }
                else{
                    if($lines[1] >= $lastsite-$BD_MERGE_OUT && $lines[1] <= $lastsite+$BD_MERGE_IN && $lines[3] >= $lastmatchsite-$BD_MERGE_IN && $lines[3] <= $lastmatchsite+$BD_MERGE_IN){
                        $lastsupportreads+=$lines[7];
                        next;                           
                    }
                }
            }
        }
        $lastchrom=$lines[0];
        $lastsite=$lines[1];
        $lastmatchchrom=$lines[2];
        $lastmatchsite=$lines[3];
        if($lastsupportreads ne ""){
            print $fh_mergefile "$lastsupportreads\n";
        }
        print $fh_mergefile "$lastchrom\t$lastsite\t$lastmatchchrom\t$lastmatchsite\t$lines[4]\t";
        $lastsupportreads=$lines[7];
    }
    print $fh_mergefile "$lastsupportreads\n";
    close $fh_filein;
    close $fh_mergefile;
    return $mergefile;
}

sub VCFConvert{
    my ($infile,$outfile)=@_;
    open(my $fh_outfile,">$outfile");
    {
my $VCF_Head = << "VCF_HEADER";
##fileformat=VCFv4.1
##source=breakdancer
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=>
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=NTLEN,Number=.,Type=Integer,Description="Number of bases inserted in place of deleted code">
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakpoint">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS,Description="Insertion of novel sequence">
##ALT=<ID=INV,Description="Inverstion">
##ALT=<ID=BND,Description="Breakends">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=SR,Number=1,Type=Integer,Description="Number of supporting reads">
VCF_HEADER
    print $fh_outfile $VCF_Head;
    }

    for(my $i=0;$i<@SeqNs;$i++){
        print $fh_outfile "\#\#contig=\<ID=$SeqNs[$i],length=$SeqLens[$i]\>\n";
    }
    print $fh_outfile "\#\#reference=$reference\n";

    print $fh_outfile "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$Sample_Name\n";

    open(my $fh_infile,"$infile");
    my $SV_count=0;
    while(<$fh_infile>){
        chomp();
        my @lines=split(/\t/,$_);
        if($lines[4] ne "ITX" && $lines[4] ne "CTX"){
            my $length=$lines[3]-$lines[1];
            if($lines[4] eq "DEL"){
                print $fh_outfile "$lines[0]\t$lines[1]\t.\t.\t\<$lines[4]\>\t.\tPASS\tEND=$lines[3];SVTYPE=$lines[4];SVLEN=\-$length\tDP\t$lines[5]\n";
            }
            elsif($lines[4] eq "INV"){
                print $fh_outfile "$lines[0]\t$lines[1]\t.\t.\t\<$lines[4]\>\t.\tPASS\tEND=$lines[3];SVTYPE=$lines[4];SVLEN=$length\tDP\t$lines[5]\n";
            }
            else{
                print $fh_outfile "$lines[0]\t$lines[1]\t.\t.\t\<$lines[4]\>\t.\tPASS\tEND=$lines[3];SVTYPE=$lines[4]\tDP\t$lines[5]\n";
            }
        }
        else{
            $SV_count++;
            my $RefA=getRef($lines[0],$lines[1]);
            my $RefB=getRef($lines[2],$lines[3]);
            my $mateidA=$SV_count;
            my $mateidB=$SV_count+1;
            print $fh_outfile "$lines[0]\t$lines[1]\tbnd_$mateidA\t$RefA\t$RefA\]$lines[2]\:$lines[3]\]\t.\tPASS\tSVTYPE=BND;MATEID=bnd_$mateidB\tDP\t$lines[5]\n";
            print $fh_outfile "$lines[2]\t$lines[3]\tbnd_$mateidB\t$RefB\t$RefB\]$lines[0]\:$lines[1]\]\t.\tPASS\tSVTYPE=BND;MATEID=bnd_$mateidA\tDP\t$lines[5]\n";
            $SV_count++;
        }
    }
    close $fh_infile;
    close $fh_outfile;
}

sub getRef{
    my ($chrom,$start)=@_;
    die "Reference file does not match the input file\n" if(!exists($RefSeq{$chrom}));
    my $seq=$RefSeq{$chrom};
    my $refseq=substr($seq,$start+1,1);
    return $refseq;
}
