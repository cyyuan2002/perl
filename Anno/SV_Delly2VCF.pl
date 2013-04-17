#!/usr/bin/perl
#===============================================================================
#
#         FILE: SV_Delly2VCF.pl
#
#        USAGE:
#
#  DESCRIPTION: This script is used to convert Delly result to VCF
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Yuan Chen
#      COMPANY: Division of Infectious Disease, DUMC
#      VERSION:
#      CREATED: 12/25/2012
#     REVISION:
#===============================================================================


use strict;
use Getopt::Long;
use File::Basename;

my %opts;
GetOptions(\%opts,"i=s","t=s","r=s","n=s","o:s","s:i","q:i");
if($opts{i} eq "" || $opts{t} eq "" || $opts{r} eq "" || $opts{n} eq ""){
    &Usage();
}

my $Inputfile=$opts{i};
my $Filetype=$opts{t};
my $Samplename=$opts{n};
my $Reffile=$opts{r};
my $Outfile=$opts{o};
my $Qualfilter=$opts{q};
my $Depthfilter=$opts{s};

$Outfile ||= $Inputfile.".vcf";
$Qualfilter ||= 20;
$Depthfilter ||=3;

my $Errorstat=0;
&CheckFile($Inputfile);
&CheckFile($Reffile);

exit(1) if($Errorstat==1);

my %RefSeq;
open(my $fh_ref,"$Reffile");
my $seqN;
my $seqS;
my @SeqNs;
my @SeqLens;
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

open(my $fh_out,">$Outfile");
my $reference=basename($Reffile);
{
my $VCF_Head = << "VCF_HEADER";
##fileformat=VCFv4.1
##source=delly
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=NTLEN,Number=.,Type=Integer,Description="Number of bases inserted in place of deleted code">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=SR,Number=1,Type=Integer,Description="Number of supporting reads">
VCF_HEADER
print $fh_out $VCF_Head;
}

for(my $i=0;$i<@SeqNs;$i++){
    print $fh_out "\#\#contig=\<ID=$SeqNs[$i],length=$SeqLens[$i]\>\n";
}

print $fh_out "\#\#reference=$reference\n";
if($Filetype ne "INV"){
    print $fh_out "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$Samplename\n";
}
else{
     print $fh_out "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n";
}

open (my $fh_delly,"$Inputfile");
while(<$fh_delly>){
    chomp();
    my @lines=split(/\t/,$_);
    if($Filetype ne "INV"){
        next if($lines[5] < $Qualfilter || $lines[4] < $Depthfilter);
    }
    else{
        next if($lines[5] < $Qualfilter );
    }
        #my $refSeq=&getRef($lines[0],$lines[1],$lines[2]);
    my $length=$lines[2]-$lines[1];
    if($Filetype eq "DEL"){
	$length="-$length";
    }
    if($Filetype ne "INV"){
        print $fh_out "$lines[0]\t$lines[1]\t.\t.\t<$Filetype>\t$lines[5]\tPASS\tEND=$lines[2];SVTYPE=$Filetype;SVLEN=$length\tDP\t$lines[4]\n";
    }
    else{
        print $fh_out "$lines[0]\t$lines[1]\t.\t.\t<$Filetype>\t$lines[5]\tPASS\tEND=$lines[2];SVTYPE=$Filetype;SVLEN=$length\n";
    }

}
close $fh_delly;

sub getRef{
    my ($chrom,$start,$end)=@_;
    die "Reference file does not match the input file\n" if(!exists($RefSeq{$chrom}));
    my $seq=$RefSeq{$chrom};
    my $seqlength=$end-$start+1;
    my $refseq=substr($seq,$start+1,$seqlength);
    return $refseq;
}

sub CheckFile{
    my $filename=shift;
    if(!-e $filename){
        print stderr "Can't find file: $filename\n";
        return 1;
        $Errorstat=1 if($Errorstat==0);
    }
    return 0;
}

sub Usage {#help subprogram
    print << "    Usage";

	Usage: -i <Delly Output> -t <File Type (DEL/DUP/INV)> -r <Reference File> -n <Sample Name>[options]

        Options:

                -o     Output File

                -q     Mapping quality filter (default: >= 20)

                -s     Supporting reads number (default: >= 3)

    Usage

    exit(1);
};
