#!/usr/bin/perl
#===============================================================================
#
#         FILE: Seq_FqExtract.pl
#
#        USAGE: ./ 
#
#  DESCRIPTION: This program is used to extract partial data from fastq file
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Yuan Chen 
#      COMPANY: Duke University Medical Center, MMRL
#      VERSION: 1.0
#      CREATED: 02/01/2012 10:22:05 AM
#     REVISION: ---
#===============================================================================


#use strict;
use Compress::Zlib;
use Getopt::Long;

my $version="1.0 Alvin Chen 2011-01-15";
my %opts;
##Read parameters
GetOptions(\%opts,"m=i","p=s","a=s","b:s","help");
if((!defined $opts{m})||(!defined $opts{a})||(!defined $opts{p})){
    &Usage();
}

my $parMode=$opts{m};
my $parSize=$opts{p};
my $parFileA=$opts{a};
my $parFileB=$opts{b};

my @FileNA;
my @FileNB;
my $outFileA;
my $outFileB;
my $fhoutfileA;
my $fhoutfileB;
my $parBlocksize=10000;

my %readsA;
my %readsB;

my $iserror=0;

if($parMode==2){
    if($parFileB eq ""){
	print stderr "Lack of fileb for paired-end mode\n";
	exit(1);
    }
}

##Check input file
&Filecheck($parFileA,0);
if($parMode==2){
    &Filecheck($parFileB,0);
}

if($parSize>1){
    print stderr "Error: -p percentage must smaller than 1\n";
    $iserror=1;
}

if($iserror==1){
    exit(1);
}

##Autodetect input filetype, generate file name of output file;
my $fileA=`basename $parFileA`;
@FileNA=split(/\./,$fileA);
pop(@FileNA);
$outFileA=join(".",@FileNA);
$outFileA="$outFileA.part.fq";
open($fhoutfileA,">$outFileA");

if($parMode==2){
	my $fileB=`basename $parFileB`;
    @FileNB=split(/\./,$fileB);
    pop(@FileNB);
    $outFileB=join(".",@FileNB);
    $outFileB="$outFileB.part.fq";
    open($fhoutfileB,">$outFileB");
}


if($parMode==1){
    my $casenum=$parBlocksize*$parSize;
    my ($fhFileA);
    open($fhFileA,"<","$parFileA");
    my %selectAlines=%{&randlines($casenum)};
    my $seqcount=0;
    my $linecount=0;
    my $readsA;
    my $readsB;
    while(my $lineA=<$fhFileA>){
        $linecount++;
        if($linecount<4){
            $readsA=$readsA.$lineA;
        }
        else{
            $linecount=0;
            $readsA=$readsA.$lineA;
            $seqcount++;
            if(exists($selectAlines{$seqcount})){
                print $fhoutfileA $readsA;
            }
            $readsA="";
        }
        if($seqcount==10000){
            $seqcount=0;
            %selectAlines={};
            %selectAlines=%{&randlines($casenum)};
        }
    }
    close $fhFileA;
    close $fhoutfileA;
}
elsif($parMode==2){
    my $casenum=$parBlocksize*$parSize;
    my ($fhFileA,$fhFileB);
    open($fhFileA,"<","$parFileA");
    open($fhFileB,"<","$parFileB");
    my %selectAlines=%{&randlines($casenum)};
    my $seqcount=0;
    my $linecount=0;
    my $readsA;
    my $readsB;
    while(my $lineA=<$fhFileA>){
        my $lineB=<$fhFileB>;
        $linecount++;
        if($linecount<4){
            $readsA=$readsA.$lineA;
            $readsB=$readsB.$lineB;
        }
        else{
            $linecount=0;
            $readsA=$readsA.$lineA;
            $readsB=$readsB.$lineB;
            $seqcount++;
            if(exists($selectAlines{$seqcount})){
                print $fhoutfileA $readsA;
                print $fhoutfileB $readsB;
            }
            $readsA="";
            $readsB="";
        }
        if($seqcount==10000){
            $seqcount=0;
            %selectAlines={};
            %selectAlines=%{&randlines($casenum)};
        }
    }
    close $fhFileA;
    close $fhFileB;
    close $fhoutfileA;
    close $fhoutfileB;
}

exit(0);

sub randlines{
    my $casenum=shift;
    my @array = 0 .. $parBlocksize;
    my $sca = \@array;
    for (my $j = @$sca-1; $j > 0; --$j) {
        my $i = int(rand($j));
        my $t = $sca->[$j]; $sca->[$j] = $sca->[$i]; $sca->[$i] = $t;
    }
    my @newarray=sort{$a<=>$b}(@array[0..$casenum]);
    my %arraynum=map{$_ => 1}@newarray;
    return \%arraynum;
}

sub Filecheck{
    my ($Filename,$mode)=@_;
    if($mode==0){
	if(!(-e $Filename)){
	    print stderr "Can't find file: $Filename\n";
	    $iserror=1;
	}
    }
    else{
	if(-e $Filename){
	    print stderr "Directory $Filename already exists\n";
	    $iserror=1;
	}
    }
}



sub Usage(){
  print << "    Usage";

	This program is used to filter the reads in Fastq format. The compressed gzip file are supported. The pair-end reads
	    files should be in compressed gzip format. The output file name is file_name.cln.gz.

	Usage:  $0 (version $version)

	<options>
		-m     Program mode, 1 single-end File, 2 pair-end Files (must given)
                -p     Percentage of the whole data
		-a     File_name of Fastq file (must given)
		-b     File_name of Fastq file, only used when mode is 2 (must given when mode is 2)

    Usage

	exit(1);
};
