#!/home/SCE/chenyuan/bin/Perl-5.12/bin/perl

use strict;
use threads;
use threads::shared;
use Thread::Queue;
use Benchmark;
my $TT0 = new Benchmark;
use Getopt::Long;
my %opts;
#my $version="1.0 Alvin Chen 2011-02-25";
my $version="1.1 Alvin Chen 2011-02-27";

GetOptions(\%opts,"p=s","d=s","r=s","t=s","o=s","n:i","help");
if((!defined $opts{p}) || (!defined $opts{d}) ||(!defined $opts{r})||(!defined $opts{t}) || (!defined $opts{o})){
    &Usage();
}

my $parProteinFile=$opts{p};
my $parDnaFile=$opts{d};
my $parRegionFile=$opts{r};
my $parOutfile=$opts{o};
my $parTempDir=$opts{t};
my $parThreadNum=(defined $opts{n}) ? $opts{n} : 8;


my %proSeq;
my @dnaRegion;
my %dnaSeq;

my $iserror=0;

&Filecheck($parProteinFile,0);
&Filecheck($parDnaFile,0);
&Filecheck($parRegionFile,0);

if($iserror==1){
    exit(1);
}

if(!(-e $parTempDir)){
    mkdir $parTempDir;
}

{
    print << "ProgramRun";

	$0 Running:

	--------------------------------------
	    Protein File:\t$parProteinFile
	    DNA File:\t$parDnaFile
	    DNA Region:\t$parRegionFile
	    Output File:\t$parOutfile
	    Temp Directory:\t$parTempDir
	    Thread Number:\t$parThreadNum
	--------------------------------------

ProgramRun
}

print "\tReading Protein Sequences...\n\n";
open(proFile,$parProteinFile);
my $seqname;
my $seq;

while(<proFile>){
    chomp();
    if(/^>(\S+)/){
	if($seqname ne ""){
	    $proSeq{$seqname}=$seq;
	    $seq="";
	}
	$seqname=$1;
    }
    else{
	$seq=$seq.$_;
    }
}
$proSeq{$seqname}=$seq;
close proFile;

print "\tReading DNA Sequences...\n\n";
open(dnaFile,$parDnaFile);
$seqname="";
$seq="";
while(<dnaFile>){
    if(/^>(\S+)/){
	if($seqname ne ""){
	    $dnaSeq{$seqname}=$seq;
	    $seq="";
	}
	$seqname=$1;
    }
    else{
	$seq=$seq.$_;
    }
}
$dnaSeq{$seqname}=$seq;
close dnaFile;

print "\tReading Region Information...\n\n";
open(regionFile,$parRegionFile);
while(<regionFile>){
    my @lines=split(/\t/,$_);
    my %dnareg;
    $dnareg{'N'}=$lines[0];
    $dnareg{'C'}=$lines[1];
    $dnareg{'S'}=$lines[3];
    $dnareg{'D'}=$lines[2];
    $dnareg{'E'}=$lines[4];
    push(@dnaRegion,\%dnareg);
}

close regionFile;
my $data_queue = new Thread::Queue;
my $result_queue = new Thread::Queue;
my $processing_count :shared = 0;
my $MAX_THREADS = $parThreadNum;
my $request_times=scalar(@dnaRegion);

my @filename;

print "\tBegin Running Genewise ... \n\n";
for (my $n = 0; $n < $MAX_THREADS; $n++)
{
    threads->create(\&thread_io);
}

open(OutFile,">$parOutfile");
my $num=0;


foreach my $data(1..$request_times)
{
    if(($data%500)==0){
	print "\t$data Sequences Done...\n\n";
    }
    if ($data_queue ->pending() > $MAX_THREADS * 2)
    {
	select(undef, undef, undef, 0.02);
        redo;
    }

    $data_queue->enqueue($data);

    if ($result_queue->pending() > 0)
    {
        while (my $result = $result_queue->dequeue_nb())
        {
            if($result) {
		 my $edtres=&editRegion($result);
		 print OutFile "$edtres";
		 #print "$edtres";
	    }
	    else {print "$num\tFailed\n";}
	    $num++;
        }
    }
}

while ($processing_count > 0 or $data_queue->pending() > 0 or $result_queue->pending() > 0)
{
    select(undef, undef, undef, 0.02);
    while (my $result = $result_queue->dequeue_nb())
    {

	if($result) {
	    my $edtres=&editRegion($result);
	    print OutFile "$edtres";
	    #print "$edtres";
	}
        else { print "$num\tFailed!\n"; }
	$num++;
    }
}

foreach my $thread (threads->list(threads::all))
{
    $thread->detach();
}

close OutFile;

my $TT1 = new Benchmark;
my $td = Benchmark::timediff($TT1, $TT0);
$td = Benchmark::timestr($td);
my ($sec) = ($td =~ /(\d+).*/);
my $speed = sprintf("%0.1f",$request_times/$sec);
print "\tAll Job Finished\n\n";
print "\t---------------------------------------\n\n";
print "\tTime expend: $td\n\tAverage Speed: $speed Times Per Second\n\n";
exit(0);

sub editRegion()
{
    my $res=shift;
    my $editres;
    my @resinfo=split(/\%/,$res);
    my $seqnum=$resinfo[0]-1;
    my @wiseline=split(/\n/,$resinfo[1]);
    for(my $i;$i<@wiseline;$i++){
	if(!($wiseline[$i]=~/\/\//)){
	    my @infors=split(/\t/,$wiseline[$i]);
	    my $reginfo=$dnaRegion[$seqnum];
	    $infors[0]=$$reginfo{'C'};
	    if($$reginfo{'D'} eq "+"){
		$infors[3]=$$reginfo{'S'}+$infors[3]-1;
		$infors[4]=$$reginfo{'S'}+$infors[4]-1;
	    }
	    else{
		$infors[3]=$$reginfo{'E'}-$infors[3]+1;
		$infors[4]=$$reginfo{'E'}-$infors[4]+1;
	    }
	    $infors[6]=$$reginfo{D};
	    $infors[8]=~s/\@\d+//g;
	    $editres=$editres.join("\t",@infors)."\n";
	}
	else{
	    $editres=$editres."//\n";
	}
    }
    return $editres;
}


sub thread_io()
{
    while (my $data = $data_queue->dequeue())
    {
        {
            lock $processing_count;
            ++$processing_count;
        }

	my $file=&CreateFile($data-1);
	my $reginfo=$dnaRegion[$data-1];
	my $dir=$$reginfo{'D'};
	my $result=`genewise -gff -silent -quiet -tfor $file`;

	if($result){
	    $result="$data\%$result";
	    $result_queue->enqueue($result);
	    unlink (split(/\s/,$file));
	}
	else{
	    $data_queue->dequeue($data);
	}

        {
            lock $processing_count;
            --$processing_count;
        }
    }
}


sub CreateFile()
{
    my $seqcount=shift;
    my $profile="$parTempDir/$seqcount.pro";
    my $dnafile="$parTempDir/$seqcount.dna";
    my $filename="$profile $dnafile";
    my $reginfo=$dnaRegion[$seqcount];
    my @seqid=split(/\@/,$$reginfo{'N'});

    open(my $pro,">$profile");
    print $pro ">$seqid[0]\n$proSeq{$seqid[0]}\n";
    close $pro;

    open(my $dna,">$dnafile");
    print $dna ">$$reginfo{'N'}\n$dnaSeq{$$reginfo{'N'}}";
    close $dna;

    return "$profile $dnafile";
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
	    print stderr "File $Filename already exists\n";
	    $iserror=1;
	}
    }
}

sub Usage(){
  print << "    Usage";

	Usage:  $0 (version $version)

	<options>
		-p     Protein sequence file
		-d     Dna sequence file
		-r     Gene region table of dna sequence
		-t     Temprary directory
		-o     Output file
		-n     Threads numbers (default:8)

    Usage
	exit(0);
}
