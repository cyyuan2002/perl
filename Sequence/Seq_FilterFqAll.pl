#!/usr/bin/perl
#===============================================================================
#
#         FILE: Seq_FilterFqAll.pl
#
#        USAGE: Seq_FilterFqAll.pl <-a fastq_1> <-b fastq_2> <-l lib_size> 
#
#  DESCRIPTION: This program will use filter_data_gz, duplication and corrector
#               to filter the paired-end fastq files
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Yuan Chen
#      COMPANY: MMRL, Duke University Medical Center
#      VERSION: 1.0
#      CREATED: 13/05/2012
#     REVISION:
#===============================================================================

use strict;
use Getopt::Long;
use Cwd;

my %opts;
GetOptions(\%opts,"a=s","s1:i","e1:i","b=s","s2:i","e2:i","u:i","n:i","j=s","p:i","l=i","f:i","q:i","c:i","r:i","k:i","d:i","e:i","help");

if((!defined $opts{a})||(!defined $opts{b})||(!defined $opts{l})||(!defined $opts{j})){
    #print "$opts{a}\t$opts{b}\t$opts{l}\t$opts{j}\n";
    &Usage();
}

my $f1=$opts{a};
my $f2=$opts{b};
my $lib_size=$opts{l};
my $jobName=$opts{j};
my $start1=(defined $opts{s1}) ? $opts{s1}: 0;
my $end1=(defined $opts{e1}) ? $opts{e1} : 0;
my $start2=(defined $opts{s2}) ? $opts{s2} : 0;
my $end2=(defined $opts{e2}) ? $opts{e2} : 0;
my $B_cutoff=(defined $opts{u}) ? $opts{u} : 40;
my $N_num=(defined $opts{n}) ? $opts{n} : 10;
my $Q_shift=(defined ($opts{q})) ? $opts{q} : 64;
my $iscorrection=(defined ($opts{c})) ? $opts{c} : 0;
my $isrmtempfile=(defined ($opts{r})) ? $opts{r} : 1;
my $cpunum=(defined($opts{p})) ? $opts{p} : 4;
my $skcutoff=(defined($opts{k})) ? $opts{k} : 2;
my $ekcutoff=(defined($opts{e})) ? $opts{e} : 2;
my $maxerror=(defined($opts{d})) ? $opts{d} : 2;
my $filetype=(defined($opts{f})) ? $opts{f} : 1;

my @delfiles;

my $filter_gz="~/bin/Filter_data_gz/filter_data_gz";
my $duplication="~/bin/Filter_data_gz/duplication";
my $kmer_freq="~/bin/correction/KmerFreq";
my $correction="~/bin/correction/Corrector";
my $merge_pair="~/bin/correction/merge_pair.pl";
my $stats_perl="~/bin/perl/Sequence/Seq_FqFilterStat.pl";
my $directory=getcwd();

#die "Can't find $filter_gz\n" if (! -e $filter_gz);
#die "Can't find $duplication\n" if (!-e $duplication);
#die "Can't find $stats_perl\n" if (!-e $stats_perl);

my $name = `basename $f1`;
chomp $name;
my $name2 = `basename $f2`;
chomp $name2;
my $tempfile1="$jobName.flt.pbs";

open (my $out,">$tempfile1");

&Titleprint($out);
print $out "echo Begin to run filter_data_gz...\n";
print $out "$filter_gz -y -z -Q $Q_shift -w $N_num -q $filetype -B $B_cutoff -l $lib_size -a $start1 -b $end1 -c $start2 -d $end2 $f1 $f2 $name.reads.stat $name.clean.tmp $name2.clean.tmp && mv $name.clean.tmp $name.clean && mv $name2.clean.tmp $name2.clean";
print $out "\necho Time is `date`\necho Begin to run dupication...\n";
print $out "$duplication $name.clean $name2.clean $name.clean.dup.clean $name2.clean.dup.clean $name.clean.dup.stat\n";
print $out "mv $name.clean.dup.clean $name.flt.dup.cln && mv $name2.clean.dup.clean $name2.flt.dup.cln\n";
print $out "perl $stats_perl $name.reads.stat $name.clean.dup.stat $jobName\n";
if($isrmtempfile==1){
    print $out "rm $name.clean $name2.clean\n";
    print $out "rm $name.reads.stat $name.clean.dup.stat\n";
}
close $out;
`qsub -l highprio $tempfile1`;

if($iscorrection==1){
    my $tempflst=time;
    $tempflst="$tempflst.lst";
    open(my $tmp,">$tempflst");
    print $tmp "$name.flt.dup.cln\n$name2.flt.dup.cln\n";
    close $tmp;
    my $tempfile2=time;
    $tempfile2="$jobName.cor.pbs";
    open(my $pbsout,">$tempfile2");
    &Titleprint($pbsout,$cpunum);
    print $pbsout "echo Building k-mer frequences...\n";
    print $pbsout "$kmer_freq -i $tempflst -o $jobName.kmer -s 17";
    print $pbsout "\necho Time is `date`\necho Correcting errors..\n";
    print $pbsout "$correction -i $tempflst -r $jobName.kmer.freq -s 17 -t $cpunum -k $skcutoff -e $ekcutoff -d maxerror";
    print $pbsout "\necho Time is `date`\necho Making pairs...\n";
    print $pbsout "$merge_pair $name.flt.dup.cln.corr $name2.flt.dup.cln.corr $jobName";
    print $pbsout "\necho All job finished\n";
    if($isrmtempfile==1){
	print $pbsout "rm $tempflst $jobName.kmer.stat\n";
	print $pbsout "rm $jobName.kmer.freq $name.flt.dup.cln.corr $name2.flt.dup.cln.corr\n";
	print $pbsout "rm $tempfile1 $tempfile2\n";
    }
    close $pbsout;
    `qsub -l highprio -W depend=afterok:$tempfile1 $tempfile2`;
}
exit(1);

sub Titleprint(){
    my $outfile=shift;
    my $cpu=shift;
    $cpu ||=1;
    print $outfile "\#!/bin/bash\n\#\$ -N $jobName.flt\n\#\$ -pe threaded $cpu\n\#\$ -q \*\@mmrl-n01\n\#\$ -o $directory\/$jobName.out\n\#\$ -e $directory\/$jobName.err\n";
    print $outfile "\#\$ -m ae\n\#\$ -M yuan.chen\@duke.edu\n";
    print $outfile "cd $directory\n";
    print $outfile "echo Working directory is $directory\necho Running on host `hostname`\necho Time is `date`\n";
}



sub Usage(){
  print << "    Usage";

	Usage:  $0 <-a read1.fq.gz> <-b read2.fq.gz> <-l lib_size> <-j Job_Name>

	<options>
		-s1 <int> trimed length at 5' end of read1, default 0
		-e1 <int> trimed length at 3' end of read1, default 0
		-s2 <int> trimed length at 5' end of read2, default 0
		-e2 <int> trimed length at 3' end of read2, default 0
		-f  <int> input file type: fq:0 ,fq.gz:1, default 1
		-b  <int> filter reads with many low quality bases,set a cutoff, default 40
		-n  <int> filter reads with >X percent base is N,set a cutoff, default 10
		-q  <int> quality shift, default:64 (illumia 1.5)
		-c  <int> correction of sequencing errors, default 0
		-p  <int> CPU usage of the correction, default 4
		-k  <int> start of kmer frequence cutoff, default 2
		-e  <int> end of kmer frequence cutoff, default 2
		-d  <int> maximum error bases allowed, default 2
		-r  <int> remove temporary files, default 1

    Usage

	exit(0);
};
