#!/usr/bin/perl
#===============================================================================
#
#         FILE: Other_PBSCreator.pl
#
#        USAGE:
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
GetOptions(\%opts,"i:s","a:s","n:s","r:s","e:s","p:s","q:i","help");

if((!defined $opts{i})||(!defined $opts{n})){
	Usage();
}

my $fileName;
my $pbsName;
my $ppNumber;
my $pbsFile;
my $directory;
my $emailnotify;
my $querymode;
my $addition;

$fileName=$opts{i};
$pbsName=$opts{n};
$addition=$opts{a};
$directory=(defined $opts{r})?$opts{r}: getcwd();
$ppNumber=(defined $opts{p}) ? $opts{p} : 1;
$emailnotify=(defined $opts{e}) ? $opts{e} : 1;
$querymode=(defined $opts{q}) ? $opts{q} : 0;

my @qsub;
my $linecount=0;
open (my $filein,$fileName) || die "Can't open file: $fileName\n";
if($querymode==0){
	while(<$filein>){
		my $outfile=$pbsName.".Job".$linecount.".PBS";
		my $pbsoutFile=(defined $opts{o}) ? $opts{o} : "$outfile.out";
		my $pbserrFile=(defined $opts{e}) ? $opts{e} : "$outfile.err";
		push (@qsub,$outfile);
		open (my $out,">$outfile");
		print $out "\#!/bin/bash\n\#\$ -v PATH\n\#\$ -v LD_LIBRARY_PATH=\"\$LD_LIBRARY_PATH:\/opt\/apps\/lib64\"\n";
		print $out "$addition\n" if($addition ne "");
		print $out "\#\$ -N $pbsName\_job$linecount\n\#\$ -pe threaded $ppNumber\n\#\$ -q \*\@mmrl-n01\n\#\$ -o $directory\/$pbsoutFile\n\#\$ -e $directory\/$pbserrFile\n";
		if($emailnotify==1){
			print $out "\#\$ -m ae\n\#\$ -M yuan.chen\@duke.edu\n";
		}
		print $out "cd $directory\n";
		print $out "echo Working directory is $directory\necho Running on host `hostname`\necho Time is `date`\n";
		print $out "$_";
		close $out;
		$linecount++;
		`qsub -l highprio $outfile`;
		print stderr "Job $outfile submitted\n";
	}

}
else{
	my $outfile=$pbsName."PBS";
	my $pbsoutFile=(defined $opts{o}) ? $opts{o} : "$outfile.out";
	my $pbserrFile=(defined $opts{e}) ? $opts{e} : "$outfile.err";
	open (my $out,">$outfile");
	print $out "\#!/bin/bash\n\#\$ -N $pbsName\_job$linecount\n\#\$ -pe threaded $ppNumber\n\#\$ -q \*\@mmrl-n01\n\#\$ -o $directory\/$pbsoutFile\n\#\$ -e $directory\/$pbserrFile\n";
	if($emailnotify==1){
		print $out "\#\$ -m ae\n\#\$ -M yuan.chen\@duke.edu\n";
	}
	print $out "cd $directory\n";
	print $out "echo Working directory is $directory\necho Running on host `hostname`\necho Time is `date`\n";
	while(<$filein>){
		chomp();
		print $out "echo Running $_\necho Time is `date`\n";
		print $out "$_\n";
		print $out "echo $_ finished\n\n";
	}
	close $out;
	`qsub -l highprio $out`;
}


sub Usage(){
	print << "    Usage";

	Version	  : 1.2
	Date	  : 2012-4-12
	Author    : Yuan,Chen (MMRL,Duke University Medical Center)
        Email     : yuan.chen\@duke.edu

	Usage: $0 <options>

		-i     input of command file , must be given (String)

		-n     name of the PBS job, must be given (String)

		-a     additional parameters of PBS

		-r     the directory of the PBS job , default: current dir(String)

		-p     the process number of PBS work, default value is 1

                -e     email notification when job finished, default 1

		-q     submit jobs in queue, default 0

		-help  Show help

    Usage

	exit(0);
}
