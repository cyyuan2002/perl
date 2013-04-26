#!/usr/bin/perl
use strict;
use Getopt::Long;
use Cwd;
my %opts;
GetOptions(\%opts,"c:s","n:s","r:s","e:s","m:s","p:s","a:s","help");

if((!defined $opts{c})||(!defined $opts{n})){
	Usage();
}

my $fileName;
my $pbsName;
my $ppNumber;
my $pbsFile;
my $directory;
my $emailnotify;
my $command;
my $nodeName;
my $addition;

$command=$opts{c};
$pbsName=$opts{n};
$addition=$opts{a};
$directory=(defined $opts{r})?$opts{r}: getcwd();
$ppNumber=(defined $opts{p}) ? $opts{p} : 1;
$emailnotify=(defined $opts{e}) ? $opts{e} : 1;
$nodeName=$opts{m};

my @qsub;
open(my $outfile,">temp.PBS");
my $pbsoutFile=(defined $opts{o}) ? $opts{o} : "$pbsName.out";
my $pbserrFile=(defined $opts{e}) ? $opts{e} : "$pbsName.err";
print $outfile "\#!/bin/bash\n\#\$ -v PATH\n\#\$ -v LD_LIBRARY_PATH=\"\$LD_LIBRARY_PATH:\/opt\/apps\/lib64\"\n\#\$ -N $pbsName\n\#\$ -pe threaded $ppNumber\n";
print $outfile "$addition\n" if ($addition ne "");
if($nodeName ne ""){
	if ($nodeName eq 'm') {
		print $outfile "\#\$ -q \*\@mmrl-n01\n";
	}
	else{
		print $outfile "\#\$ -q \*\@$nodeName\n";
	}
}
print $outfile "\#\$ -o $directory\/$pbsoutFile\n\#\$ -e $directory\/$pbserrFile\n";
if($emailnotify==1){
	print $outfile "\#\$ -m ae\n\#\$ -M yuan.chen\@duke.edu\n";
}
print $outfile "cd $directory\n";
print $outfile "echo $command\necho Working directory is $directory\necho Running on host `hostname`\necho Time is `date`\n";
print $outfile "$command\n";
close $outfile;
`qsub -l highprio temp.PBS`;
unlink "temp.PBS";
print "Job submitted!\n";
exit(1);

sub Usage(){
	print << "    Usage";

	Version	  : 1.0
	Date	  : 2012-5-5
	Author    : Yuan,Chen (MMRL,Duke University Medical Center)
        Email     : yuan.chen\@duke.edu

	Usage: $0 <options>

		-c     command, must be given (String)

		-n     name of the PBS job, must be given (String)

		-r     the directory of the PBS job , default: current dir(String)

		-a     additional parameters for PBS

		-p     the process number of PBS work, default value is 1

		-e     email notification when job finished, default 1

		-m     job nodes, default: mmrl

		-help  Show help

    Usage

	exit(0);
}
