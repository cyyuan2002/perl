#!/usr/bin/perl
##This program is used to submit jobs to servers

use strict;
my ($shellfile)=shift;
open(shellfile,$shellfile);
my @shells=<shellfile>;
my %jobids;

for(my $i=0;$i<@shells;$i++){
    my $jobname=`bsub -q Q_Serial $shells[$i]`;
#my $jobname=`bsub $shells[$i]`;
    $jobname=~/^Job\s\<(\d+)\>/;
    $jobids{$1}=1;
}

my $jobcount=scalar(@shells);

while($jobcount>0){
    sleep 600;
    my $jobstat=`bjobs`;
    if($jobstat eq "No unfinished job found\n"){
	$jobcount=0;
    }
    else{
	my @jobstats=split(/\n/,$jobstat);
	if(@jobstats<2){
	    $jobcount=1;
	}
	else{
	    $jobcount=0;
	    for(my $i=1;$i<@jobstats;$i++){
		my @jobinfo=split(/\s+/,$jobstats[$i]);
		if(exists($jobids{$jobinfo[0]})){
		    $jobcount++;
		}
	    }
	}
    }
}
