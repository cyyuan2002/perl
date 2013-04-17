#!/usr/bin/perl
use strict;

my ($parFilein)=shift;
open(my $fh_parfilein,$parFilein) || die "Can't open file: $parFilein";
my $lastID;
<$fh_parfilein>;
my $first=0;
my $isoutput=0;
while(<$fh_parfilein>){
	chomp();
	my @lines=split(/\t/,$_);
	if($lastID ne $lines[0]){
		print "\n" if($lastID ne "" && $first==1);
		$lastID=$lines[0];
		#print "$lines[0]\t";
		$first=0;
		if($lines[1] ne ""){
			if($lines[1]=~/EC/){
				next;
			}
			else{
				print "$lastID\t$lines[1]";
				$first=1;
			}
		}
	}
	else{
		if($lines[1] ne ""){
			if($first==0){
				if($lines[1]=~/EC/){
					next;
				}
				else{
					print "$lastID\t$lines[1]";
					$first=1;
				}
			}
			else{
				next if($lines[1]=~/EC/);
				print ", $lines[1]";
			}
		}
	}
}
close $fh_parfilein;
exit(0);