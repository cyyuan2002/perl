#!/usr/bin/perl
use strict;

my $parFasfile=shift;

open(my $fh_fasfile,$parFasfile) || die "Can't open file: $parFasfile\n";
my $seq;
my $seqcount;
my $lengthsum;
my $lastseqid;

while(<$fh_fasfile>){
	chomp();
	if(/^>(\S+)/){
		$seqcount++;
		$lengthsum+=length($seq);
		if($lastseqid ne ""){
			my $seqlen=length($seq);
			print "$lastseqid\t$seqlen\n";
		}
		$seq="";
		$lastseqid=$1;
	}
	else{
		$seq.=$_;
	}
}
close $fh_fasfile;
{
	my $seqlen=length($seq);
	print "$lastseqid\t$seqlen\n";
}
my $averlength=$lengthsum/$seqcount;
print "Total Sequences: $seqcount\n";
print "Total Sequences Length: $lengthsum\n";
print "Average Sequence Length: $averlength\n";