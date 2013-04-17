#!/usr/bin/perl
use strict;

my ($parAnnofile,$parGOfile)=@ARGV;
my %gidinfo;

open(my $fh_annofile,$parAnnofile) || die "Can't open file: $parAnnofile\n";
while(<$fh_annofile>){
	chomp();
	my @lines=split(/\t/,$_);
	my @gnames=split(/,/,$lines[7]);
	my $gid=$lines[8];
	foreach my $i(@gnames){
		$gidinfo{$i}=$gid;
	}
} 
close $fh_annofile;

open(my $fh_gofile,$parGOfile) || die "Can't open file: $parGOfile\n";
my %goinfo;
my $lastgid;
my @GOs;
while(<$fh_gofile>){
	chomp();
	my @lines=split(/\t/,$_);
	if($lastgid ne $lines[0]){
		if(!exists($gidinfo{$lastgid})){
			my @GOnums=@GOs;
			$goinfo{$lastgid}->{"GOs"}=\@GOnums;
			#print "ADD1:$lastgid\n";
		}
		else{
			if(!exists($goinfo{$gidinfo{$lastgid}}->{"GOs"})){
				my @GOnums=@GOs;
				$goinfo{$gidinfo{$lastgid}}->{"GOs"}=\@GOnums;
				#print "ADD2:$gidinfo{$lastgid}\n";
			}
			else{
				my @GOnums=(@GOs,@{$goinfo{$gidinfo{$lastgid}}->{"GOs"}});
				my %GOhash;
				$GOhash{$_}++ for @GOnums;
				@GOnums=keys %GOhash;
				$goinfo{$gidinfo{$lastgid}}->{"GOs"}=\@GOnums;
				#print "ADD3:$gidinfo{$lastgid}\n";
			}
		}
		
		if(!exists($gidinfo{$lines[0]})){
			if($lines[2] ne ""){
				$goinfo{$lines[0]}->{"des"}=$lines[2];
			}
		}
		else{
			if(!exists($goinfo{$gidinfo{$lines[0]}})){
				if($lines[2] ne ""){
					$goinfo{$gidinfo{$lines[0]}}->{"des"}=$lines[2];
				}
			}
		}
		$lastgid=$lines[0];
		@GOs=();
	}
	push (@GOs,$lines[1]);
}
close $fh_gofile;

foreach my $key(keys %goinfo){
	my @gos=@{$goinfo{$key}->{"GOs"}};
	if($goinfo{$key}->{"des"} ne ""){
		print "$key\t$gos[0]\t$goinfo{$key}->{'des'}\n";
	}
	else{
		print "$key\t$gos[0]\n";
	}
	for(my $i=1;$i<@gos;$i++){
		print "$key\t$gos[$i]\n"
	}
}