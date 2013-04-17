#!/usr/bin/perl
use strict;

my ($parIDFile,$parMatchField,$parIDField,$parTransFile,$parTransField,$parOutFile,$parIsfilter)=@ARGV;
if(@ARGV<6){
	print "Usage:$0 ID_File ID_MatchField_1 ID_Field Tran_File ID_MatchField_2 Output_File [filter_mode 0]\n";
	exit(1);
}

$parIsfilter ||=0;

if($parOutFile eq ""){
	$parOutFile="$parTransFile.trs";
}
my %idinfos;
open(my $fh_idfile,$parIDFile) || die "Can't open file: $parIDFile\n";
while(<$fh_idfile>){
	chomp();
	my @lines=split(/\t/,$_);
	my @ids=split(/\,/,$lines[$parMatchField]);
	foreach my $id(@ids){
		$idinfos{$id}=$lines[$parIDField];
	}
}
close $fh_idfile;

open(my $fh_transfile,$parTransFile) || die "Can't open file: $parTransFile\n";
open(my $fh_outfile, ">$parOutFile");
while(<$fh_transfile>){
	chomp();
        $_=~s/\"//g;
	my @lines=split(/\t/,$_);
	if(exists($idinfos{$lines[$parTransField]})){
		$lines[$parTransField]=$idinfos{$lines[$parTransField]};
		my $line=join("\t",@lines);
		print $fh_outfile "$line\n";
	}
	else{
                if($parIsfilter==0){
                    print $fh_outfile "$_\n";
                }
	}
}
close $fh_outfile;
close $fh_transfile;
exit(0)
