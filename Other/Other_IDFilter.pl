#!/usr/bin/perl
use strict;
use Getopt::Long;

my $version="1.0 Alvin Chen 2011-7-7";

my %opts;
GetOptions(\%opts,"a=s","af:i","b=s","bf:i","bi:s","m:i","help");
if((!defined $opts{a})||(!defined $opts{b})){
    &Usage();
}

my $parFileA=$opts{a};
my $parFileB=$opts{b};
my $parFileAField=(defined $opts{af}) ? $opts{af} : 0;
my $parFileBField=(defined $opts{bf}) ? $opts{bf} : 0;
my $parMode=(defined $opts{m}) ? $opts{m} : 0;
my $parFileBAdd=(defined $opts{bi}) ? $opts{bi} : "";
my %idinfo;


my $iserror=0;
if($parMode>0 && $parFileBAdd eq ""){
	print stderr "Error: No parameters for '-bi' when -m > 0\n";
	$iserror=1;
}

Filecheck($parFileA,0);
Filecheck($parFileB,0);

exit(1) if ($iserror==1);
open(my $fh_fileb,$parFileB);
my @FieldsAdd;
if($parMode>0){
	@FieldsAdd=split(/,/,$parFileBAdd);
}
while(<$fh_fileb>)
{
	chomp();
	my @lines=split(/\t/,$_);
	if($parMode==0){
		#print "$lines[$parFileBField]\n";
		$idinfo{$lines[$parFileBField]}=1;
	}
	else{
		my $Addinfo=$lines[$FieldsAdd[0]];
		if(scalar(@FieldsAdd>1)){
			for(my $i=1;$i<@FieldsAdd;$i++){
				$Addinfo.="\t$lines[$FieldsAdd[$i]]";
			}
		}
		$idinfo{$lines[$parFileBField]}=$Addinfo;
	}
}
close $fh_fileb;

open(my $fh_filea,$parFileA);
while(<$fh_filea>){
	chomp();
	my @lines=split(/\t/,$_);
	if(exists($idinfo{$lines[$parFileAField]})){
		if($parMode==0){
			print "$_\n";
		}
		elsif($parMode==1){
			for(my $i=0;$i<@lines;$i++){
				if($i==$parFileAField){
					print "$idinfo{$lines[$i]}\t";
				}
				else{
					print "$lines[$i]\t";
				}
			}
			print "\n";
		}
		else{
			print "$_\t$idinfo{$lines[$parFileAField]}\n";
		}
	}
}
close $fh_filea;

exit(0);

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

	This program is used to filtered information in fileA by IDs in fileB

	Usage:  $0 (version $version)

	<options>
		-a     Filename of FileA
		-af    Field of ID in FileA (default:0)
		-b     Filename of FileB
		-bf    Field of ID in FileB (default:0)
		-bi    Other fields need to be added (exam:1,2,3), only works when \"-r 1/2\"
		-m     Modes of matched IDs (default 0:keep info in FileA; 1:change ID info by -bi fields; 2:Add info at back)

    Usage

	exit(0);
};
