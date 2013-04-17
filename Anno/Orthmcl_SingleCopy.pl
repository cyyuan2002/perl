#!/usr/bin/perl
use strict;

=head
This program is used to search the single copy Orthologs from the finial results of orthomcl. The input files are "groups.txt", which created by mcl
and the species file contains all the species that you want to use in the search, which should be consistant with the name used in orthomcl
sample:
coer
conc
cont
muco
=cut

my ($groupfile,$speciefile,$orthologNums)=@ARGV;
if(@ARGV<2){
    print "Usage:$0 <Group_File> <Species_File> <Orthlog_Num>\n";
    exit(0);
}

my %spes;
open(my $fh_spe,$speciefile) || die "Can't open file: $speciefile\n";
while(<$fh_spe>){
    chomp();
    $spes{$_}=1;
}
close $fh_spe;

if($orthologNums==""){
    $orthologNums=scalar(keys %spes);
}

my $specount=scalar(keys(%spes));
open(my $fh_group,$groupfile) || die "Can't open file: $groupfile\n";
while(<$fh_group>){
    chomp();
    my @orthologs=split(/\s/,$_);
    my %spefound;
    my $skip=0;
    my $spefoundnum;
    for(my $i=1;$i<@orthologs;$i++){
	my @info=split(/\|/,$orthologs[$i]);
	my $spe=$info[0];
	if(!exists($spefound{$spe})){
	    $spefound{$spe}=$orthologs[$i];
            $spefoundnum++;
	}
        else{
            $skip=1;
            last;
        }
    }
    next if($spefoundnum<$orthologNums || $skip==1);
    print "$orthologs[0]\t";
    foreach my $key (keys %spefound){
	print "$spefound{$key},";
    }
    print "\n";
}
close $fh_group;
exit(1);
