#!/usr/bin/perl
#Uncompleted
use strict;

my ($orthfile,$proteinFas,$outfolder)=@ARGV;
##----orthfile format----
#Group_4241:	gpro|52205,conc|77988,coer|49509,mver|MVEG_08461,bden|28647,rory|RO3G_06472,cont|g115.t1,muco|144714,pbla|129589,
#Group_4242:	gpro|95630,conc|2850,coer|80027,mver|MVEG_03816,bden|28658,rory|RO3G_06567,cont|g1237.t1,muco|26375,pbla|116798,
#Group_4243:	gpro|59696,conc|79896,coer|101489,mver|MVEG_03797,bden|28663,rory|RO3G_12350,cont|g3519.t1,muco|144611,pbla|123890,
#Group_4244:	gpro|135322,conc|71579,coer|79539,mver|MVEG_02761,bden|35737,rory|RO3G_06675,cont|g423.t1,muco|139407,pbla|133421,
#Group_4245:	gpro|499228,conc|87715,coer|16045,mver|MVEG_07087,bden|26848,rory|RO3G_14433,cont|g4905.t1,muco|148767,pbla|95960,


if(@ARGV<3){
    print "Usage:$0 <ortholog_info> <goodProteins.fasta> <output_folder>\n";
    exit(0);
}

my @spes;

my %fasseqs;
my $seqN;
my $seqS;
my $spe;
open (my $fh_proteinfas,$proteinFas) || die "Can't find fasta file: $proteinFas\n";
while(<$fh_proteinfas>){
    chomp();
    if(/^>(\S+)/){
        if($seqN ne ""){
            $fasseqs{$spe}->{$seqN}=$seqS;
        }
        my @seqName=split(/\|/,$1);
        $spe=$seqName[0];
        $seqN="$1";
        $seqS="";
    }
    else{
        $seqS.="$_";
    }
}
$fasseqs{$spe}->{$seqN}=$seqS;
close $fh_proteinfas;

if(!(-e $outfolder)){
    mkdir ($outfolder);
}

open (my $fh_orthfile,$orthfile) || die "Can't open file: $orthfile\n";
while(<$fh_orthfile>){
    chomp();
    my @lines=split(/\t/,$_);
    my $groupN=$lines[0];
    $groupN=~s/\://g;
    my @genes=split(/,/,$lines[1]);
    open(my $fh_out,">$outfolder\/$groupN.fas");
    foreach my $gid(@genes){
        my @ids=split(/\|/,$gid);
        my $speN=$ids[0];
        if(!exists($fasseqs{$speN}->{$gid})){
            print stderr "Can't find sequence: $gid\n";
        }
        else{
            my $seq=$fasseqs{$speN}->{$gid};
            print $fh_out ">$speN\n$seq\n";
        }
    }
    close $fh_out;
}
close $fh_orthfile;

exit(1);
