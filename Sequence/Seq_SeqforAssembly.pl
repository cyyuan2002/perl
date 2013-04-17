#!/usr/bin/perl
#this program is used to extract sequences from anno sequence for assembly

use strict;

my ($parAnnofile,$parFasfile,$parOutfile)=@ARGV;
if(@ARGV<3){
    print "Usage:$0 Anno_file Fasta_file Output_file\n";
    exit(1);
}
my $iserror=0;

&Filecheck($parAnnofile,0);
&Filecheck($parFasfile,0);

my %seqs;

if($iserror==1){
    exit(1);
}

open(my $fh_fasfile,$parFasfile);
my $seqName;
my $seq;
while(<$fh_fasfile>){
    chomp();
    if(/^>(\S+)/){
		if($seq ne ""){
		    $seqs{$seqName}=$seq;
		}
		$seqName=$1;
		$seq="";
    }
    else{
		$seq=$seq.$_;
    }
}
$seqs{$seqName}=$seq;
close $fh_fasfile;

open (my $fh_annofile,"$parAnnofile");
open (my $fh_outfile, ">$parOutfile");
my %seqids;

while(<$fh_annofile>){
	chomp();
	my @lines=split(/\t/,$_);
	print $fh_outfile "\#$lines[2]\t$lines[1]\n";
	my @fasseqs=&getseq($lines[0]);
	my $fastaseqs;
	for(my $i=0;$i<@fasseqs;$i++){
	    my $seqN=$fasseqs[$i]->{"id"};
	    my $seq=$fasseqs[$i]->{"seq"};
	    $fastaseqs.=">$seqN\n$seq\n";
	}
	print $fh_outfile "$fastaseqs";
}
close $fh_outfile;
exit(0);

sub getseq{
    my $seqids=shift;
    my @seqid=split(/,/,$seqids);
    my @seqfastas;
    for(my $i=0;$i<@seqid;$i++){
    	if($seqid[$i]=~/\_/){
			my @seqN=split(/\_/,$seqid[$i]);
			$seqid[$i]=$seqN[0];
		}
		if(exists($seqs{$seqid[$i]})){
		    #$seq=&checkfas($seqid[$i],$seqs{$seqid[$i]});
		    my %seqinfo=("id"=>$seqid[$i],"seq"=>$seqs{$seqid[$i]});
		    push(@seqfastas,\%seqinfo);
		}
		else{
		    print "Error: Cant' find sequence $seqid[$i]\n";
		    #exit(1);
		}
    }
    return (@seqfastas);
}

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
	    print stderr "File $Filename already exists\n";
	    $iserror=1;
	}
    }
}

