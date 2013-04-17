#!/usr/bin/perl
##mapping the tags to the cutted file, defined the sequence direction by the max tag mapping direction 4X4

use strict;

my ($taglib,$tagcount,$mismatch) =@ARGV;
if(@ARGV<2){
    print stderr "Usage:$0 tag_library tag_count mismatch(0)\n";
    exit;
}

if($mismatch==""){
    $mismatch=0;
}
elsif($mismatch>3){
    print stderr "Mismatch must less than 4!\n";
    exit;
}

open(taglib,$taglib) || die "Can't open file:$taglib\n";
my %taggene;
my %tagdir;
my %tagpos;

##read taglib files;
while(<taglib>){
    chomp();
    my @lines=split(/\t/,$_);
    $taggene{$lines[2]}=$lines[0];
    $tagdir{$lines[2]}=$lines[3];
    $tagpos{$lines[2]}=$lines[4];
}
close taglib;

##create index table;
my @taglibs=sort(keys %taggene);
my @bases=("A","C","G","T");
my $indexcount=0;
my %indexbases;
my %indexbasesE;
my $lastbases="";

if($mismatch>0){
    for(my $a=0;$a<4;$a++){
	for(my $b=0;$b<4;$b++){
	    for(my $c=0;$c<4;$c++){
		for(my $d=0;$d<4;$d++){
		    my $indexstr="$bases[$a]$bases[$b]$bases[$c]$bases[$d]";
		    for(my $i=$indexcount;$i<@taglibs;$i++){
			my $top4=substr($taglibs[$i],0,4);
			if($indexstr eq $top4){
			    if($lastbases ne ""){
				$indexbasesE{$lastbases}=$i;
			    }
			    $indexbases{$indexstr}=$i;
			    $lastbases=$indexstr;
			    $indexcount=$i;
			    last;
			}
		    }
		}
	    }
	}
    }
    $indexbasesE{$lastbases}=scalar(@taglibs);
}

open(tagcount,$tagcount) || die "Can't open file:$tagcount\n";
my $outfile=$tagcount.".map";
open(outfile,">$outfile");

if($mismatch==0){
    while(<tagcount>){
        chomp();
	my @lines=split(/\t/,$_);
	if(exists($taggene{$lines[0]})){
	    print outfile "$lines[0]\t$taggene{$lines[0]}\t$tagdir{$lines[0]}\t0\t$tagpos{$lines[0]}\n";
	}
    }
}
else{
    while(<tagcount>){
	chomp();
	my @lines=split(/\t/,$_);
	if(exists($taggene{$lines[0]})){
	    print outfile "$lines[0]\t$taggene{$lines[0]}\t$tagdir{$lines[0]}\t0\t$tagpos{$lines[0]}\n";
	}
	else{
	    my $tagseq=$lines[0];
	    my $top4s=substr($tagseq,0,4);
	    my @top4=split("",$top4s);
	    foreach my $basekey (keys %indexbases){
		my @basekeyL=split("",$basekey);
		my $miscount=0;
		for(my $i=0;$i<4;$i++){
		    if($top4[$i] ne $basekeyL[$i]){
			$miscount++;
		    }
		}
		if($miscount<=$mismatch){
		    &tagmatch($tagseq,$basekey);
		}
	    }
	}
    }
}


sub tagmatch{
    my ($tag,$basekey)=@_;
    my $indexstart=$indexbases{$basekey};
    my $indexend=$indexbasesE{$basekey};
    my @tags=split("",$tag);
    my %tagmatch;  ##save multimatch
    my %tagmatchscore; ## save multimatch score
    my @scores;
    
    for(my $i=$indexstart;$i<$indexend;$i++){
	my $libtag=$taglibs[$i];
	my @libtags=split("",$libtag);
	my $miscount=0;
	for(my $i=0;$i<@libtags;$i++){
	    if($tags[$i] ne $libtags[$i]){
		my $misscore=&IUPACCheck($tags[$i],$libtags[$i]);
		$miscount=$miscount+$misscore;
	    }
	}
	if($miscount<=$mismatch){
	    $tagmatch{$libtag}="$tag\t$taggene{$libtag}\t$tagdir{$libtag}\t$miscount\t$tagpos{$libtag}";
	    $tagmatchscore{$libtag}=$miscount;
	    push(@scores,$miscount);
	}
    }
    
    my @sortscores=sort {$a <=> $b} @scores;
    for my $key (keys %tagmatch){
	if($tagmatchscore{$key}==$scores[0]){
	    print outfile "$tagmatch{$key}\n";
	}
    }
}

sub IUPACCheck{
    my ($baseA,$baseB)=@_;
    ##$baseA tag;
    if($baseB=~/[ATCG]/){
	return 1;
    }
    else{
	my $baseBX=&Trans($baseB);
	if($baseA=~/$baseBX/){
	    return 0; 
	}
	else{
	    return 1;
	}
    }
    
}

sub Trans{
    my $base=shift;
    if($base eq "M"){ return "AC";}
    elsif($base eq "R"){ return "AG";}
    elsif($base eq "W"){ return "AT";}
    elsif($base eq "S"){ return "CG";}
    elsif($base eq "Y"){ return "CT";}
    elsif($base eq "K"){ return "GT";}
    elsif($base eq "V"){ return "ACG";}
    elsif($base eq "H"){ return "ACT";}
    elsif($base eq "D"){ return "AGT";}
    elsif($base eq "B"){ return "CGT";}
    elsif($base eq "N"){ return "ATCG";}
}