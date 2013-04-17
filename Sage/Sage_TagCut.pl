#!/usr/bin/perl
#this program is used to find CATG sites in the target database, all the sites are searched in bidirectory,the Redundent sequences are removed from the tag data
#output format:

#seqname	seqlength	tagseq	strand	startsite	redundant
#fdsafdsa	63	AAAAAAAAAAAAACCGA	+	22	0
#fdsafdsa	63	TTTTTTTTTTTTTCATG	-	1	1


use strict;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"o:s","i:s","l:s","r:s","help");

my $version="1.1 Alvin Chen 2010-12-27";
if(!defined $opts{i}){
    &Usage();
}

my $inputfile=$opts{i};
my $outputfile=(defined $opts{o}) ? $opts{o} :"$inputfile.tags";
my $taglength=(defined $opts{l}) ? $opts{l} : 17;
my $redunt=(defined $opts{r}) ? $opts{r} : 1;

open(input,$inputfile) || die "Can't open file: $inputfile\n";
my $seq;
my $seqname;
if($outputfile eq ""){
	$outputfile=$inputfile.".tags";
}

if($taglength eq ""){
    $taglength=17;
}

my %taginfo;
my @tagorder;
my %redtag;
my %redgene;
my %redgenerm;
my %redmulti;

while(<input>){
    chomp();
    if(/^>(\S+)/){
	if($seq ne ""){
	    &sitesearch($seqname,$seq);
	    $seqname=$1;
	    $seq="";
	}
	else{
	    $seqname=$1;
	    $seq="";
	}
    }
    else{
	$seq=$seq.$_;
    }
}
&sitesearch($seqname,$seq);
close input;


#my %redgeneinfo;

open (output,">$outputfile");

for(my $i=0;$i<@tagorder;$i++){
    if($taginfo{$tagorder[$i]} ne ""){
	print output "$taginfo{$tagorder[$i]}0\n";
    }
}
if($redunt>0){
    foreach my $key(keys %redgene){
	foreach my $keyin(keys %{$redgene{$key}}){
	    if(!exists($redgenerm{$keyin})){ #"$seqN\t$seqlength\t$tagrev\t-\t$revS\t";
		print output "$key\t$redgene{$key}->{$keyin}->{'length'}\t$keyin\t$redgene{$key}->{$keyin}->{'dir'}\t$redgene{$key}->{$keyin}->{'pos'}\t1\n";
	    }
	}
    }
    if($redunt==2){
	foreach my $key(keys %redgenerm){
	    foreach my $keyin(keys %{$redgenerm{$key}}){
		print output "$keyin\$";
	    }
	    print output "\tv\t$key\tx\tm\t2\n";
	}
    }
}
#for(my $i=0;$i<@redtags;$i++){
    #print output "$redtags[$i]1\n";
#}

close output;

sub sitesearch(){
    my ($seqN,$seqS)=@_;
    my $seqSU=uc($seqS);
    my $seqlength=length($seqSU);
    while($seqSU=~/CATG/g){
	my $currentpos=pos($seqSU);
	my $tagsequence=substr($seqSU,$currentpos,$taglength);
	if(!($tagsequence=~/NN/)){
	    my $accpos=$currentpos+1;
	    if(length($tagsequence)==$taglength){
		if(!exists($taginfo{$tagsequence})){
		    my $taginfo="$seqN\t$seqlength\t$tagsequence\t+\t$accpos\t";
		    $taginfo{$tagsequence}=$taginfo;
		    push (@tagorder,$tagsequence);
		}
		else{
		    my @infors=split(/\t/,$taginfo{$tagsequence});
		    $taginfo{$tagsequence}="";
		    my $dir;
		    if($redunt>0){
			if($infors[3] ne "+"){
			    $dir="x";
			}
			else{
			    $dir="+";
			}
			&checktag($seqN,$seqlength,$tagsequence,$dir,"m");
		    }
		}
	    }
	}
	my $revstart=$currentpos-4-$taglength;
	if($revstart>=0){
	    my $tagrev=substr($seqSU,$revstart,$taglength);
	    if(!($tagrev=~/N/)){
		if(length($tagrev)==$taglength){
		    $tagrev=reverse($tagrev);
		    $tagrev=~ tr/[ATGC]/[TACG]/;
		    my $revS=$revstart+1;
		    if(!exists($taginfo{$tagrev})){
			my $taginfo="$seqN\t$seqlength\t$tagrev\t-\t$revS\t";
			$taginfo{$tagrev}=$taginfo;
			push (@tagorder,$tagrev);
		    }
		    else{
			my @infors=split(/\t/,$taginfo{$tagrev});
			$taginfo{$tagrev}="";
			if($redunt>0){
			    my $dir;
			    if($redunt>0){
				if($infors[3] ne "-"){
				    $dir="x";
				}
				else{
				    $dir="-";
				}
			    }
			    &checktag($seqN,$seqlength,$tagrev,$dir,"m");
			    #if(!exists($redtag{$tagrev})){
			    #    my $dir="$infors[3]|-";
			    #    my $tpos="$infors[4]|$revS";
			    #    &checktag($seqN,$seqlength,$tagrev,$dir,$tpos);
			    #}
			    #else{
			    #    &checktag($seqN,$seqlength,$tagrev,"-",$revS);
			    #}
			}
		    }
		}
	    }
	}
    }
}

sub checktag{
    my ($Tgene,$TgeneL,$Tseq,$Tdir,$Tpos)=@_;
    if(!exists($redtag{$Tseq})){
	$redtag{$Tseq}=$Tgene;
	$redgene{$Tgene}->{$Tseq}->{'length'}=$TgeneL;
	$redgene{$Tgene}->{$Tseq}->{'dir'}=$Tdir;
	$redgene{$Tgene}->{$Tseq}->{'pos'}=$Tpos;
    }
    else{
	if(!exists($redgenerm{$Tseq})){
	    if($redtag{$Tseq} eq $Tgene){
		#$redgene{$Tgene}->{$Tseq}->{'dir'}="$redgene{$Tgene}->{$Tseq}->{'dir'}|$Tdir";
		#$redgene{$Tgene}->{$Tseq}->{'pos'}="$redgene{$Tgene}->{$Tseq}->{'pos'}|$Tpos";
		if($redgene{$Tgene}->{$Tseq}->{'dir'} ne $Tdir){
		    $redgene{$Tgene}->{$Tseq}->{'dir'}="x";
		}
	    }
	    else{
		$redgenerm{$Tseq}->{$Tgene}=1;
		$redgenerm{$Tseq}->{$redtag{$Tseq}}=1;
	    }
	}
	else{
	    if(!exists($redgenerm{$Tseq}->{$Tgene})){
		$redgenerm{$Tseq}->{$Tgene}=1;
	    }
	}
    }
}



sub Usage() #help subprogram
{
    print << "    Usage";
    
	Usage: $0 version $version <options>

		-i     Name of input fasta file , must be given
					
		-o     File for output the tags (default: xxxxx.tags) 
		
		-l     Length of the tags (default: 17)
	
		-r     Output the redundent tags (0 (No redundent tags), 1 (Including gene unique redundent) or 2 (All) default: 1)

		-help  Show help

    Usage

	exit(0);
};		

