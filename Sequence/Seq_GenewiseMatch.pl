#!/usr/bin/perl
use strict;

my ($parAnnoFile,$parGWFile,$parOutfile,$parMatchPercent)=@ARGV;
if(@ARGV<4){
    print "Usage:$0 Annotation_file GeneWise_file Output_file Match_percent\n";
    exit(1);
}
if(!(-e $parAnnoFile)){
    print stderr "Can't find file: $parAnnoFile\n";
    exit(1);
}
if(!(-e $parGWFile)){
    print stderr "Can't find file: $parGWFile\n";
    exit(1);
}
open(my $annofile,$parAnnoFile);
my @geneAnno;
my %chromNum;
my $linecount=0;
my $lastchrom;
while(<$annofile>){
    chomp();
    my @lines=split(/\t/);
    my %info;
    if($lastchrom ne $lines[0]){ $chromNum{$lines[0]}=$linecount; }
    $info{"chrom"}=$lines[0];
    $info{"start"}=$lines[1];
    $info{"end"}=$lines[2];
    $info{"dir"}=$lines[3];
    $info{"block"}=$lines[4];
    $info{"sizes"}=$lines[5];
    $info{"starts"}=$lines[6];
    $info{"gid"}=$lines[8];
    push(@geneAnno,\%info);
    $linecount++;
}
close annofile;
my $annoNum=$linecount;
$linecount=0;


my %wiseinfo;
my @starts;
my @ends;
my @noinfos;
my $cdslength;
open(outfile,">$parOutfile");
my $parLostfile=$parOutfile.".los";
open (lostfile,">$parLostfile");
open(my $gwfile,$parGWFile);
while(<$gwfile>){
    chomp();
    my @lines=split(/\t/,$_);
    if(/match/){
	$wiseinfo{"chrom"}=$lines[0];
	$wiseinfo{"score"}=$lines[5];
	$wiseinfo{"dir"}=$lines[6];
	if($lines[6] eq "+"){
	    $wiseinfo{"start"}=$lines[3];
	    $wiseinfo{"end"}=$lines[4];
	}
	else{
	    $wiseinfo{"start"}=$lines[4];
	    $wiseinfo{"end"}=$lines[3];
	}
	my @geneid=split(/\-/,$lines[8]);
	$wiseinfo{"gene"}=$geneid[0];
    }
    elsif(/cds/){
	push(@starts,$lines[3]);
	push(@ends,$lines[4]);
	$cdslength=$cdslength+abs($lines[4]-$lines[3])+1;
	#my $cdsR=$lines[3]." ".$lines[4];
	#push(@cds,$cdsR);
    }
    elsif(/\/\//){
	#print "gene: $wiseinfo{'gene'}\n";
	if($wiseinfo{"dir"} eq "+"){
	    $wiseinfo{"starts"}=\@starts;
	    $wiseinfo{"ends"}=\@ends;
	}
	else{
	    my @revstarts=reverse(@starts);
	    my @revends=reverse(@ends);
	    $wiseinfo{"starts"}=\@ends;
	    $wiseinfo{"ends"}=\@starts;
	}
	$wiseinfo{"cdslength"}=$cdslength;
	my $check=&findRegion();
	#print "check:$check\n";
	if($check ne "NULL"){
	    my %annoinfo=%{$geneAnno[$check]};
    	    print outfile "$wiseinfo{'gene'}\t$annoinfo{'gid'}\t$wiseinfo{'chrom'}\t$wiseinfo{'dir'}\t";
	    my @sts=@{$wiseinfo{"starts"}};
	    my @es=@{$wiseinfo{"ends"}};
	    my @ss;
	    for(my $i=0;$i<@sts;$i++){
		my $size=$es[$i]-$sts[$i]+1;
		print outfile "$size,";
	    }
	    print outfile "\t";
	    for(my $i=0;$i<@sts;$i++){
		print outfile "$sts[$i],";
	    }
	    print outfile "\n";
	}
	else{
	    print lostfile "$wiseinfo{'gene'}\t$wiseinfo{'chrom'}\t$wiseinfo{'dir'}\t";
	    my @sts=@{$wiseinfo{"starts"}};
	    my @es=@{$wiseinfo{"ends"}};
	    my @ss;
	    for(my $i=0;$i<@sts;$i++){
		my $size=$es[$i]-$sts[$i]+1;
		print lostfile "$size,";
	    }
	    print lostfile "\t";
	    for(my $i=0;$i<@sts;$i++){
		print lostfile "$sts[$i],";
	    }
	    print lostfile "\n";
	}
	$linecount++;
	if(!($linecount%500)){
	    print "$linecount line finished\n";
	}
	$cdslength=0;
	%wiseinfo={};
	@starts=();
	@ends=();
    }
}
close $gwfile;
close outfile;
close lostfile;


sub findRegion{
    my $chrom=$wiseinfo{"chrom"};

    for(my $i=0;$i<$annoNum;$i++){
	my $annoinfo=$geneAnno[$i];
	#print "$$annoinfo{'chrom'}\t$wiseinfo{'chrom'}\n";
	if($$annoinfo{"chrom"} ne $wiseinfo{"chrom"}){next;}

	#if($$annoinfo{"chrom"} ne $wiseinfo{"chrom"}){return "NULL";}
	my $annostart=$$annoinfo{"start"};
	my $annoend=$$annoinfo{"end"};
	my $annodir=$$annoinfo{"dir"};
	my $wisestart=$wiseinfo{"start"};
	my $wiseend=$wiseinfo{"end"};
	my $wisedir=$wiseinfo{"dir"};
	my $annoid=$$annoinfo{"gid"};
	#if($annoid eq "gene_28028"){
	#    print "gene_28028\n";
	#}
	if($annodir ne $wisedir){ next;}
	if(($annostart<=$wisestart && $annoend>=$wisestart) || ($wisestart<=$annostart && $wiseend>=$annostart)){
	    my @annostarts=split(/\,/,$$annoinfo{"starts"});
	    my @annosizes=split(/\,/,$$annoinfo{"sizes"});
	    my @wisestarts=@{$wiseinfo{"starts"}};
	    my @wiseends=@{$wiseinfo{"ends"}};
	    my @annoends;
	    for(my $j=0;$j<@annostarts;$j++){
		my $end=$annostarts[$j]+$annosizes[$j]-1;
		push(@annoends,$end);
	    }
	    my $lengthmatch;
	    for(my $m=0;$m<@annostarts;$m++){
		for(my $n=0;$n<@wisestarts;$n++){
		    if($wisestarts[$n]<=$annostarts[$m] && $wiseends[$n]>=$annostarts[$m]){
			if($wiseends[$n]<=$annoends[$m]){
			    my $overlength=$wiseends[$n]-$annostarts[$m]+1;
			    $lengthmatch+=$overlength;
			}
			else{
			    my $overlength=$annoends[$m]-$annostarts[$m]+1;
			    $lengthmatch+=$overlength;
			}
			last;
		    }
		    elsif($annostarts[$m]<=$wisestarts[$n] && $annoends[$m]>=$wisestarts[$n]){
			if($annoends[$m]<=$wiseends[$n]){
			    my $overlength=$annoends[$m]-$wisestarts[$n]+1;
			    $lengthmatch+=$overlength;
			}
			else{
			    my $overlength=$wiseends[$n]-$wisestarts[$n]+1;
			    $lengthmatch+=$overlength;
			}
			last;
		    }
		}
	    }

	    if(($lengthmatch/$wiseinfo{"cdslength"})*100>=$parMatchPercent){
		return $i;
	    }
	    elsif($wisestart>$annoend){
		return "NULL";
	    }
	}

    }
    return "NULL";

}
