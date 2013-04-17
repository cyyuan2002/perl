#!/usr/bin/perl
use strict;

my ($parAnnofile,$parBlastregion,$parOutfile)=@ARGV;

if(@ARGV<3){
    print "Useage:$0 Anno_File Blastlocus_File Output_File\n";
    exit(1);
}

my $parCoverage=60;
my $iserror=0;
my %chromregion;
my $parFilterOutfile="$parOutfile.fout";

&checkfile($parAnnofile);
&checkfile($parBlastregion);
if($iserror==1){
    exit(1);
}

open(my $fh_annofile,$parAnnofile);
my $chrom;
my @generegions;

my $count;
while(<$fh_annofile>){
    chomp();
    my @lines=split(/\t/);
    if($chrom ne $lines[0]){
	$chrom=$lines[0];
	$count=0;
    }
    my %geneinfo;
    $geneinfo{'chrS'}=$lines[2];
    $geneinfo{'chrE'}=$lines[3];
    $geneinfo{'strand'}=$lines[4];
    $geneinfo{'exonnum'}=$lines[5];
    $geneinfo{'exons'}=$lines[6];
    $geneinfo{'refgene'}=$lines[7];
    $geneinfo{'geneid'}=$lines[8];
    $chromregion{$chrom}->[$count]=\%geneinfo;
    $count++;
}
#$chromregion{$chrom}=\@generegions;
close $fh_annofile;

open(my $fh_blastregion,$parBlastregion);
open(my $fh_output,">$parOutfile");
open(my $fh_filterout,">$parFilterOutfile");
while(<$fh_blastregion>){
    chomp();
    my @lines=split(/\t/,$_);
    if($lines[2]<$parCoverage){
	print $fh_filterout "$lines[0]\n";
	next;
    }
    my @annogenes=@{$chromregion{$lines[3]}};
    for(my $i=0;$i<@annogenes;$i++){
	my %geneinfo=%{$annogenes[$i]};
	if(($geneinfo{'chrS'}<=$lines[5] && $geneinfo{'chrE'}>=$lines[5] ) ||($lines[5]<=$geneinfo{'chrS'} && $lines[6] >=$geneinfo{'chrS'})){
		my @annoexons=&exonarray($geneinfo{'exons'},"+");
		my @geneexons=&exonarray($lines[9],$lines[4]);
		my $isoverlap=&checkisoform(\@geneexons,\@annoexons,$lines[1]);
		if($isoverlap==1){
		    if($geneinfo{'strand'} eq $lines[4]){
			print $fh_output "$lines[0]\t2\t$geneinfo{'refgene'}\t+\n";
		    }
		    else{
			print $fh_output "$lines[0]\t2\t$geneinfo{'refgene'}\t-\n";
		    }
		    last;
		}
		elsif($isoverlap==2){
		    if($geneinfo{'strand'} eq $lines[4]){
			print $fh_output "$lines[0]\t1\t$geneinfo{'refgene'}\t+\n";
		    }
		    else{
			print $fh_output "$lines[0]\t1\t$geneinfo{'refgene'}\t-\n";
		    }
		    last;
		}
		else{
		    print $fh_output "$lines[0]\t0\n";
		    last;
		}
	    }
	#}
	elsif($geneinfo{'chrS'}>$lines[6]){
	    print $fh_output "$lines[0]\t0\n";
	    last;
	}
    }

}

sub exonarray{
    my ($exoninfo,$dir)=@_;
    my @exon=split(/,/,$exoninfo);
    my @exons;
    for(my $i=0;$i<@exon;$i++){
	my @info=split(/:/,$exon[$i]);
	my %infor=('ts'=>$info[0],'te'=>$info[1]);
	push(@exons,\%infor);
    }
    if($dir eq "-"){
	@exons=reverse(@exons);
    }

    @exons=&exonmerge(@exons);
    return @exons;
}

sub exonmerge{
    my @exons=@_;
    for(my $i=0;$i<@exons-1;$i++){
	if($exons[$i]->{'te'}>=$exons[$i+1]->{'ts'}){
	    $exons[$i]->{'te'}=$exons[$i+1]->{'te'};
	    splice(@exons,$i+1,1);
	    $i--;
	}
    }
    return @exons;
}

sub checkisoform{  ##check the isoform of different blat hit
    my ($tinfo,$sinfo,$qsize)=@_;
    my @texons=@{$tinfo};
    my @sexons=@{$sinfo};

    my $exonmatch=0;
    my $lengthmatch=0;

    #my $exonnum=0;
    #my $querylength=0;

    my $texoncount=0;
    my $sexoncount=0;

    while($texoncount<@texons && $sexoncount<@sexons){
	my %texonregion=%{$texons[$texoncount]};
	my %sexonregion=%{$sexons[$sexoncount]};
	my $sexonlength=$sexonregion{'te'}-$sexonregion{'ts'}+1;
	if($sexonlength<1){ ##
	    $sexoncount++;
	    next;
	}

	my $texonlength=$texonregion{'te'}-$texonregion{'ts'}+1;
	if($texonlength<1){
	    $texoncount++;
	    next;
	}
	#my @texonregion=split(/\:/,$texons[$texoncount]);
	#my @sexonregion=split(/\:/,$sexons[$sexoncount]);
	if($texonregion{'ts'}<=$sexonregion{'ts'} && $texonregion{'te'}>=$sexonregion{'ts'}){
	    if($texonregion{'te'}<=$sexonregion{'te'}){
	        my $overlength=$texonregion{'te'}-$sexonregion{'ts'}+1;
		$lengthmatch+=$overlength;
		if($overlength/$sexonlength>0.9 || $overlength/$texonlength>0.9){
		    $exonmatch++;
		}
		$texoncount++;
	    }
	    else{
		my $overlength=$sexonregion{'te'}-$sexonregion{'ts'}+1;
		$lengthmatch+=$overlength;
		if($overlength/$sexonlength>0.9 || $overlength/$texonlength>0.9){
		    $exonmatch++;
		}
		$sexoncount++;
	    }
	}
	elsif($sexonregion{'ts'}<=$texonregion{'ts'} && $sexonregion{'te'}>=$texonregion{'ts'}){
	    if($texonregion{'te'}<=$sexonregion{'te'}){
		my $overlength=$texonregion{'te'}-$texonregion{'ts'}+1;
		$lengthmatch+=$overlength;
		if($overlength/$sexonlength>0.9 || $overlength/$texonlength>0.9){
		    $exonmatch++;
		}
		$texoncount++;
	    }
	    else{
		my $overlength=$sexonregion{'te'}-$texonregion{'ts'}+1;
		$lengthmatch+=$overlength;
		if($overlength/$sexonlength>0.9 || $overlength/$texonlength>0.9){
		    $exonmatch++;
		}
		$sexoncount++;
	    }
	}
	else{
	    if($texonregion{'ts'}>$sexonregion{'te'}){
		$sexoncount++;
	    }
	    elsif($sexonregion{'ts'}>$texonregion{'te'}){
		$texoncount++;
	    }
	}
    }

    my $qmp=$lengthmatch/$qsize;

    if($qmp > 0.6 || $exonmatch > 3){
	return 1;
    }
    elsif($qmp>0.3 || $exonmatch > 1){
	return 2;
    }
    else{
	return 0;
    }
}


sub checkfile{
    my $filename=shift;
    if(!(-e $filename)){
	print stderr "Can't find file: $filename\n";
	$iserror=1;
    }
}
