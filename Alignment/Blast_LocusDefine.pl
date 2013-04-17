#!/usr/bin/perl
use strict;
use Bio::SearchIO;

my ($FblastN,$FTblastX,$FRegion,$FOut)=@ARGV;
if (@ARGV <4){
    print "Usage:$0 BlastN_File TBlastX_File Region_File Out_File\n";
    exit(1);
}
my $iserror=0;
&checkfile($FblastN);
&checkfile($FTblastX);
&checkfile($FRegion);
if ($iserror==1){
    exit(1);
}

my %reginfors;
my $parMinLength=30;
my $parMinIdentity=0.5;
my $parRdnOverpercent=0.90;
my $parFlankingLength=10000;
my %genemapping;

open(my $regionfile,"$FRegion");
while(<$regionfile>){
    chomp();
    my @lines=split(/\t/,$_);
    my %infor;
    $infor{'chrom'}=$lines[0];
    $infor{'length'}=$lines[1];
    $infor{'start'}=$lines[2];
    $infor{'end'}=$lines[3];
    $infor{'strand'}=$lines[4];
    $reginfors{$lines[9]}=\%infor;
}
close $regionfile;

my $searchBN=Bio::SearchIO->new(-format => 'blast', -file=>$FblastN);
my $searchTBX=Bio::SearchIO->new(-format => 'blast', -file => $FTblastX);

my $qName;
while(my $resBN=$searchBN->next_result()){
    #my $resTBN=$searchTBN->next_result();
    #only record the first hit, length coverage 50bp, identity 50%, the overlap region according to high identity region.
    #overlap: 1.sort by region first, blastn  2. cut overlap region, 3. use high identity region, 4. same identity use first.

    $qName=$resBN->query_name();
    #print "blastn qName:$qName";
    if($resBN->hits==0){
	next;
    }

    my $blastnhits=$resBN->hits;
    #print "blastn hits:$blastnhits\t";

    my $hitBN=$resBN->next_hit();
    my $sName=$hitBN->name();

    my @exons;
    my @hspBNs;
    my $hspcount=0;

    ##filter hsp by parameters
    my $strand;
    while(my $hspBN=$hitBN->next_hsp()){
	if($strand eq ""){
	    $strand=$hspBN->strand('hit');
	}
	if($hspBN->length('total')<$parMinLength || $hspBN->frac_identical('total')<$parMinIdentity || $hspBN->strand('hit') ne $strand){
	    next;
	}
	if($hspBN->strand('hit')<0){
	    my $tempS=$hspBN->{'HIT_START'};
	    $hspBN->{'HIT_START'}=$hspBN->{'HIT_END'};
	    $hspBN->{'HIT_END'}=$tempS;
	}
	$hspBNs[$hspcount]=$hspBN;
	$hspcount++;
    }
    my $genestrand;
    if($strand eq "1"){$genestrand="+";}
    else{$genestrand="-";}

    if(scalar(@hspBNs)==0){
	next;
    }

    my @sorted_hspBNs;
    if(scalar(@hspBNs)<2){
	@sorted_hspBNs=@hspBNs;
    }
    else{
	@sorted_hspBNs=hspsort(@hspBNs);
    }

    my @refExon=overlapcheck($sName,$genestrand,0,@sorted_hspBNs);
    my %geneinfo=%{$reginfors{$sName}};
    my $chrom=$geneinfo{'chrom'};
    my @checkedexons=&exoncheck($genestrand,@refExon);
    my $coverlength=shift(@checkedexons);
    if(scalar(@checkedexons)<1){
	next;
    }
    my $coverpercent=substr((($coverlength/$resBN->query_length())*100),0,4);
    my ($chrS,$chrE,$exonnum)=&chromregion(@checkedexons);
    $genemapping{$qName}->{'chrom'}=$chrom;
    $genemapping{$qName}->{'strand'}=$genestrand;
    $genemapping{$qName}->{'geneid'}=$sName;
    $genemapping{$qName}->{'exons'}=\@refExon;
    $genemapping{$qName}->{'chrS'}=$chrS;
    $genemapping{$qName}->{'chrE'}=$chrE;
    $genemapping{$qName}->{'exonnum'}=$exonnum;
    $genemapping{$qName}->{'length'}=$resBN->query_length();
    $genemapping{$qName}->{'coverage'}=$coverpercent;
}

while(my $resTBX=$searchTBX->next_result()){
    $qName=$resTBX->query_name();
    if($resTBX->hits==0){
	next;
    }
    my $blasthits=$resTBX->hits;
#    print "tblastx qName:$qName\n";

    my $hitTBX=$resTBX->next_hit();
    my $sName=$hitTBX->name();
    #k33:1021589
#    if($qName eq "k33:1020637u"){
#	sleep 1;
#    }
    my @exons;
    my @hspTBXs;
    my $hspcount=0;
    my $genestrand;
    while(my $hspTBX=$hitTBX->next_hsp()){
	my $strand;
	my $qstrand=$hspTBX->strand('query');
	my $tstrand=$hspTBX->strand('hit');
	$strand=$qstrand*$tstrand;
	if($genestrand eq ""){
	    $genestrand=$strand;
	}
	if(($hspTBX->length('total')*3)<$parMinLength || $hspTBX->frac_identical('total')<$parMinIdentity || $genestrand ne $strand){
	    next;
	}
	if($hspTBX->{QUERY_FRAME}<0){
	    my $tempS=$hspTBX->{QUERY_START};
	    $hspTBX->{QUERY_START}=$hspTBX->{QUERY_END};
	    $hspTBX->{QUERY_END}=$tempS;
	}
	if($hspTBX->{HIT_FRAME}<0){
	    my $tempS=$hspTBX->{HIT_START};
	    $hspTBX->{HIT_START}=$hspTBX->{HIT_END};
	    $hspTBX->{HIT_END}=$tempS;
	}
	$hspTBXs[$hspcount]=$hspTBX;
	$hspcount++;
    }

    #my $genestrand;
    if($genestrand eq "1"){$genestrand="+";}
    else{$genestrand="-";}
    #$genestrand==1 ? $genestrand='+' : $genestrand='-';

    if(scalar(@hspTBXs)==0){
	next;
    }

    my @sorted_hspTBXs;
    if(scalar(@hspTBXs)<2){
	@sorted_hspTBXs=@hspTBXs;
    }
    else{
	@sorted_hspTBXs=hspsort(@hspTBXs);
    }

    my @refExon=overlapcheck($sName,$genestrand,1,@sorted_hspTBXs);
    my %geneinfo=%{$reginfors{$sName}};
    my $chrom=$geneinfo{'chrom'};
    my @checkedexons=&exoncheck($genestrand,@refExon);
    my $coverlength=shift(@checkedexons);
    if(scalar(@checkedexons)<1){
	next;
    }
    my $coverpercent=substr((($coverlength/$resTBX->query_length())*100),0,4);
    my ($chrS,$chrE,$exonnum)=&chromregion(@checkedexons);
    if(!(exists($genemapping{$qName}))){
	$genemapping{$qName}->{'chrom'}=$chrom;
	$genemapping{$qName}->{'exons'}=\@checkedexons;
	$genemapping{$qName}->{'strand'}=$genestrand;
	$genemapping{$qName}->{'geneid'}=$sName;
	$genemapping{$qName}->{'chrS'}=$chrS;
	$genemapping{$qName}->{'chrE'}=$chrE;
	$genemapping{$qName}->{'exonnum'}=$exonnum;
	$genemapping{$qName}->{'length'}=$resTBX->query_length();
	$genemapping{$qName}->{'coverage'}=$coverpercent;
    }
    else{
	my %tbxinfo;
	$tbxinfo{'chrom'}=$chrom;
	$tbxinfo{'exons'}=\@refExon;
	$tbxinfo{'strand'}=$genestrand;
	$tbxinfo{'geneid'}=$sName;
	$tbxinfo{'chrS'}=$chrS;
	$tbxinfo{'chrE'}=$chrE;
	$tbxinfo{'exonnum'}=$exonnum;
	&blastcompare($qName,\%tbxinfo);
    }
}

open (my $fileout,">$FOut");
##output result;
foreach my $key(keys %genemapping){
    my %geneinfo=%{$genemapping{$key}};
    my @exons=@{$geneinfo{'exons'}};
    my $length=$geneinfo{'length'};
    my $chrom=$geneinfo{'chrom'};
    my $genestrand=$geneinfo{'strand'};
    my $geneid=$geneinfo{'geneid'};
    my $chrS=$geneinfo{'chrS'};
    my $chrE=$geneinfo{'chrE'};
    my $exonnum=$geneinfo{'exonnum'};
    my $coverage=$geneinfo{'coverage'};
    my $querys;
    my $targets;
    for(my $j=0;$j<@exons;$j++){
	my %exoninfo=%{$exons[$j]};
	my $qs=$exoninfo{'qs'};
	my $qe=$exoninfo{'qe'};
	my $ts=$exoninfo{'ts'};
	my $te=$exoninfo{'te'};
	$querys=$querys.$qs.":".$qe.",";
	$targets=$targets.$ts.":".$te.",";
    }
    print $fileout "$key\t$length\t$coverage\t$chrom\t$genestrand\t$chrS\t$chrE\t$exonnum\t$querys\t$targets\n";
}
close $fileout;
exit(0);



##sub functions
sub blastcompare{
    my ($qName,$reftbxinfo)=@_;
    my %blninfo=%{$genemapping{$qName}};
    my %tbxinfo=%{$reftbxinfo};
    if($blninfo{'strand'} ne $tbxinfo{'strand'}){
	#print "Error:DDDD\n";
	return 0;
    }
    my $strand=$blninfo{'strand'};
    if(($blninfo{'chrS'}>=$tbxinfo{'chrS'} && $blninfo{'chrS'}<=$tbxinfo{'chrE'}) || ($tbxinfo{'chrS'}>=$blninfo{'chrS'} && $tbxinfo{'chrS'}<=$blninfo{'chrE'})){
	my @blninfos=@{$blninfo{'exons'}};
	my @tbxinfos=@{$tbxinfo{'exons'}};
	my $lengthbln=scalar(@blninfos);
	my $lengthtbx=scalar(@tbxinfos);
	my $blncount=0;
	my $tbxcount=0;
	my $exoncount=0;
	my @exons;
	while($blncount<$lengthbln && $tbxcount<$lengthtbx){
	    my $blnstart=${$blninfos[$blncount]}{'qs'};
	    my $blnend=${$blninfos[$blncount]}{'qe'};
	    my $blnhitstart=${$blninfos[$blncount]}{'ts'};
	    my $blnhitend=${$blninfos[$blncount]}{'te'};
	    my $tbxstart=${$tbxinfos[$tbxcount]}{'qs'};
	    my $tbxend=${$tbxinfos[$tbxcount]}{'qe'};
	    my $tbxhitstart=${$tbxinfos[$tbxcount]}{'ts'};
	    my $tbxhitend=${$tbxinfos[$tbxcount]}{'te'};
	    if(scalar(@exons)>0){
		my %lastexon=%{$exons[-1]};
		my $lastqs=$lastexon{'qs'};
		my $lastqe=$lastexon{'qe'};
		if($lastqe>=$tbxend){
		    $tbxcount++;
		    next;
		}
		else{
		    if($tbxstart<=$lastqe){
			my $overlength=$lastqe-$tbxstart+1;
			if($strand eq "+"){
			    $tbxstart=$lastqe+1;
			    $tbxhitstart=$tbxhitstart+$overlength;
			}
			else{
			    $tbxstart=$lastqe+1;
			    $tbxhitend=$tbxhitend-$overlength;
			}
		    }
		}
	    }

	    if(($blnstart<=$tbxstart && $blnend>=$tbxstart) || ($tbxstart<=$blnstart && $tbxend>=$blnstart)){ ##overlap 1
		if(($blnhitstart<=$tbxhitstart && $tbxhitend>=$tbxhitstart) || ($tbxhitstart<=$blnhitstart && $tbxhitend>=$blnhitstart)){ ##determing whether hit overlap
		    if($blnstart<=$tbxstart && $tbxend>=$tbxstart){ ##blastn less blastx
			if($blnend>=$tbxend){ ##blastx included in blastn
			    my %exoninfo;
			    $exoninfo{'qs'}=$blnstart;
			    $exoninfo{'qe'}=$blnend;
			    $exoninfo{'ts'}=$blnhitstart;
			    $exoninfo{'te'}=$blnhitend;
			    push(@exons,\%exoninfo);
			}
			else{ ##tbxend > blnend , overlap check next blninfo;
			    if(exists($blninfos[$blncount+1])){
				my $blnstart2=${$blninfos[$blncount+1]}{'qs'};
				if($blnstart2<$tbxend){
				    my $overlength=$tbxend-$blnstart2+1;
				    $tbxend=$blnstart2-1;
				    my %exoninfo;
				    if($strand eq "+"){
					$tbxhitend=$tbxhitend-$overlength;
					$exoninfo{'qs'}=$blnstart;
					$exoninfo{'qe'}=$tbxend;
					$exoninfo{'ts'}=$blnhitstart;
					$exoninfo{'te'}=$tbxhitend;
				    }
				    else{
					$tbxhitstart=$tbxhitstart+$overlength;
					$exoninfo{'qs'}=$blnstart;
					$exoninfo{'qe'}=$tbxend;
					$exoninfo{'ts'}=$tbxhitstart;
					$exoninfo{'te'}=$blnhitend;
				    }
				    push(@exons,\%exoninfo);

				}
				else{
				    my %exoninfo;
				    if($strand eq "+"){##strand
					$exoninfo{'qs'}=$blnstart;
					$exoninfo{'qe'}=$tbxend;
					$exoninfo{'ts'}=$blnhitstart;
					$exoninfo{'te'}=$tbxhitend;
					push(@exons,\%exoninfo);
				    }
				    else{
					$exoninfo{'qs'}=$blnstart;
					$exoninfo{'qe'}=$tbxend;
					$exoninfo{'ts'}=$tbxhitstart;
					$exoninfo{'te'}=$blnhitend;
					push(@exons,\%exoninfo);
				    }
				}
			    }
			    else{
				my %exoninfo;
				if($strand eq "+"){##strand
				    $exoninfo{'qs'}=$blnstart;
				    $exoninfo{'qe'}=$tbxend;
				    $exoninfo{'ts'}=$blnhitstart;
				    $exoninfo{'te'}=$tbxhitend;
				    push(@exons,\%exoninfo);
				}
				else{
				    $exoninfo{'qs'}=$blnstart;
				    $exoninfo{'qe'}=$tbxend;
				    $exoninfo{'ts'}=$tbxhitstart;
				    $exoninfo{'te'}=$blnhitend;
				    push(@exons,\%exoninfo);
				}
			    }
			}
		    }
		    elsif($tbxstart<=$blnstart && $tbxend>=$blnstart){##tblastx less blastn
			my %exoninfo;
			if($blnend>=$tbxend){ ##blastn longer than tblastx
			    if(($blncount-1)>=0){
				my $blnend2=${$blninfos[$blncount-1]}{'qe'};
				if($blnend2>$tbxstart){ ##overlap with previous blastn
				    my $overlength=$tbxstart-$blnend2+1;
				    $tbxstart=$blnend2+1;
				    if($strand eq "+"){
					$tbxhitstart=$tbxhitstart+$overlength;
					$exoninfo{'qs'}=$tbxstart;
					$exoninfo{'qe'}=$blnend;
					$exoninfo{'ts'}=$tbxhitstart;
					$exoninfo{'te'}=$blnhitend;
				    }
				    else{
					$tbxhitend=$tbxhitend-$overlength;
					$exoninfo{'qs'}=$tbxstart;
					$exoninfo{'qe'}=$blnend;
					$exoninfo{'ts'}=$blnhitstart;
					$exoninfo{'te'}=$tbxhitend;
				    }
				}
				else{
				    if($strand eq "+"){
					$exoninfo{'qs'}=$tbxstart;
					$exoninfo{'qe'}=$blnend;
					$exoninfo{'ts'}=$tbxhitstart;
					$exoninfo{'te'}=$blnhitend;
				    }
				    else{
					$exoninfo{'qs'}=$tbxstart;
					$exoninfo{'qe'}=$blnend;
					$exoninfo{'ts'}=$blnhitstart;
					$exoninfo{'te'}=$tbxhitend;
				    }
				}
			    }
			    else{
				if($strand eq "+"){
				    $exoninfo{'qs'}=$tbxstart;
				    $exoninfo{'qe'}=$blnend;
				    $exoninfo{'ts'}=$tbxhitstart;
				    $exoninfo{'te'}=$blnhitend;
				}
				else{
				    $exoninfo{'qs'}=$tbxstart;
				    $exoninfo{'qe'}=$blnend;
				    $exoninfo{'ts'}=$blnhitstart;
				    $exoninfo{'te'}=$tbxhitend;
				}
			    }
			}
			else{##blastn included in tblastx
			    if(($blncount-1)>=0){
				my ($blnstart2,$blnend2);
				$blnend2=${$blninfos[$blncount-1]}{'qe'};
				if(exists($blninfos[$blncount+1])){
				    $blnstart2=${$blninfos[$blncount-1]}{'qs'};
				    if($blnend2>$tbxstart){
					my $overlength1=$blnend2-$tbxstart+1;
					if($blnstart2<$tbxend){
					    my $overlength2=$tbxend-$blnstart2+1;
					    $tbxstart=$blnend2+1;
					    $tbxend=$blnstart2-1;
					    if($strand eq "+"){
						$tbxhitstart=$tbxhitstart+$overlength1;
						$tbxhitend=$tbxhitend-$overlength2;
					    }
					    else{
						$tbxhitstart=$tbxhitstart+$overlength2;
						$tbxhitend=$tbxhitend-$overlength1;

					    }
					    $exoninfo{'qs'}=$tbxstart;
					    $exoninfo{'qe'}=$tbxend;
					    $exoninfo{'ts'}=$tbxhitstart;
					    $exoninfo{'te'}=$tbxhitend;
					}
					else{ ##blnstart2 > tbxend
					    $tbxstart=$blnend2+1;
					    if($strand eq "+"){
						$tbxhitstart=$tbxhitstart+$overlength1;
					    }
					    else{
						$tbxhitend=$tbxhitend-$overlength1;
					    }
					    $exoninfo{'qs'}=$tbxstart;
					    $exoninfo{'qe'}=$tbxend;
					    $exoninfo{'ts'}=$tbxhitstart;
					    $exoninfo{'te'}=$tbxhitend;
					}
				    }
				    else{ ## blnend2< tbxend -> no overlap
					if($blnstart2<$tbxend){
					    my $overlength2=$tbxend-$blnstart2+1;
					    $tbxend=$blnstart2;
					    if($strand eq "+"){
						$tbxhitend=$tbxhitend-$overlength2;
					    }
					    else{
						$tbxhitstart=$tbxhitstart+$overlength2;
					    }
					    $exoninfo{'qs'}=$tbxstart;
					    $exoninfo{'qe'}=$tbxend;
					    $exoninfo{'ts'}=$tbxhitstart;
					    $exoninfo{'te'}=$tbxhitend;
					}
					else{
					    $exoninfo{'qs'}=$tbxstart;
					    $exoninfo{'qe'}=$tbxend;
					    $exoninfo{'ts'}=$tbxhitstart;
					    $exoninfo{'te'}=$tbxhitend;
					}
				    }
				}
				else{
				    if($blnend2>$tbxstart){
					my $overlength=$tbxstart-$blnend2+1;
					$tbxstart=$blnend2+1;
					if($strand eq "+"){
					    $tbxhitstart=$tbxhitstart+$overlength;
					}
					else{
					    $tbxhitend=$tbxhitend-$overlength;
					}
					$exoninfo{'qs'}=$tbxstart;
					$exoninfo{'qe'}=$tbxend;
					$exoninfo{'ts'}=$tbxhitstart;
					$exoninfo{'te'}=$tbxhitend;
				    }
				    else{
					$exoninfo{'qs'}=$tbxstart;
					$exoninfo{'qe'}=$tbxend;
					$exoninfo{'ts'}=$tbxhitstart;
					$exoninfo{'te'}=$tbxhitend;
				    }
				}
			    }
			    else{
				if(exists($blninfos[$blncount+1])){
				    my $blnstart2=${$blninfos[$blncount-1]}{'qs'};
				    if($blnstart2<$tbxend){
					my $overlength=$tbxend-$blnstart2+1;
					$tbxend=$blnstart2-1;
					if($strand eq "+"){
					    $tbxhitend=$tbxhitend-$overlength;
					}
					else{
					    $tbxhitstart=$tbxhitstart+$overlength;
					}
					$exoninfo{'qs'}=$tbxstart;
					$exoninfo{'qe'}=$tbxend;
					$exoninfo{'ts'}=$tbxhitstart;
					$exoninfo{'te'}=$tbxhitend;
				    }
				    else{
					$exoninfo{'qs'}=$tbxstart;
					$exoninfo{'qe'}=$tbxend;
					$exoninfo{'ts'}=$tbxhitstart;
					$exoninfo{'te'}=$tbxhitend;
				    }
				}
				else{
				    $exoninfo{'qs'}=$tbxstart;
				    $exoninfo{'qe'}=$tbxend;
				    $exoninfo{'ts'}=$tbxhitstart;
				    $exoninfo{'te'}=$tbxhitend;
				}
			    }

			}
			push(@exons,\%exoninfo);
		    }
		}
		else{ ##blast hit not overlap return blastn result;
		    my %exoninfo;
		    $exoninfo{'qs'}=$blnstart;
		    $exoninfo{'qe'}=$blnend;
		    $exoninfo{'ts'}=$blnhitstart;
		    $exoninfo{'te'}=$blnhitend;
		    push(@exons,\%exoninfo);
		}
		$blncount++;
		$tbxcount++;
	    }
	    else{
		if($tbxend<=$blnstart){
		    if(($blncount-1)>=0){
			my $blnend2=${$blninfos[$blncount-1]}{'qe'};
			if($blnend2>$tbxstart){ ##overlap with previous one
			    my $overlength=$blnend2-$tbxstart+1;
			    $tbxhitend=$tbxhitend-$overlength;
			    my %exoninfo;
			    $exoninfo{'qs'}=$blnend2+1;
			    $exoninfo{'qe'}=$tbxend;
			    $exoninfo{'ts'}=$tbxhitstart;
			    $exoninfo{'te'}=$tbxhitend;
			    push(@exons,\%exoninfo);
			}
			else{
			    my %exoninfo;
			    $exoninfo{'qs'}=$tbxstart;
			    $exoninfo{'qe'}=$tbxend;
			    $exoninfo{'ts'}=$tbxhitstart;
			    $exoninfo{'te'}=$tbxhitend;
			    push(@exons,\%exoninfo);
			}
		    }
		    else{
			my %exoninfo;
			$exoninfo{'qs'}=$tbxstart;
			$exoninfo{'qe'}=$tbxend;
			$exoninfo{'ts'}=$tbxhitstart;
			$exoninfo{'te'}=$tbxhitend;
			push(@exons,\%exoninfo);
		    }
		    $tbxcount++;
		}
		elsif($blnstart<=$tbxend){
		    my %exoninfo;
		    $exoninfo{'qs'}=$blnstart;
		    $exoninfo{'qe'}=$blnend;
		    $exoninfo{'ts'}=$blnhitstart;
		    $exoninfo{'te'}=$blnhitend;
		    push(@exons,\%exoninfo);
		    $blncount++;
		}
		else{
		    print "Error:FFFF\n";
		}
	    }
	}
	if($blncount<$lengthbln){
	    while($blncount<$lengthbln){
		my %exoninfo;
		#my %blninfo=%{$blninfos[$blncount]};
		#$exoninfo{'start'}=${$blninfos[$blncount]}{'qs'};
		$exoninfo{'qs'}=${$blninfos[$blncount]}{'qs'};
		$exoninfo{'qe'}=${$blninfos[$blncount]}{'qe'};
		$exoninfo{'ts'}=${$blninfos[$blncount]}{'ts'};
		$exoninfo{'te'}=${$blninfos[$blncount]}{'te'};
		push(@exons,\%exoninfo);
		$blncount++;
	    }
	}

	if($tbxcount<$lengthtbx){
	    my %lastexon=%{$exons[-1]};
	    my $lastexonE=$lastexon{'qe'};
	    while($tbxcount<$lengthtbx){
		my %exoninfo;
		my $tbxstart=${$tbxinfos[$tbxcount]}{'qs'};
		my $tbxend=${$tbxinfos[$tbxcount]}{'qe'};
		my $tbxhitstart=${$tbxinfos[$tbxcount]}{'ts'};
		my $tbxhitend=${$tbxinfos[$tbxcount]}{'te'};
		if($lastexonE>=$tbxend){
		    $tbxcount++;
		    next;
		}
		if($lastexonE>=$tbxstart){
		    my $overlength=$lastexonE-$tbxstart+1;
		    if($strand eq "+"){
			$tbxhitend=$tbxhitend-$overlength;
		    }
		    else{
			$tbxhitstart=$tbxhitstart+$overlength;
		    }

		    $exoninfo{'qs'}=$lastexonE+1;
		    $exoninfo{'qe'}=$tbxend;
		    $exoninfo{'ts'}=$tbxhitstart;
		    $exoninfo{'te'}=$tbxhitend;
		}
		else{
		    $exoninfo{'qs'}=$tbxstart;
		    $exoninfo{'qe'}=$tbxend;
		    $exoninfo{'ts'}=$tbxhitstart;
		    $exoninfo{'te'}=$tbxhitend;
		}
		push(@exons,\%exoninfo);
		$tbxcount++;
	    }
	}
	my @checkedexons=&exoncheck($strand,@exons);
	my $coverlength=shift(@checkedexons);
	if(scalar(@checkedexons)<1){
	    return 0;
	}
	my $coverpercent=substr((($coverlength/$genemapping{$qName}->{'length'})*100),0,4);
	my ($chrS,$chrE,$exonnum)=&chromregion(@checkedexons);
	$genemapping{$qName}->{'exons'}=\@checkedexons;
	$genemapping{$qName}->{'chrS'}=$chrS;
	$genemapping{$qName}->{'chrE'}=$chrE;
	$genemapping{$qName}->{'exonnum'}=$exonnum;
	$genemapping{$qName}->{'coverage'}=$coverpercent;
    }
    else{
	#print "Error:EEEE\n";
	return 0;
    }
}

sub chromregion{
    my @exons=@_;
    my %exon1=%{$exons[0]};
    my $chrS=$exon1{'ts'};
    my $arraylength=scalar(@exons);
    my %exonE=%{$exons[$arraylength-1]};
    my $chrE=$exonE{'te'};
    for(my $i=0;$i<@exons;$i++){
	my $start=${$exons[$i]}{'ts'};
	my $end=${$exons[$i]}{'te'};
	if($start<$chrS){
	    $chrS=$start;
	}
	if($end>$chrE){
	    $chrE=$end;
	}
    }

    return ($chrS,$chrE,$arraylength);
}

sub overlapcheck{
    my @hsps=@_;
    my $geneid=shift(@hsps);
    my $strand=shift(@hsps);
    my $mode=shift(@hsps);
    my $exoncount=0;
    #my %exoninfo;

    ##remove all redundant exons;
    for(my $i=0;$i<@hsps-1;$i++){
	#Reloop:
	#print "$i\n";
	my $reg1s=$hsps[$i]->query->start;
	my $reg1e=$hsps[$i]->query->end;
	my $reg2s=$hsps[$i+1]->query->start;
	my $reg2e=$hsps[$i+1]->query->end;
	my $overregion;
	my $reg1length=$reg1e-$reg1s;
	my $reg2length=$reg2e-$reg2s;
	my $overS;
	my $overE;
	my $overmode=0;
	my $larger;
	if($reg1s<=$reg2s && $reg1e>=$reg2s){ ##overlap 1
	    if($reg1e<=$reg2e){
		$overregion=$reg1e-$reg2s;
		$overmode=1;
	    }
	    else{
		$overmode=2;
		$larger=1;
	    }
	}
	elsif($reg2s<=$reg1s && $reg2e>=$reg1s){ ##overlap 2
	    if($reg2e<=$reg1e){
		$overregion=$reg2e-$reg1s;
		$overmode=1;
	    }
	    else{
		$overmode=2;
		$larger=2;
	    }
	}
	my $overper1=$overregion/$reg1length;
	my $overper2=$overregion/$reg2length;
	if($overmode==2){
	    if($larger==1){
		splice(@hsps,$i+1,1);
		$i=-1; ##remove second record, and move i to 0;
		#goto Reloop;
	    }
	    elsif($larger==2){
		splice(@hsps,$i,1);
		$i=-1;
		#goto Reloop;
	    }
	}
	elsif($overmode==1){
	    if($overper1 >= $parRdnOverpercent || $overper2 >= $parRdnOverpercent){
		my %geneinfo=%{$reginfors{$geneid}};
		my $center=int(($geneinfo{'end'}-$geneinfo{'start'}+2*$parFlankingLength)/2);
		my $hitstart1=$hsps[$i]->hit->start;
		my $hitend1=$hsps[$i]->hit->end;
		my $hitstart2=$hsps[$i+1]->hit->start;
		my $hitend2=$hsps[$i+1]->hit->end;
		my $hitcenter1=int(($hitstart1+$hitend1)/2);
		my $hitcenter2=int(($hitstart2+$hitend2)/2);

		if(($hitstart1<=$hitstart2 && $hitend1>=$hitstart2) || ($hitstart2<=$hitstart1 && $hitend2>=$hitstart1)){
		    if($hsps[$i]->percent_identity >= $hsps[$i+1]->percent_identity){
			splice(@hsps,$i+1,1);
			$i=-1;
			#goto Reloop;
		    }
		    else{
			splice(@hsps,$i,1);
			$i=-1;
			#goto Reloop;
		    }
		}
		else{
		    if (abs($hitcenter1-$center)<=abs($hitcenter2-$center)){
		        splice(@hsps,$i+1,1);
		        $i=-1;
		        #goto Reloop;
		    }
		    else{
		        splice(@hsps,$i,1);
		        $i=-1;
		        #goto Reloop;
		    }
		}
	    }
	}
    }

    ##output exons and deal with overlap
    my %exonregion;
    my @exons;
    my %geneinfo=%{$reginfors{$geneid}};
    my $chromS=$geneinfo{'start'}-$parFlankingLength;
    if($chromS <0 ){$chromS=1;}

    for(my $i=0;$i<@hsps;$i++){
	my ($reg1s,$reg1e,$hit1s,$hit1e);
	if(exists $exonregion{$i}){
	    $reg1s=$exonregion{$i}->{'start'};
	    $reg1e=$exonregion{$i}->{'end'};
	    $hit1s=$exonregion{$i}->{'hitstart'};
	    $hit1e=$exonregion{$i}->{'hitend'};
	}
	else{
	    $reg1s=$hsps[$i]->query->start;
	    $reg1e=$hsps[$i]->query->end;
	    $hit1s=$hsps[$i]->hit->start;
	    $hit1e=$hsps[$i]->hit->end;
	}

	if($i == (scalar(@hsps)-1)){
	    my %exoninfo;
	    $exoninfo{'qs'}=$reg1s;
	    $exoninfo{'qe'}=$reg1e;
	    $exoninfo{'ts'}=$hit1s+$chromS-1;
	    $exoninfo{'te'}=$hit1e+$chromS-1;
	    push(@exons,\%exoninfo);
	    last;
	}

	my $reg2s=$hsps[$i+1]->query->start;
	my $reg2e=$hsps[$i+1]->query->end;
	my $overregion;
	my $reg1length=$reg1e-$reg1s;
	my $reg2length=$reg2e-$reg2s;
	my $overS;
	my $overE;

	if($reg1s<=$reg2s && $reg1e>=$reg2s){ ##overlap
	    if($reg1e<=$reg2e){
		$overS=$reg2s; #128
		$overE=$reg1e; #144
		my $alignhsp1=$hsps[$i]->get_aln();
		my $alignhsp2=$hsps[$i+1]->get_aln();
		my $matchline1=$alignhsp1->match_line();
		my $matchline2=$alignhsp2->match_line();
		#my $alignStart1=$overS-$reg1s;
		#my $alignStart2=$overE-$reg2s;
		my $overlength=$overE-$overS+1;
		my ($matchsub1,$matchsub2);
		if($mode==1){
		    my $overlengtht1=int($overlength/3);
		    my $overlengtht2=0-$overlengtht1;
		    if($hsps[$i]->{QUERY_FRAME}<0){
			$matchline1=reverse($matchline1);
		    }
		    if($hsps[$i+1]->{QUERY_FRAME}<0){
			$matchline2=reverse($matchline2);
		    }
		    $matchsub1=substr($matchline1,$overlengtht2);
		    $matchsub2=substr($matchline2,0,$overlengtht1);
		}
		else{
		    my $overlengtht=0-$overlength;
		    $matchsub1=substr($matchline1,$overlengtht);
		    $matchsub2=substr($matchline2,0,$overlength);
		}
		my $greaterscore=&scorepare($matchsub1,$matchsub2);
		if($greaterscore==1){
		    my $j=$i+1;
		    $reg2s=$reg2s+$overlength;
		    if($strand eq "+"){
			my $hitStart2=$hsps[$i+1]->hit->start+$overlength;
			$exonregion{$j}->{'start'}=$reg1e+1;
			$exonregion{$j}->{'end'}=$reg2e;
			$exonregion{$j}->{'hitstart'}=$hitStart2;
			$exonregion{$j}->{'hitend'}=$hsps[$i+1]->hit->end;
		    }
		    else{
			my $hitEnd2=$hsps[$i+1]->hit->end-$overlength;
			$exonregion{$j}->{'start'}=$reg1e+1;
			$exonregion{$j}->{'end'}=$reg2e;
			$exonregion{$j}->{'hitstart'}=$hsps[$i+1]->hit->start;
			$exonregion{$j}->{'hitend'}=$hitEnd2;
		    }
		    my %exoninfo;
		    $exoninfo{'qs'}=$reg1s;
		    $exoninfo{'qe'}=$reg1e;
		    $exoninfo{'ts'}=$chromS+$hit1s-1;
		    $exoninfo{'te'}=$chromS+$hit1e-1;
		    push(@exons,\%exoninfo);
		}
		else{
		    $reg1e=$reg1e-$overlength;
		    if($strand eq "+"){
			$hit1e=$hit1e-$overlength;
		    }
		    else{
			$hit1s=$hit1s+$overlength;
		    }
		    my %exoninfo;
		    $exoninfo{'qs'}=$reg1s;
		    $exoninfo{'qe'}=$reg1e;
		    $exoninfo{'ts'}=$chromS+$hit1s-1;
		    $exoninfo{'te'}=$chromS+$hit1e-1;
		    push(@exons,\%exoninfo);
		}
	    }
	    else{
		print "Error!AAAA\n";
	    }
	}
	else{
	    if($reg2s>=$reg1e){
		my %exoninfo;
		$exoninfo{'qs'}=$reg1s;
		$exoninfo{'qe'}=$reg1e;
		$exoninfo{'ts'}=$chromS+$hit1s-1;
		$exoninfo{'te'}=$chromS+$hit1e-1;
		push(@exons,\%exoninfo);
	    }
	    else{
		my %exoninfo;
		if($reg1e>=$reg2e){
		    $exoninfo{'qs'}=$reg1s;
		    $exoninfo{'qe'}=$reg1e;
		    $exoninfo{'ts'}=$chromS+$hit1s-1;
		    $exoninfo{'te'}=$chromS+$hit1e-1;
		    push(@exons,\%exoninfo);
		    $i++;
		}
		else{
		    $exoninfo{'qs'}=$reg1s;
		    $exoninfo{'qe'}=$reg1e;
		    $exoninfo{'ts'}=$chromS+$hit1s-1;
		    $exoninfo{'te'}=$chromS+$hit1e-1;
		    push(@exons,\%exoninfo);
		    my $overlength=$reg2e-$reg1e+1;
		    my $hitEnd2=$hsps[$i+1]->hit->end-$overlength;
		    my $j=$i+1;
		    $exonregion{$j}->{'start'}=$reg1e+1;
		    $exonregion{$j}->{'end'}=$reg2e;
		    $exonregion{$j}->{'hitstart'}=$hsps[$i+1]->hit->start+$overlength;
		    $exonregion{$j}->{'hitend'}=$hitEnd2;
		}
	    }
	}
    }

    return @exons;

}


sub scorepare{
    my ($line1,$line2)=@_;
    my @match1=split("",$line1);
    my @match2=split("",$line2);
    my ($score1,$score2);
    foreach my $chr(@match1){
	if($chr eq "*"){$score1=$score1+2;}
	elsif($chr eq ":"){$score1=$score1+1;}
	elsif($chr eq "."){$score1=$score1+0.5;}
    }
    foreach my $chr(@match2){
	if($chr eq "*"){$score2=$score2+2;}
	elsif($chr eq ":"){$score2=$score2+1;}
	elsif($chr eq "."){$score2=$score2+0.5;}
    }
    if($score1>=$score2){
	return 1;
    }
    else{
	return 2;
    }
}

sub hspsort{
    my @hsps=@_;
    my @sorted;
    while(scalar(@hsps)>0){
	my $minindex;
	my $minstart;
	$minstart=$hsps[0]->query->start;
	for(my $i=1;$i<@hsps;$i++){
	    my $qstart=$hsps[$i]->query->start;
	    if($qstart<$minstart){
		$minstart=$qstart;
		$minindex=$i;
	    }
	}
	push(@sorted,$hsps[$minindex]);
	splice(@hsps,$minindex,1);
    }
    return @sorted;
}

sub exoncheck{
    my @exons=@_;
    my $strand=shift(@exons);

    my $coverlength=0;
    for(my $i=0;$i<@exons;$i++){
	my %exoninfo=%{$exons[$i]};
	my $qs=$exoninfo{'qs'};
	my $qe=$exoninfo{'qe'};
	my $ts=$exoninfo{'ts'};
	my $te=$exoninfo{'te'};
	if(($qe-$qs+1)<$parMinLength){
	    splice(@exons,$i,1);
	    $i=$i-1;
	}
	elsif($te<$ts){
	    splice(@exons,$i,1);
	    $i=$i-1;
	}
	else{
	    $coverlength=$coverlength+($qe-$qs+1);
	}
    }

    return ($coverlength,@exons);
}

sub checkfile{
    my $filename=shift;
    if(!(-e $filename)){
	print stderr "Can't find file: $filename\n";
	$iserror=1;
    }
}
