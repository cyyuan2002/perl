#!/home/SCE/chenyuan/bin/Perl-5.12/bin/perl

use strict;
use Bio::Tools::Run::RemoteBlast;
use Bio::SearchIO;
use Bio::SeqIO;
use Getopt::Long;
my %opts;
my $version="1.0 Alvin Chen 2011-07-24";

GetOptions(\%opts,"i=s","d:s","e:s","p=s","r:n","o=s","help");
if((!defined $opts{i}) || (!defined $opts{p}) ||(!defined $opts{o})){
    &Usage();
}

my $parInputFile=$opts{i};
my $parDatabase=(defined $opts{d}) ? $opts{d} : "nr";
my $parProgram=$opts{p};
my $parOutputFile=$opts{o};
my $parEvalue=(defined $opts{e}) ? $opts{e} : "1e-6";
my $parResultNum=(defined $opts{r}) ? $opts{r}:5;
my $iserror=0;

&Filecheck($parInputFile,0);
if($iserror==1){
	exit(1);
}

my @params=('-prog' => $parProgram, '-data'=>$parDatabase,'-expect'=>$parEvalue,'-readmethod'=>'SearchIO');
my $factory=Bio::Tools::Run::RemoteBlast->new(@params);
$Bio::Tools::Run::RemoteBlast::retrieve_parameter{'DESCRIPTIONS'}=$parResultNum;

my $v=1;
my $sequences=Bio::SeqIO->new(-file=>$parInputFile,-format=>'fasta');
open(my $fh_fileout,">$parOutputFile");
my $blastcount=1;

while(my $seq=$sequences->next_seq()){
	my $result=$factory->submit_blast($seq);
	print STDERR "Num:$blastcount  Query: ",$seq->id," wait...";
	$blastcount++;
	while(my @rids=$factory->each_rid){
		foreach my $rid(@rids){
			my $rc=$factory->retrieve_blast($rid);
			if(!ref($rc)){
				if($rc<0){
					$factory->remove_rid($rc);
				}
				print STDERR "." if ( $v > 0 );
				sleep 3;
			}
			else{
				print STDERR "\n";
				my $res=$rc->next_result();				
				my @blastinfo=&FormatResult($res);
				if(@blastinfo){
					for(my $i=0;$i<@blastinfo;$i++){
						my %blastRes=%{$blastinfo[$i]};
						my $qName=$blastRes{'query'};
						my $sName=$blastRes{'subject'};
						my $qLength=$blastRes{'qlen'};
						my $sLength=$blastRes{'slen'};
						my $qStart=$blastRes{'queryS'};
						my $qEnd=$blastRes{'queryE'};
						my $sStart=$blastRes{'chrS'};
						my $sEnd=$blastRes{'chrE'};
						my $evalue=$blastRes{'evalue'};
						my $des=$blastRes{'description'};
						my $strand=$blastRes{'strand'};
						my $identity=$blastRes{'identity'};
						my $score=$blastRes{'score'};
						print $fh_fileout "$qName\t$qLength\t$qStart\t$qEnd\t$sStart\t$sEnd\t$sLength\t$score\t$evalue\t$strand\t$identity\t$sName\t$des\n";	
					
					}
				}
				$factory->remove_rid($rid);
			}
		}
	}
}
exit(0);

sub FormatResult{
	my $resBN=shift;
	my $qName=$resBN->query_name();
	if(scalar($resBN)->hits()==0){
		return;
	}
	my $hitscount=0;
	my @blastRes;
	while (my $hitBN=$resBN->next_hit()){
		last if($hitscount==$parResultNum);
	    my $sName=$hitBN->name();
	    my $hspcount=0;
	    my $strand;
	    my $evalue;
	    my $identity;
	    my $score;
	    
	    my $hspBN=$hitBN->next_hsp();
	    if($hspBN->strand('hit')==0){
	    	$strand=$hspBN->strand('query');
	    }
		elsif($hspBN->strand('query')==0){
			$strand=$hspBN->strand('hit');
		}
		else{
			$strand=$hspBN->strand('query')*$hspBN->strand('hit');
		}
		my $genestrand;
		
		if($strand > 0){$genestrand="+";}
    	else{$genestrand="-";}
		
		my $queryS=$hspBN->query->start;
		my $queryE=$hspBN->query->end;
		my $chrS=$hspBN->hit->start;
		my $chrE=$hspBN->hit->end;
		
		$evalue=$hspBN->evalue();
		$identity=$hspBN->frac_identical('total');
		$identity=sprintf("%.2f",$identity);
		$score=$hspBN->score();
  
	    $blastRes[$hitscount]->{'query'}=$qName;
	    $blastRes[$hitscount]->{'subject'}=$sName;
	    $blastRes[$hitscount]->{'description'}=$hitBN->description();
	    $blastRes[$hitscount]->{'strand'}=$genestrand;
		$blastRes[$hitscount]->{'queryS'}=$queryS;
		$blastRes[$hitscount]->{'queryE'}=$queryE;
	    $blastRes[$hitscount]->{'chrS'}=$chrS;
	    $blastRes[$hitscount]->{'chrE'}=$chrE;
	    $blastRes[$hitscount]->{'evalue'}=$evalue;
	    $blastRes[$hitscount]->{'identity'}=$identity;
	    $blastRes[$hitscount]->{'score'}=$score;
	    $blastRes[$hitscount]->{'qlen'}=$resBN->query_length();
	    $blastRes[$hitscount]->{'slen'}=$hitBN->length();
	    $hitscount++;
	}
	return @blastRes;
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


sub Usage(){
  print << "    Usage";

	Usage:  $0 (version $version)

	<options>
		-i     Input Fasta sequences File
		-d     NCBI Blast database (default: nr)
		-p     Blast Program
		-o     Output File
		-r     Number of the Blast Hits

    Usage
	exit(0);
}
