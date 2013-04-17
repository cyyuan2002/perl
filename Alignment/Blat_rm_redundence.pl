##this program is used to remove the redundance data from the predicted transcript sequence
##if two hit map on the same scaffold and share a overlap longer than 50% of the total length, the shorter one will seem as redundence sequence.


#!/usr/bin/perl -w
use strict;
use warnings;

#match	mis- 	rep. 	N's	Q gap	Q gap	T gap	T gap	strand	Q        	Q   	Q    	Q  	T        	T   	T    	T  	block	blockSizes 	qStarts	 tStarts
#     	match	match	   	count	bases	count	bases	      	name     	size	start	end	name     	size	start	end	count
#819	48	0	0	1	9	7	6105	-	jgi|Brafl1|70551|fgenesh2_pg.scaffold_29000001	876	0	876	scaffold_3	6048608	3018922	3025894	9	93,9,53,145,162,165,53,4,183,	0,102,111,164,309,471,636,689,693,	3018922,3019015,3019036,3022962,3023492,3024164,3024889,3025028,3025711,

my ($filein,$fileout1,$fileout2)=@ARGV;
open(filein,"$filein");
########this step is to found the best hit of each query###########
my $temp1=$filein.".tmp1";
open(temp1,">$temp1");
my $QueryN="";
my $Iden=0;
my $lineInfo=0;
while(<filein>){
	chomp();
	my @lines=split(/\t/,$_);
	my $score=$lines[0];
	my @queryname=split(/\|/,$lines[9]);
	my $qname=$queryname[2];	
	my $qlength=$lines[10];
	my $identity=$score/$qlength;
	if ($qname eq $QueryN){
		if ($identity>$Iden){
			$lineInfo=$_;
			$Iden=$identity;
			#print "$qname,$identity\n";
		}
	}
	else{
		if (!($QueryN eq "")){
			print temp1 "$lineInfo\n";	
		}
		$lineInfo=$_;
		$QueryN=$qname;
		$Iden=$identity;
	}
}
close temp1;
##########################
#generate a temp file for sort by scaffold
my $temp2=$filein.".tmp2";
my $command="sort -k14,14 $temp1 > $temp2";
system($command);
###--
open (temp2,"$temp2");

my (@qName,@qIde,@qLen,@str,@sStart,@sEnd,@sLen);
my $subName="";
open (fileout1,">$fileout1");
open (fileout2,">$fileout2");
while(<temp2>){
	chomp();
	my @lines=split(/\t/,$_);
	my $score=$lines[0];
	my $strand=$lines[8];
	my @queryname=split(/\|/,$lines[9]);
	my $qname=$queryname[2];	
	my $qlength=$lines[10];
	my $subname=$lines[13];
	my $substart=$lines[15];
	my $subend=$lines[16];
	my $sublength=$subend-$substart;
	#print "$score,$qlength\n";
	my $identity=$score/$qlength;
	if ($subName eq $subname){
		my ($state,$id)=&checksite($substart,$subend,$strand,$sublength);
		# result=0 need to change data;result=1 jump to next;result=2 add new data; 
		if ($state==0){
			$qName[$id]=$qname;
			$qIde[$id]=$identity;
			$qLen[$id]=$qlength;
			$str[$id]=$strand;
			$sStart[$id]=$substart;
			$sEnd[$id]=$subend;
			$sLen[$id]=$sublength;			
		}
		elsif ($state==1){
			next;
			print fileout2 "$qname\n";
		}
		elsif ($state==2){
			push(@qName,$qname);
			push(@qIde,$identity);
			push(@qLen,$qlength);
			push(@str,$strand);
			push(@sStart,$substart);
			push(@sEnd,$subend);
			push(@sLen,$sublength);	
		}
	}
	else{
		if (!($subName eq "")){
			print "Output $subName\n";
			for(my $i=0;$i<@qName;$i++){
				print fileout1 "$qName[$i]\t$qLen[$i]\t$qIde[$i]\t$str[$i]\t$subName\t$sLen[$i]\t$sStart[$i]\t$sEnd[$i]\n";
			}
		}
		@qName=();
		@qIde=();
		@qLen=();
		@str=();
		@sStart=();
		@sEnd=();
		@sLen=();
		$subName=$subname;
		push(@qName,$qname);
		push(@qIde,$identity);
		push(@qLen,$qlength);
		push(@str,$strand);
		push(@sStart,$substart);
		push(@sEnd,$subend);
		push(@sLen,$sublength);
	}
}

close temp2;
unlink $temp1;
unlink $temp2;

sub checksite(){
	my ($substart,$subend,$strand,$sublength)=@_;
	for(my $i=0;$i<@sStart;$i++){
		if ($str[$i] eq $strand){
			if($substart<=$sStart[$i] && $subend>=$sStart[$i]){
				if ($subend>=$sEnd[$i]){
					return (0,$i);
				}
				else {
					my $totalLen=$sEnd[$i]-$substart;
					my $overlap=$subend-$sStart[$i];
					my $overper=$overlap/$totalLen;
					if ($overper>0.5){
						if ($sublength>$sLen[$i]){
							return (0,$i);
						}
						else {return (1,0);}
					}
				}
			}
			elsif($substart<=$sEnd[$i] && $subend>=$sEnd[$i]){
				my $totalLen=$subend-$sStart[$i];
				my $overlap=$sEnd[$i]-$substart;
				my $overper=$overlap/$totalLen;
				if ($overper>0.5){
					if ($sublength>$sLen[$i]){
						return (0,$i);
					}
					else {return (1,0);}
				}
			}
			elsif($sStart[$i]<=$substart && $sEnd[$i]>=$subend){
				return (1,0);
			}
		}
	}
	return (2,0);
}