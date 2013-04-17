use strict;
#---------------------------------------get the distribution of error----------------------------------------#
sub error_distribution{		
	(my $seq_len,my $rd_len,my $er_rate,my $rd_num,my $rd_info)=@_;
	print STDERR "error rate: $er_rate\n";
	my (%hash,$total);
	my $basic_err=0.001;
	my $total_err=($er_rate-$basic_err)*$rd_len;
	for (my $i=1;$i<=$rd_len;$i++) {
		$hash{$i}=0.00001*$i**3;
		$total+=$hash{$i};
	}
	for (my $circle=1;$circle<=$rd_len;$circle++) {
		$hash{$circle}=$basic_err+$hash{$circle}/$total*$total_err;
		$hash{$circle}=$hash{$circle}*$rd_num*2;
		$hash{$circle}=int $hash{$circle};
	}
	
	for(my $circle=1;$circle<=$rd_len;$circle++) {
		while ($hash{$circle}>0) {
			my $j=int (rand $rd_num);
			push @{$rd_info->{$j}},$circle;
			$hash{$circle}--;
		}
	}
	undef %hash;
}

#----------------------------------------------get Sequence from fa file-----------------------------------#
sub getSequence(){
   my $ref;
   while(my $line = <INPUT>){
      if($line =~ />/){
         seek(INPUT,length($line)+1,1);
         last;
      }else{
         chomp($line);
         $line =~ tr/atgc/ATGC/;
         $ref .= $line;
      }
   }
   return $ref;
}

#------------------------------------get possion value-----------------------------------------------------#
sub getPValue{
	my($k,$u) = @_;
	my $value = exp(-1*$u)*($u**$k)/N($k);
}
sub N{
	my($n)=@_;
	if($n<=0){
		return 1;
	}
	my $sum = 1;
	for(my $i=0;$i<$n;$i++){
		$sum *= ($i+1);
	}
	return $sum;
}

#-------------------------simulate heta distribution------------------------------------------------------#
sub genHexaRef{
   my($heta_rate,$ref) = @_;
   print STDERR "heta rate: $heta_rate\n";
   my $ref_len = length(${$ref});
   my $count = 0;
   my %heta;
   my $n = $heta_rate*$ref_len;
   print STDERR "begin to generate HexaRef:times $n\n";
   while($count <= $n){
      my $position = int(rand($ref_len));
      unless(exists $heta{$position}){
         $heta{$position}=1;
         substr(${$ref},$position,1) =~ tr/ATGC/GCAT/;
         $count++;
      }
   }
   print STDERR "HexaRef generated\n";
}

#-------------------------simulate reads------------------------------------------------------------------#
sub genReads{
   my($ref,$h_ref,$error_info,$insertsize,$read_len,$num)=@_;
   my $count = 0;
   my $length = length(${$ref});
   print STDERR "ref_len:$length\n";
   print STDERR "num:$num\n";
   while($count<$num){
      my $pos = rand($length);
      my $isH = (int rand(10))%2;
      my $line;
      #deside which sequence to be process
      if($isH==0){
         $line = substr(${$ref},$pos,$insertsize);
      }else{
         $line = substr(${$h_ref},$pos,$insertsize);
      }
      #get reads from the slice
      if(length($line)<$insertsize){
         next;
      }else{
         $count++;
         my $read1 = substr($line,0,$read_len);
         my $read2 = substr($line,length($line)-$read_len,$read_len);
         $read2 = reverse($read2);
         $read2 =~ tr/ATGC/TACG/;
         #generate error
         if(exists $error_info->{$count}){
            foreach my $j (@{$error_info->{$count}}) {
		my $k=int ((rand 10)%2);
		if ($k) {
		substr($read1,$j-1,1) = tr/ATGC/GCAT/;
		}else{
                  substr($read2,$j-1,1) = tr/ATGC/GCAT/;
		}
            }
         }
         #generate reads
         my $isOrigin = int(rand(10))%2;
         if($isOrigin){
            print OUT1 ">${insertsize}_${read_len}_$count\n${read1}\n";
            print OUT2 ">${insertsize}_${read_len}_$count\n${read2}\n";    
         }else{
            print OUT2 ">${insertsize}_${read_len}_$count\n${read1}\n";
            print OUT1 ">${insertsize}_${read_len}_$count\n${read2}\n";
         }
      }
   }
}
#---------------------------------main----------------------------------------------------------------------------#

die "perl ReadGener.pl <error_rate> <h_rate> <coverage> <insertsize> <read_length> <data.fa> <dir>\n" if(@ARGV==0);
my($error_rate,$h_rate,$coverage,$insertsize,$read_len,$data,$dir)=@ARGV;
mkdir $dir unless(-e $dir);
open INPUT,$data or die "can not open file $data\n";
open OUT1,">$dir/read1.fa";
open OUT2,">$dir/read2.fa";
while(<INPUT>){
   if(/>/){
      my $ref = &getSequence();
      my $length = length $ref;
      my $num = $coverage*$length/2/$read_len;
      my $h_ref = $ref;
      my %error_info;
      &genHexaRef($h_rate,\$h_ref);
      if($error_rate>0){
         &error_distribution($length,$read_len,$error_rate,$num,\%error_info);
      }
      print STDERR "begin to generate reads\n";
      &genReads(\$ref,\$h_ref,\%error_info,$insertsize,$read_len,$num);
      print STDERR "all reads generate\n";
   }
}
close INPUT;
close OUT1;
close OUT2;

