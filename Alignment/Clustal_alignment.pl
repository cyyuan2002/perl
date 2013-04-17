#!/usr/bin/perl -w 

use Data::Dumper;
use strict;
use Getopt::Long;

my $Version="1.0";
my $Date="2004-9-16";
my $update="2004-9-16";
my $Author="Jiang Hui-feng";
my $function="calculation indel from clustalw alignment and exacting alignment sequences cut off all gaps for ka/ks estimation";
my $Contact="jianghf\@genomics.org.cn";

my %opts;
GetOptions(\%opts,"i:s","indel:s","o:s","help");#"seq_in2_ja2:s","seq_in1_in2:s","seq_ja1_ja2:s","help");

if ((!defined $opts{i})|| (!defined $opts{indel}) || (!defined $opts{o})) {
	Usage();
}
our (%hash_indel,%hash_seq);
my %hash_aln={};
my @split_by_codon_seq=();
my @sequence=();
my $start_codon='';
my @name=();
my $id='';
my $num=0;
			#	

	get_align_from_clustalw ($opts{i},\%hash_aln);

	foreach (keys %hash_aln) {
		if ($hash_aln{$_} ne '') {
			if ($_=~/^Y(\S+)/) {
				my $len=length($hash_aln{$_});
				my $gap=$hash_aln{$_};
				$gap=~s/^\-+//g;
				my $gap_len=length($gap);
				$start_codon=$len-$gap_len+1;
			}
			my $sequence=$hash_aln{$_};
			push (@name,$_);
			$num=length ($hash_aln{$_});
	#		print "$_\t$hash_aln{$_}\n";
			push (@sequence,$sequence);
		}
	}

	@split_by_codon_seq=&split_seq (@sequence);
	count_indel (@split_by_codon_seq);
	
	open(INdel,">>$opts{indel}");
	foreach  (keys %hash_indel) {
		print INdel "$_\n$hash_indel{$_}\n";
	}

	my $id_1=0;
	my $seq_1='';
	my $id_2=0;
	my $seq_2='';
	my $id_ja_1=0;
	my $seq_ja_1='';
	my $id_ja_2=0;
	my $seq_ja_2='';
	my $id_ho_1=0;
	my $seq_ho_1='';

	foreach  (keys %hash_seq) {
		#print "$_\n$hash_seq{$_}\n";
		my $len=length ($hash_seq{$_});
	#	if (/(\S+)_indica_chr11/) {
		if (/S(\S+)/) {
			$id_1=$_;
			$seq_1=$hash_seq{$_};
		}
		elsif (/Y(\S+)/) {
			$id_2=$_;
			$seq_2=$hash_seq{$_};
		}
	#	elsif (/(\S+)_japonica_chr11/) {
		elsif (/(\S+)_japonica/) {
			$id_ja_1=$_;
			$seq_ja_1=$hash_seq{$_};
		}
		elsif (/(\S+)_japonica_chr12/) {
			$id_ja_2=$_;
			$seq_ja_2=$hash_seq{$_};
		}
		elsif (/(\S+)_homolog/) {
			$id_ho_1=$_;
			$seq_ho_1=$hash_seq{$_};
		}

	}
#	open (Seq_in2_1,">>$opts{seq_in2_1}");
	open (OUT,">>$opts{o}");
#	open (Seq_in1_in2,">>$opts{seq_in1_in2}");
#	open (Seq_ja1_ja2,">>$opts{seq_ja1_ja2}");
#	print Seq_in2_ja1 "$id_ho_1\_$id_ja_2\n$seq_ho_1\n$seq_ja_2\n\n";
	print OUT "$id_1\_$id_2\n$seq_1\n$seq_2\n\n";
#	print Seq_in1_in2 "$id_in_1\_$id_in_2\n$seq_in_1\n$seq_in_2\n\n";
#	print Seq_ja1_ja2 "$id_ja_1\_$id_ja_2\n$seq_ja_1\n$seq_ja_2\n\n";

sub get_align_from_clustalw () {
	my ($file,$r_return)=@_;
	open (I,$file)||die "can't open *.aln";
	while (<I>) {
		chomp;
			if (/(\S+)\s+(\S+)/) {

				if ($1=~/\*/ || $1=~/CLUSTAL/) {
					next;
				}
				if (exists $$r_return{$1}) {
					$id=$1;
					$$r_return{$1}.=$2;
				}
				else {
					$$r_return{$1}=$2;
				}
			}
	}
	close I;
}

sub split_seq () {
	my (@seq)=@_;
	my $line='';
	my @codon_seq=();

	for (my $j=$start_codon;$j<=$num ;$j+=3) {
		for (my $i=0;$i<@seq ;$i++) {
			my $tri=substr ($seq[$i],$j-1,3);
			if (length($tri)==3) {
				$line=$line."\t".$tri;
			}
		}
		push (@codon_seq,$line);
		$line='';
	}
	return (@codon_seq);
}

sub count_indel () {
	my (@codon)=@_;
	my $gap=0;
	for (my $i=0;$i<@codon ;$i++) {
#		print "$codon[$i]\n";
		my @info=split(/\s+/,$codon[$i]);
		my $g=0;
		for (my $j=1;$j<=@name ;$j++) {
		#	print "$name[$j-1]\n";
			if ($info[$j]=~/-/) {
				$g=1;
				$info[$j]=~s/-//g;
				$gap=3-length($info[$j]);
				my $line="$i\t$gap\n";
				if (exists $hash_indel{$name[$j-1]}) {
					$hash_indel{$name[$j-1]}.=$line;
				}
				else {
					$hash_indel{$name[$j-1]}=$line;
				}
			}
			elsif ($info[$j]=~/TAA/i || $info[$j]=~/TGA/i) {
				$g=1;
			}
		}
		if ($g==0) {
			for (my $r=1;$r<=@name ;$r++) {
				if (exists $hash_seq{$name[$r-1]}) {
					$hash_seq{$name[$r-1]}.=$info[$r];
	#				print "$name[$r-1]\t$hash_seq{$name[$r-1]}\n";
				}
				else {
					$hash_seq{$name[$r-1]}=$info[$r];
				}
			}
			$g=0;
		}
	}
}





sub Usage #help subprogram
{
    print << "    Usage";
    
	Version: $Version
	Date   : $Date
	Update : $update
	Author : $Author
	Contact: $Contact
	
	Function description :
	
		$function

	Usage: $0 <options>

		-i             input of clustalw result file , must be given (String)
		 
			
		-indel         output of indel count, must be given (String)
		

		-o             output of alignment sequence filter out gaps


		-help          Show help , have a choice

    Usage

	exit(0);
};		




