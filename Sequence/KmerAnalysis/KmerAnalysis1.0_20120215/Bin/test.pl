use strict;
use warnings;
use FindBin qw($Bin);

die "perl $0 <size> <read_len> <error_rate> <kmer_size> <coverage> <repeat_rate> <hete_rate> <outDir>
example:\nperl $0 10000000 100 0.01 17 20 0.1 0.003 Test\n" if(@ARGV==0);

my($size,$rd_len,$err_rate,$kmer,$coverage,$repeat,$hete,$outDir) = @ARGV;
&printPara();#debug
mkdir("$outDir") unless(-e $outDir);
chdir($outDir);

# generate sequence with repeat
my $number = int($size*$repeat/300);
my $len = int($size*$repeat);
system("perl $Bin/GenRepeat.pl --size $size --STR_len $len --STR_node $number");

# record the sequence with jellyfish
mkdir("Data");
chdir("Data");
system("echo ../lardo.fa >kmer.lst");
system("perl $Bin/kmer_counter.pl --k $kmer --memory 1 --prefix seq kmer.lst");
system("sh seq.sh"); 
chdir("..");

# generate reads with error and hete rate
system("perl $Bin/ReadGener.pl $err_rate $hete $coverage 500 $rd_len lardo.fa Reads");

# kmer Analysis with jellyfish
mkdir("Kmer") unless(-e "Kmer");
chdir("Kmer");
system("ls ../Reads/*.fa >kmer.lst");	#get read list
system("perl $Bin/kmer_counter.pl --k $kmer --memory 1 --prefix kmer kmer.lst");
system("sh kmer.sh"); 

# clean the data to save space
chdir("..");
system("mv Kmer/kmer.stats Kmer/kmer.histo Data/seq.stats Data/seq.histo . && rm -r Reads Kmer Data lardo.*");

# convert the jellyfish files to the format like "count number percent";
my $total = `grep Distinct kmer.stats`;
chomp $total;
$total = (split(/\s+/,$total))[1];
open FI,"kmer.histo";
open FO,">kmer.xls";
my $line = 0;
while(<FI>){
	chomp;
	my($count,$num) = split;
	my $percent = $num/$total;
	print FO "$count\t$percent\t$num\n";
	$line++;
	last if($line>=500);
}
close FO;
close FI;

#########################################################################################################################################################
############################################### subroutine ##############################################################################################
#########################################################################################################################################################

#print parameters for debug
sub printPara(){
	print STDERR "genome size: $size\nread len:$rd_len\nerror rate:$err_rate\nkmer size:$kmer\ncoverage:$coverage\nrepeat rate:$repeat\nhete rate:$hete\noutDir:$outDir\n";
}