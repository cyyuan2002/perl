use strict;
use warnings;

# <size> <read_len> <error_rate> <kmer_size> <coverage> <repeat_rate> <hete_rate> <outDir>

my($x,$e,$r,$h) = (0,0,0,0);
my $program = "/ifshk1/BC_gag/Group/ZHANDL/Bin/KmerAnalysis/test.pl";
my $dir = "/ifshk1/BC_gag/Group/ZHANDL/Temp/Kmer";

for($x=53;$x<=53;$x+=1){# coverage
	for($e=0.005;$e<=0.006;$e+=0.001){# error rate of reads
		for($r=0.8;$r<=0.9;$r+=0.1){# repeat rate
			for($h=0.018;$h<=0.020;$h+=0.001){
				my $E=$e*1000;
				my $R=$r*100;
				my $H=$h*1000;
				print "perl $program 100000 100 $e 17 $x $r $h $dir/HE${E}R${R}H${H}X${x}\n";
			}
		}
	}
}
