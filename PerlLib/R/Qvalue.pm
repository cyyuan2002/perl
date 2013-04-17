package R::Qvalue;

use strict;
use warnings;
use Carp;

sub new {
	my $packagename = shift;
	my (@pvalues) = @_;
	
	my $self = { 
		pvalues => \@pvalues,
		qvalues => undef,
	};

	bless ($self, $packagename);

	return($self);
}

####
sub get_Qvalues {
	my $self = shift;

	my $tmp_pvalues = "/tmp/pval.$$.txt";
	open (my $ofh, ">$tmp_pvalues") or die "Error, cannot write to $tmp_pvalues";
	print $ofh join("\n", @{$self->{pvalues}}) . "\n";
	close $ofh;
	
	my $tmp_Rscript = "/tmp/qval.$$.Rscript";
	my $tmp_Qvalues = "/tmp/qval.$$.txt";
	
	open ($ofh, ">$tmp_Rscript") or die "Error, cannot write to $tmp_Rscript";
	print $ofh "pvalues = read.table(file=\"$tmp_pvalues\", header=F)\n";
	print $ofh "library(qvalue)\n";
	print $ofh "qobj = qvalue(pvalues[,1])\n";
	print $ofh "sink(file=\"$tmp_Qvalues\")\n";
	print $ofh "for (i in 1:length(qobj\$qvalues)) {\n";
	print $ofh "   cat (paste(qobj\$qvalues[i], '\\n'))\n";
	print $ofh "}\n";
	print $ofh "sink()\n";
	close $ofh;


	my $cmd = "R --vanilla -q --slave < $tmp_Rscript";

	my $ret = system($cmd);
	if ($ret) {
		confess "Error, cmd: $cmd died with ret $ret";
	}

	my @qvalues = `cat $tmp_Qvalues`;
	chomp @qvalues;

	#unlink ($tmp_pvalues, $tmp_Rscript, $tmp_Qvalues);

	return(@qvalues);
}

#######


sub __test {

	my @pvalues;
	for (1..1000) {
		if (rand(3) > 2) {
			my $pvalue = rand(1)/10000;
			push (@pvalues, $pvalue);
		}
		else {
			push (@pvalues, rand(1));
		}
	}
	
	my $test_obj = new R::Qvalue(@pvalues);

	my @qvalues = $test_obj->get_Qvalues();
	
	print "Results:\n" . join("\n", @qvalues) . "\n";

	exit(0);
}


		
		


1; #EOM
