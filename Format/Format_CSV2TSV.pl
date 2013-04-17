#/usr/bin/perl
use strict;

my($parFileName) =shift;
open(my $fh_infile,$parFileName) || die "Can't open file:$parFileName\n";
my @content=<$fh_infile>;
close $fh_infile;

my $length=scalar(@content);
my @lastline=split(/\,/,$content[$length-1]);
my $wides=scalar(@lastline);
for(my $i=0;$i<@content;$i++){
	my $info=$content[$i];
	$info=~s/\"//g;
	my @lines=split(/,/,$info);
	next if(scalar(@lines)<$wides);
	my $lineinfo=join("\t",@lines);
	print "$lineinfo";
}
exit(0);