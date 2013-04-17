#!/usr/bin/perl
use Getopt::Long;
use strict;

my %opts;
my $version;
my $jellyfish;

GetOptions(\%opts,"i=s","k:i","c:i","s:i","t:i","r:i","g:i","p:s");
if(!defined $opts{i}){
	&Usage();
}

my $inputfile=$opts{i};
my $kmerlength=(defined $opts{k}) ? $opts{k} : 31;
my $hashsize=(defined $opts{s}) ? $opts{s} : 10000000;
my $threads=(defined $opts{t}) ? $opts{t} : 4;
my $jellyfish=(defined $opts{p}) ?  $opts{p} : "";
my $counterlength=(defined $opts{c}) ? $opts{c} : 5;
my $graphoutput=(defined $opts{g}) ? $opts{g} : 1;
my $removetemp=(defined $opts{r}) ? $opts{r} : 1;

my @filename=split(/\./,$inputfile);
my $filetype=pop(@filename);
my $iszipped=0;
if($filetype eq "gz"){
	$iszipped=1;
}

if($jellyfish eq ""){
	$jellyfish=`which jellyfish`;
	if($jellyfish eq ""){
		print stderr "Please tell me where you have installed jellyfish!";
		exit(0);
	}
	$jellyfish=~s/\n//g;
}

my $tempfile;

print STDERR "\nRunning jellyfish...\n\n";
if($iszipped==1){
	$tempfile=time;
	`gunzip -c $inputfile > $tempfile`;
	`$jellyfish count -m $kmerlength -o $inputfile.jf -c $counterlength -C -s $hashsize -t $threads $tempfile`;
	`rm $tempfile`;
}
else{
	`$jellyfish count -m $kmerlength -o $inputfile.jf -c $counterlength -C -s $hashsize -t $threads $inputfile`;
}

print STDERR "Merging the results...\n\n";
`$jellyfish merge -o $inputfile.jf $inputfile.jf_*`;
print STDERR "Stats the results...\n\n";
`$jellyfish stats $inputfile.jf > $inputfile.stat`;
print STDERR "Output graphic distribution...\n\n";
`$jellyfish histo $inputfile.jf > $inputfile.hist`;
if($graphoutput>=1){
	open (my $scriptfile,">$inputfile.gnu");
	print $scriptfile "set terminal png nocrop size 1280,1024\nset format y \'10\^\(\%\.0f\)\'\nset title \'Kmer Distrbution of $inputfile\'\n";
	print $scriptfile "set xlabel \'multiplicity\'\nset ylabel \'Number of distict K-mers with given multiplicity\'\nset xrange[2:100]\n";
	print $scriptfile "plot \'$inputfile.hist\' using 1\:2 with line title \"kmer_$kmerlength\"\;\n";
	close $scriptfile;
}
`gnuplot < $inputfile.gnu > $inputfile.hist.png`;

print "Remove the temporary files...\n\n";
if($removetemp==1){
	`rm $inputfile.jf`;
	`rm -rf $inputfile.jf_*`;
	`rm $inputfile.gnu`;
}
print "All jobs finished!...\n";

sub Usage(){
	print << "	Usage";
	
	Usage: $0 version $version <options>
	
	-i      Name of the fastq file, must be given
	
	-k      K-mer length (default:31)
	
	-c      counter length (default:5)
	         
	-s      Hash size (default:10000000)
	
	-t      Threads use (default:4)
	
	-r      Remove all the Temporary file (default:1 (Yes))
	
	-p      Path of Jellyfish (example:/Users/Alvin/bin/jellyfish/bin/jellyfish)
	
	-g      Graph output (default:1)
	
	Usage
	
	exit(0);
}