#!/usr/bin/perl

use Getopt::Long;  
use strict;

my %opts;
GetOptions(\%opts,"i:s","o:s","d:s");


if ((!defined $opts{i})|| (!defined $opts{o})||(!defined $opts{d})) {
	Usage();
}

my $Input;
my $Output;
my $Database;


$Input = $opts{i};
$Output = $opts{o};
$Database = $opts{d};

open(Input,$Input) || die "can't open file $Input\n";
open(Output,">$Output") || die "can't open output file $Output\n";
print "Running...\n";
while(<Input>){
	chomp();
	my $title=$_;
	getcontent($title);
}

sub getcontent(){
	open(Database,$Database) || "can't open database $Database\n";
	my $linename=@_[0];
	my $ismatch=0;
	#print $linename."\n";
	while(<Database>){
		chomp();
		if(/^>(\S+)/){
			if($linename eq $1){
				$ismatch=1;
				print Output "\>$linename\n";
				print "Out put sequence $linename\n";
			}
			elsif($ismatch==1){ 
				last;
			}			
		}
		elsif(/(\w+)/){
			if($ismatch==1){
				print Output $1."\n";
			}
		}
	}
	close Database;
}

close Input;
close Output;

sub Usage() #help subprogram
{
    print << "    Usage";
    
	Usage: $0 <options>

		-i     input of fasta file , must be given (String)
					
		-o     output of fasta , must be given (String)

		-d     file name of the database
		
    Usage

	exit(0);
};	
