##this program is used to remove linker/primer sequence from solexa result

#!/usr/bin/perl
use strict;
use Getopt::Long;

#----------------user option--------------

my %opts;

GetOptions(\%opts,"i:s","o:s","v:s","l:s","help");

if(!defined($opts{i} && !defined($opts{v})))
{
    Usage();
    exit();
}

my $input = $opts{i};
my $output = defined $opts{o} ? $opts{o} : $input.".clean";
my $linker = $opts{v};
##my $minilength = defined $opts{m} ? $opts{m} : 11;
##my $gap = defined $opts{g} ? $opts{g} : 1;
my $log = defined $opts{l} ? $opts{l} : $input.".log";
#------------------------------------------

#--------------variable list---------------
my %linkers;
my %rev_linkers;
#------------------------------------------

#----------------read linker file----------
open (linker,$linker) || die "Can't open linker file: $linker, please check again!\n";
my $seqname;
my $seq;
while(<linker>){
    chomp();
    if (/>(\S+)/){
        $seqname=$1;
    }
    else{
        $linkers{$seqname}=$_;
    }
}
&revseq;
#------------------------------------------

#-------------read file--------------------
open (input,"$input") || die "Can't open file: $input, please check again!\n";
open (logfile,">$log");
open (outfile,">$output");

my $seqcount;
my $line;
my $site;
my $seqname;
my $isoutput=1;
my @oldtime=localtime(time);
my $oldsecond=$oldtime[2]*3600+$oldtime[1]*60+$oldtime[0];

while(<input>){
    chomp();
    $line++;
    if($line==1){
        $seqname=$_;
	$isoutput=1;
        #print outfile "$_\n";
        my @now = localtime(time);
        my $nowsecond=$now[2]*3600+$now[1]*60+$now[0];
        if($nowsecond-$oldsecond>5){
            print "$seqcount sequences filtered\n";
            $oldsecond=$nowsecond;
        }
    }
    elsif($line==2){
        $site=&match($_);
        if($site eq "null"){
            print outfile "$seqname\n$_\n";
        }
        else{
            my @sites=split(/,/,$site);
            if($sites[1] eq "for"){
                my $seq=substr($_,$sites[3]);
		if (length($seq)>=30){
                	print outfile "$seqname\n$seq\n";
		}
		else{
			$isoutput=0;
		}
            }
            elsif($sites[1] eq "rev"){
                my $seq=substr($_,0,$sites[2]);
		if (length($seq)>=30){
                        print outfile "$seqname\n$seq\n";
                }
		else{
			$isoutput=0;
		}
            }
            my $start=$sites[2]+1;
            my $end=$sites[3]+1;
            print logfile "$seqname\t$sites[0]\t$start\t$end\n";
        }    
    }
    elsif($line==4){
        if($site eq "null"){
            print outfile "$_\n";
        }
        else{
		if($isoutput==1){
            		my @sites=split(/,/,$site);
            		if($sites[1] eq "for"){
                		my $seq=substr($_,$sites[3]);
                		print outfile "$seq\n";
            		}
            		elsif($sites[1] eq "rev"){
                		my $seq=substr($_,0,$sites[2]);
                		print outfile "$seq\n";
            		}
		}
        }    
    }
    else{
	if($isoutput==1){
        	print outfile "$_\n";
	}
    }
    if($line==4){
        $line=0;
        $site="";
        $seqcount++;
    }
}
print "total $seqcount sequences filtered\n";
#############################################

#----------------match linker---------------
## only match 5' or 3' region
## if lack of 5' region, must start with middle region
## if lack of 3' region, must end with middle region
## N will match all

sub match{
    my $seq=$_;
    foreach my $key(keys %linkers){
        if($seq=~/$linkers{$key}/){
            my $site="$key,for,$-[0],$+[0]";
            return $site;
        }
    }
    foreach my $key(keys %rev_linkers){
        if($seq=~/$rev_linkers{$key}/){
            my $site="$key,rev,$-[0],$+[0]";
            return $site;
        }
    }
    return "null";
}
###########################################

#----------------get reverse linker---------

sub revseq{
    foreach my $key (keys %linkers){
        my $rev_seq=reverse($linkers{$key});
        $rev_seq=~ tr/[atgc]/[tacg]/;
        $rev_seq=~ tr/[ATGC]/[TACG]/;
        my $seqn=$key."_rev";
        $rev_linkers{$seqn}=$rev_seq;
    }
}

#----------------------------------
sub Usage{
    print "Usage:Program -i sequence(fastq) -v primers(fasta) [-l logfile] [-o output file]\n";
}
