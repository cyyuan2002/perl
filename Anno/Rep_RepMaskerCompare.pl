#!/usr/bin/perl
use strict;

my ($RepeatControl, $RepeatNew)=@ARGV;

my $overlap=0.6;

my $tempfile="Repeat.out.tmp";
open(my $fh_temp,">$tempfile");

open(my $fh_repctl, "$RepeatControl") || die "Can't open file: $RepeatControl\n";
while(<$fh_repctl>){
    $_=&trim($_);
    next if(!/^\d/);
    my @lines=split(/\s+/,$_);
    print $fh_temp join("\t",@lines),"\n";
}
close $fh_repctl;

open(my $fh_repnew, "$RepeatNew") || die "Can't open file: $RepeatNew\n";
while(<$fh_repnew>){
    $_=&trim($_);
    next if(!/^\d/);
    my @lines=split(/\s+/,$_);
    print $fh_temp join("\t",@lines),"\n";
}
close $fh_repnew;
close $fh_temp;

`cat $tempfile | sort -k 5,5 -k 6,6n -k 7,7n > Repeat.out.tmp2`;

my %repeatcounts;
my %repeatsinfo;
my %lastRep;
my %lastDenovoRep;

open(my $fh_temp2,"Repeat.out.tmp2");
while(<$fh_temp2>){
    chomp();
    my @lines=split(/\t/,$_);
    if($lines[10] eq "Denovo"){
        if(exists($lastRep{'id'})){
            if($lastRep{'e'} > $lines[5]){
                if($lastRep{'e'} >= $lines[6]){
                    if(exists($repeatsinfo{$lines[9]})){
                        $repeatsinfo{$lines[9]}.=",$lastRep{'id'}";
                    }
                    else{
                        $repeatsinfo{$lines[9]}=$lastRep{'id'};
                    }
                }
                else{
                    my $overlength=$lastRep{'e'}-$lines[5]+1;
                    my $replength=$lines[6]-$lines[5]+1;
                    my $overpercent=$replength/$overlength;
                    if($overpercent/$replength >= $overlap){
                        if(exists($repeatsinfo{$lines[9]})){
                            $repeatsinfo{$lines[9]}.=",$lastRep{'id'}";
                        }
                        else{
                            $repeatsinfo{$lines[9]}=$lastRep{'id'};
                        }
                    }
                }
            }
        }
        $lastDenovoRep{'chrom'}=$lines[4];
        $lastDenovoRep{'s'}=$lines[5];
        $lastDenovoRep{'e'}=$lines[6];
        $lastDenovoRep{'id'}=$lines[9];
        if(exists($repeatcounts{$lines[9]})){
            $repeatcounts{$lines[9]}++;
        }
        else{
            $repeatcounts{$lines[9]}=1;
        }
    }
    else{
        if(exists($lastDenovoRep{'id'})){
            if($lastDenovoRep{'e'} > $lines[5]){
                my $overlength;
                my $replength=$lastDenovoRep{'e'}-$lastDenovoRep{'s'}+1;
                if($lastDenovoRep{'e'} >= $lines[6]){
                    $overlength=$lines[6]-$lines[5]+1;
                }
                else{
                    $overlength=$lastDenovoRep{'e'}-$lines[5]+1;
                }
                my $overpercent=$replength/$overlength;
                if($overpercent/$replength >= $overlap){
                    if(exists($repeatsinfo{$lines[9]})){
                        $repeatsinfo{$lines[9]}.=",$lastRep{'id'}";
                    }
                    else{
                        $repeatsinfo{$lines[9]}=$lastRep{'id'};
                    }
                }
            }
        }
        $lastRep{'chrom'}=$lines[4];
        $lastRep{'s'}=$lines[5];
        $lastRep{'e'}=$lines[6];
        $lastRep{'id'}=$lines[9];
        if($lines[10] ne "Simple_repeat"){
            if(exists($repeatcounts{$lines[9]})){
                $repeatcounts{$lines[9]}++;
            }
            else{
                $repeatcounts{$lines[9]}=1;
            }
        }
    }
}


foreach my $key(sort {$a<=>$b} keys %repeatsinfo){
    print "$key\t$repeatsinfo{$key}\n";
}

print "##\n";
foreach my $key(sort {$a<=>$b} keys %repeatcounts){
    print "$key\t$repeatcounts{$key}\n";
}

exit(0);

sub trim($)
{
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}