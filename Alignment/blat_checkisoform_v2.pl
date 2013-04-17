#!/usr/bin/perl 
use strict;

my $file=shift;
open(file,$file) || die "Can't open file $file\n";
my @blat=<file>;
close file;

my %isoforms;
my %checked;

for(my $i=0;$i<@blat;$i++){
    my @lines=split(/\t/,$blat[$i]);
    my $tname=$lines[9];
    my $tlength=$lines[0]+$lines[1];
    my $tstrand=$lines[8];
    my $tchroms=$lines[13];
    my $tstart=$lines[15];
    my $tend=$lines[16];
    
    my @texonlengths=split(/,/,$lines[18]);
    my @tchromstarts=split(/,/,$lines[20]);
    my @texons;
    for(my $k=0;$k<@texonlengths;$k++){
        my $tends=$tchromstarts[$k]+$texonlengths[$k]-1;
        my $texon=$tchromstarts[$k].":".$tends;
        push(@texons,$texon);
    }
    
    for(my $j=$i+1;$j<@blat;$j++){
        my @line=split(/\t/,$blat[$j]);

        if(exists $checked{$line[9]}){
            next;
        }

        if($line[13] eq $tchroms){
            if($line[8] eq $tstrand){
                if(($tstart<=$line[15] && $tend>$line[15]) || $tstart<=$line[16] && $tend>$line[16]){
                    
                    my $slength=$line[0]+$line[1];
                    my $sname=$line[9];
                    my @sexonlengths=split(/,/,$line[18]);
                    my @schromstarts=split(/,/,$line[20]);
                    my @sexons;
                    for(my $k=0;$k<@sexonlengths;$k++){
                        my $sends=$schromstarts[$k]+$sexonlengths[$k]-1;
                        my $sexon=$schromstarts[$k].":".$sends;
                        push(@sexons,$sexon);
                    }
                    
                    my $exonmatch=0;
                    my $lengthmatch=0;
                    my $exonnum=0;
                    my $querylength=0;
                                        
                    for(my $m=0;$m<@texons;$m++){
                        my @texonregion=split(/\:/,$texons[$m]);
                        
                        for(my $n=$exonnum;$n<@sexons;$n++){
                            my @sexonregion=split(/\:/,$sexons[$n]);
                            if($texonregion[0]<=$sexonregion[0] && $texonregion[1]>=$sexonregion[0]){
                                if($texonregion[1]<=$sexonregion[1]){
                                    my $overlength=$texonregion[1]-$sexonregion[0]+1;
                                    $lengthmatch+=$overlength;
                                }
                                else{
                                    my $overlength=$sexonregion[1]-$sexonregion[0]+1;
                                    $lengthmatch+=$overlength;
                                }
                                $exonmatch++;
                                $exonnum=$n;
                                last;
                            }
                            elsif($texonregion[0]<=$sexonregion[1] && $texonregion[1]>=$sexonregion[1]){
                                if($texonregion[1]<=$sexonregion[1]){
                                    my $overlength=$texonregion[1]-$texonregion[0];
                                    $lengthmatch+=$overlength;
                                }
                                else{
                                    my $overlength=$sexonregion[1]-$texonregion[0];
                                    $lengthmatch+=$overlength;
                                }
                                $exonmatch++;
                                $exonnum=$i;
                                last;
                            }
                        }
                        
                    }              
                    
                    my $submatchpercent=$lengthmatch/$tlength;
                    my $querymatchpercent=$lengthmatch/$slength;
                
                    if($exonmatch>=3 || $submatchpercent>=0.5 || $querymatchpercent>=0.5){
                        if(!exists $isoforms{$tname}){
                            $isoforms{$tname}=$tname."\t".$sname;
                            $checked{$tname}=1;
                            $checked{$sname}=1;
                        }
                        else{
                            $isoforms{$tname}=$isoforms{$tname}."\t".$sname;
                            $checked{$sname}=1;
                        }
                    } 
                }
            }
        }
    }
}

my $outfile=$file.".isoform2";
open(outfile,">$outfile");

for my $key (keys %isoforms){
    print outfile "$isoforms{$key}\n";
}
close outfile;
