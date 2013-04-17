#!/usr/bin/perl
use strict;

my $filename=shift;

open (filein,$filename) || die "Can't open file $filename \n";
my @annos=<filein>;
close filein;

my @checked;
my @exons;
my $rnaid;
my $isrecord;
my $rnaregion;
my $strand;
my $cdslength;
my %isoforms;
my $chroms;

for(my $i=0;$i<@annos;$i++){
    my @lines=split(/\t/,$annos[$i]);
    if($lines[2] eq "mRNA"){
        if($rnaid ne "" && $isrecord==1){
            &check_trans($i);
            @exons=();
            $rnaid="";
            $strand="";
            $chroms="";
            $cdslength=0;
        }
        
        my @rnaids=split(/\;/,$lines[8]);
        my @rnanum=split(/\s/,$rnaids[0]);
        $rnaid=$rnanum[1];
        $chroms=$lines[0];
        $rnaregion=$lines[3].":".$lines[4];
        $strand=$lines[6];
        
        $isrecord=1;
        
        foreach my $key(@checked){
            if($rnaid eq $key){
                $isrecord=0;
                last;
            }
        }
        
    }
    else{
        if($isrecord==1){
            my $exon=$lines[3].":".$lines[4];
            push (@exons,$exon);
            $cdslength=$cdslength+($lines[4]-$lines[3]+1);
        }
    }
}

my $outfile=$filename.".isoform";
open(outfile,">$outfile");
foreach my $key (keys %isoforms){
    print outfile "$isoforms{$key}\n";
}




sub check_trans{
    my $startline=shift;
    my $checkdetail=0;
    my $exonmatch=0;
    my $lengthmatch=0;
    my $exonnum=0;
    my $querylength=0;
    my $queryid;
    
    
    for(my $j=$startline;$j<@annos;$j++){
        my @lines=split(/\t/,$annos[$j]);
        if($lines[2] eq "mRNA"){
            
            if($checkdetail==1){
                my $submatchpercent=$lengthmatch/$cdslength;
                my $querymatchpercent=$lengthmatch/$querylength;
                
                if($exonmatch>=3 || $submatchpercent>=0.5 || $querymatchpercent>=0.5){
                    if(!exists $isoforms{$rnaid}){
                        $isoforms{$rnaid}=$rnaid."\t".$queryid;
                        push (@checked,$rnaid);
                        push (@checked,$queryid);
                        #print "$rnaid\t$queryid\n";
                    }
                    else{
                        $isoforms{$rnaid}=$isoforms{$rnaid}."\t".$queryid;
                        push (@checked,$queryid);
                        #print "$rnaid\t$queryid\n";
                    }
                }
            }
            my @rnaids=split(/\;/,$lines[8]);
            my @rnanum=split(/\s/,$rnaids[0]);
            $queryid=$rnanum[1];
            $checkdetail=0;
            
            if($lines[0] eq $chroms){
                if($lines[6] eq $strand){
                    my @queryregion=split(/\:/,$rnaregion);
                    if(($queryregion[0]<=$lines[3] && $queryregion[1]>$lines[3]) || $queryregion[0]<=$lines[4] && $queryregion[1]>$lines[4]){
                        $checkdetail=1;
                    }
                }
            }
            
            $exonmatch=0;
            $lengthmatch=0;
            $exonnum=0;
            $querylength=0;
        }
        else{
            if($checkdetail==1){
                $querylength=$querylength+($lines[4]-$lines[3]+1);
                for(my $i=$exonnum;$i<@exons;$i++){
                    my @exonregion=split(/\:/,$exons[$i]);
                    if($exonregion[0]<=$lines[3] && $exonregion[1]>$lines[3]){
                        if($exonregion[1]<=$lines[4]){
                            my $overlength=$exonregion[1]-$lines[3];
                            $lengthmatch+=$overlength;
                        }
                        else{
                            my $overlength=$lines[4]-$lines[3];
                            $lengthmatch+=$overlength;
                        }
                        $exonmatch++;
                        $exonnum=$i;
                        last;
                    }
                    elsif($exonregion[0]<=$lines[4] && $exonregion[1]>=$lines[4]){
                        if($exonregion[1]<=$lines[4]){
                            my $overlength=$exonregion[1]-$exonregion[0];
                            $lengthmatch+=$overlength;
                        }
                        else{
                            my $overlength=$lines[4]-$exonregion[0];
                            $lengthmatch+=$overlength;
                        }
                        $exonmatch++;
                        $exonnum=$i;
                        last;
                    }
                }
            }
        }
    }
    
    if($checkdetail==1){
        my $submatchpercent=$lengthmatch/$cdslength;
        my $querymatchpercent=$lengthmatch/$querylength;
                
        if($exonmatch>=3 || $submatchpercent>=0.5 || $querymatchpercent>=0.5){
            if(!exists $isoforms{$rnaid}){
                $isoforms{$rnaid}=$rnaid."\t".$queryid;
                #print "$rnaid\t$queryid\n";
            }
            else{
                $isoforms{$rnaid}=$isoforms{$rnaid}."\t".$queryid;
                #print "$rnaid\t$queryid\n";
            }
        }
    }
}


