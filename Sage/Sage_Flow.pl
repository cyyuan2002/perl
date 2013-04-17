#This program is used for analyzing SAGE data


#!/usr/bin/perl

use strict;
my $version="1.0 Alvin Chen 2010-12-28";

use Getopt::Long;
my %opts;
#c configure file, i input file, s score cut for tagreads, a R space file, d Database to match,t Tag length, r redundance model
GetOptions(\%opts,"c:s","i=s","s:s","r:i","t:i","d=s","tl:i","a:s","cl:i","f:i","g:i","m:i","rt:s","l:s","td:s","st:i","help"); 
if((!defined $opts{i})||(!defined $opts{d})){
    &Usage();
}
##parameters
my $timestart=time;

my $parSourcefile=$opts{i};
my $parDBfile=$opts{d};
my $parRdatafile=(defined $opts{a}) ? $opts{a} : "$parSourcefile.RData";
my $parConffile=(defined $opts{c}) ? $opts{c} : "config.txt";
my $parScorecut=(defined $opts{s}) ? $opts{s} : 20;
my $parRdmode=(defined $opts{r}) ? $opts{r} : 1;
my $parTaglength=(defined $opts{tl}) ? $opts{tl} : 17;
my $parCountCutoff=(defined $opts{f}) ? $opts{f} : 5;
my $parCleanfile=(defined $opts{e}) ? $opts{e} : 0;
my $parSegnumber=(defined $opts{g}) ? $opts{g} : 30000;
my $parMismatch=(defined $opts{m}) ? $opts{m} : 0;
my $parTopnumber=(defined $opts{t}) ? $opts{t} : 200;
my $parRTemp=(defined $opts{rt}) ? $opts{rt} : "temple.R";
my $parTargetDir=(defined $opts{td}) ? $opts{td} : "$parSourcefile.dir";
my $parLogfile=(defined $opts{l}) ? $opts{l} : "$parSourcefile.$parTargetDir.log";
my $parClearFile=(defined $opts{cl}) ? $opts{cl} :0;
my $parSorttype=(defined $opts{st}) ? $opts{st} : 0;

if($parSorttype==0){
    $parSorttype="p.value";
}
else{
    $parSorttype="logFC";
}

open (logfile,">$parLogfile");

print logfile "\nLoading configuration...\n";
##programs
my $proSageFilter;
my $proTagMapping;
my $proTagCut;
my $proReadsNumb;
my $proR;
my $proJobMan;
my $proExactTest;
my $proTagMapGene;
##SAGE Infors
my @sageFile;
my @sageGroup;
my @sageDesc;
my @sageTitle;
my %sageGroupHash;

##Error collection
my $iserror=0;
##All the created file;
my @fileCreated;
my @fileCompair;
##


&Filecheck($parSourcefile,0);
&Filecheck($parDBfile,0);
&Filecheck($parConffile,0);
&Filecheck($parRTemp,0); 
&Filecheck($parTargetDir,1);

if($iserror==1){
    exit(1);
}

open(conffile,$parConffile);
while(<conffile>){
    #SageFilter	/home/SCE/chenyuan/Sage/bin/Sage_Filter.pl
    chomp();
    my @lines=split(/\t/,$_);
    if($lines[0] eq "SageFilter"){
	$proSageFilter=$lines[1];
    }
    elsif($lines[0] eq "TagMapping"){
	$proTagMapping=$lines[1];
    }
    elsif($lines[0] eq "TagCut"){
	$proTagCut=$lines[1];
    }
    elsif($lines[0] eq "ReadsNumb"){
	$proReadsNumb=$lines[1];
    }
    elsif($lines[0] eq "R"){
	$proR=$lines[1];
    }
    elsif($lines[0] eq "JobManager"){
	$proJobMan=$lines[1];
    }
    elsif($lines[0] eq "ExactTest"){
	$proExactTest=$lines[1];
    }
    elsif($lines[0] eq "TagMapGene"){
	$proTagMapGene=$lines[1];
    }
}
close conffile;

if($proSageFilter eq "" || $proTagMapping eq "" || $proTagCut eq "" || $proReadsNumb eq "" || $proR  eq "" || $proJobMan eq ""){
    print stderr"Configure information error...Please check $parConffile file\n";
    exit(1);
}

&Filecheck($proSageFilter,0);
&Filecheck($proTagMapping,0);
&Filecheck($proTagCut,0);
&Filecheck($proReadsNumb,0);
&Filecheck($proR,0);
&Filecheck($proJobMan,0);
&Filecheck($proExactTest,0);

if($iserror==1){
    exit(1);
}


open(sourcefile,$parSourcefile) || die "Can't open source file:$parSourcefile\n";
my $linecount=0;
while(<sourcefile>){
    ##soucefile format:
    ##filename group description
    chomp();
    my @lines=split(/\t/,$_);
    &Filecheck($lines[0],0);
    push(@sageFile,$lines[0]);
    push(@sageGroup,$lines[1]);
    if(exists($sageGroupHash{$lines[1]})){
	$sageGroupHash{$lines[1]} = "$sageGroupHash{$lines[1]}#$linecount";
    }
    else{
	$sageGroupHash{$lines[1]}=$linecount;
    }
    #print stderr "$lines[1]\t$sageGroupHash{$lines[1]}\n";
    push(@sageDesc,$lines[2]);
    push(@sageTitle,$lines[3]);
    $linecount++;
}
if($iserror==1){
    exit(1);
}

##make directory for the analysis and link the data files to the directory
print logfile "\nCreate Directory $parTargetDir\n";
mkdir ($parTargetDir);
print logfie "\nCreate link files of original data files\n";
symlink ("..\/$parDBfile","$parTargetDir\/$parDBfile");
$parDBfile="$parTargetDir\/$parDBfile";
for(my $i=0;$i<@sageFile;$i++){
    symlink ("..\/$sageFile[$i]","$parTargetDir\/$sageFile[$i]");
    $sageFile[$i]="$parTargetDir\/$sageFile[$i]";
}
$parSourcefile="$parTargetDir\/$parSourcefile";
$parRdatafile="$parTargetDir\/$parRdatafile";

print logfile "\nFilter Sage Files...\n\n";
for(my $i=0;$i<@sageFile;$i++){
    print logfile "perl $proSageFilter -i $sageFile[$i] -s $parScorecut -l $parTaglength\n";
    `perl $proSageFilter -i $sageFile[$i] -s $parScorecut`;
    push(@fileCreated,"$sageFile[$i].flt");
}
print logfile "\nFilter Finished\n";

print logfile "\nCount Tag Numbers...\n\n";
for(my $i=0;$i<@sageFile;$i++){
    print logfile "perl $proReadsNumb $sageFile[$i].flt > $sageFile[$i].flt.no\n";
    `perl $proReadsNumb $sageFile[$i].flt > $sageFile[$i].flt.no`;
    push(@fileCreated,"$sageFile[$i].flt.no");
}
print logfile "\nCount Finished\n";

##cut tags from database file;
print logfile "\nCut Tags from Fasta sequences\n\n";
print logfile "perl $proTagCut -i $parDBfile -l $parTaglength -r $parRdmode\n";
`perl $proTagCut -i $parDBfile -l $parTaglength -r $parRdmode`;
push(@fileCreated,"$parDBfile.tags");

print logfile "\nCreate R Command file\n";
##create sage target file for R
open(Rtarget,">$parSourcefile.target");
push(@fileCreated,"$parSourcefile.target");
print Rtarget "files	group	description\n";
for(my $i=0;$i<@sageFile;$i++){
    print Rtarget "$sageFile[$i].flt.no\t$sageGroup[$i]\t$sageDesc[$i]\n";
}
close Rtarget;

##create R source file and running R;

my $strRCmd=&dealRTemp();
open(Rscript,">$parSourcefile.R");
push (@fileCreated,"$parSourcefile.R");
print Rscript $strRCmd;
close Rscript;

print logfile "\nRunning R to analyze data...\n\n";
print logfile "$proR CMD BATCH --no-save $parSourcefile.R $parSourcefile.Rout\n";
system ("$proR CMD BATCH --no-save $parSourcefile.R $parSourcefile.Rout");
print logfile "\nAnalysis Finished\n\n";

print logfile "\n------------------R Report-------------------\n\n";
open (Rreport,"$parSourcefile.Rout");
while(<Rreport>){
    print logfile $_;
}
close Rreport;
print logfile "\n\n------------------Finished-------------------\n";
unlink ("$parSourcefile.Rout");

print logfile "\nMapping SAGE Tags to Database...\n";

### Mapping tags to genes;
`more $parSourcefile.tab | cut -f1 | sed 's/"//g' > $parSourcefile.tab.tags`;
if($parMismatch==0){
    print logfile "perl ../bin/Sage_TagMapping.pl $parDBfile.tags $parSourcefile.tab.tags 0";
    `perl ../bin/Sage_TagMapping.pl $parDBfile.tags $parSourcefile.tab.tags 0`;
    print logfile "\nMapping Finished\n";
}
else{
    `split -l $parSegnumber $parSourcefile.tab.tags $parSourcefile.tab.tags.`;
    my $tagfile=`ls $parSourcefile.tab.tags.[a-z][a-z]`;
    my @tagfiles=split(/\n/,$tagfile);
    push(@fileCreated,@tagfiles);
    open (MapShell,">$parSourcefile.sh");
    push(@fileCreated,"$parSourcefile.sh");
    for(my $i=0;$i<@tagfiles;$i++){
	print logfile "perl $proTagMapping $parDBfile.tags $tagfiles[$i] $parMismatch\n";
	print MapShell "perl $proTagMapping $parDBfile.tags $tagfiles[$i] $parMismatch\n";
	push (@fileCreated,"$tagfiles[$i].map");
    }
    close MapShell;
    system ("perl $proJobMan $parSourcefile.sh");
    `cat $parSourcefile.tab.tags.[a-z][a-z].map > $parSourcefile.tab.tags.map`;
    print logfile "Mapping Finished\n";
}

## Match tagmap to the significant changed genes;

print logfile "\nMapping Tag Information to the Selected Tags\n";
for(my $i=0;$i<@fileCompair;$i++){
    print logfile "perl $proTagMapGene $parSourcefile.tab.tags.map $fileCompair[$i]\n";
    `perl $proTagMapGene $parSourcefile.tab.tags.map $fileCompair[$i]`;
}

if($parClearFile==1){
    print logfile "\nClean the Temporary Files\n";
    for (my $i=0;$i<@fileCreated;$i++){
	unlink $fileCreated[$i];
    }
}

## End of program
my $timecost=time-$timestart;
print logfile "\n\nAll Job Finished\n";

printf logfile ("\n\nTotal running time: %02d:%02d:%02d\n\n", int($timecost / 3600), int(($timecost % 3600) / 60), int($timecost % 60));
close logfile;
exit(0);


sub dealRTemp(){
    open (rtemp,$parRTemp);
    my @arrRTemp=<rtemp>;
    my $strRCmd=join("",@arrRTemp);
    my $colnames;

    $colnames="\"$sageTitle[0]\"";
    for(my $i=1;$i<@sageTitle;$i++){
	$colnames=$colnames.",\"$sageTitle[$i]\"";
    }
    
    ##replace the markers with real parameter;
    $strRCmd=~s/<Sourcefile>/$parSourcefile/g;
    $strRCmd=~s/<CountCutoff>/$parCountCutoff/g;
    $strRCmd=~s/<RDatafile>/$parRdatafile/g;
    $strRCmd=~s/<ColNames>/$colnames/g;
    my $strInsertCmd;
    if($strRCmd=~/<ExactTest>/){
	$strInsertCmd="d <- estimateCommonDisp\(d\)\nsource\(\"$proExactTest\"\)\n";
	my @arrayGroups=sort(keys %sageGroupHash);
	for(my $i=0;$i<@arrayGroups-1;$i++){
	    my $strInsertCmd2="exact$arrayGroups[$i]v$arrayGroups[$i+1] <- exactTesta\(d,";
	    my $strgroups="c\(\"$arrayGroups[$i]\",\"$arrayGroups[$i+1]\"\),c\(";
	   
	    my @groupA=split("#",$sageGroupHash{$arrayGroups[$i]});
	    my @groupB=split("#",$sageGroupHash{$arrayGroups[$i+1]});
	    print "@groupA\t@groupB\n";
	    $strgroups=$strgroups.($groupA[0]+1); ##group number from 0, colnumber from 1
	    for(my $j=1;$j<@groupA;$j++){
		$strgroups=$strgroups.",".($groupA[$j]+1);
	    }
	    for(my $j=0;$j<@groupB;$j++){
		$strgroups=$strgroups.",".($groupB[$j]+1);
    	    }
	    $strgroups="$strgroups\)";
	    $strInsertCmd="$strInsertCmd$strInsertCmd2$strgroups,$parTopnumber,\"$parSorttype\"\)";
	    my $filename="$parTargetDir\/exact$arrayGroups[$i]-$arrayGroups[$i+1].tab";
	    push(@fileCompair,$filename);
	    $strInsertCmd="$strInsertCmd\nwrite.table\(exact$arrayGroups[$i]v$arrayGroups[$i+1],file=\"$filename\",col.names=FALSE,row.names=FALSE,sep=\"\\t\"\)\n";
	}
	$strRCmd=~s/<ExactTest>/$strInsertCmd/;
    }
    
    return $strRCmd;
}


sub Filecheck{
    my ($Filename,$mode)=@_;
    if($mode==0){
	if(!(-e $Filename)){
	    print stderr "Can't find file: $Filename\n";
	    $iserror=1;
	}
    }
    else{
	if(-e $Filename){
	    print stderr "Directory $Filename already exists\n";
	    $iserror=1;
	}
    }
}



sub Usage(){
  print << "    Usage";
    
	Usage:  $0 (version $version)
		
	<options>
		-i     File contains all SAGE fastq filenames (must be given)
		-d     Fasta sequences for tag mapping (must be given)
		-a     Filename of R space file (default: sourcefile.Rdata)
		-c     Configure file, contains all the programs path (default: config.txt)
		-s     Scores for filtering the SAGE raw data (default: 20)
		-r     Model of output the redundent tags
		         0: remove all redundent tags
			 1: including gene unique redundent tags (default)
			 2: all tags
		-f     Cutoff the total number of the Tags (default: 5)
		-t     Numbers top significant result  (default:200)
		-g     Numbers of every segmental tag files, important for >0 missmatch mapping (default: 30000)
		-m     Maximum missmatch number of tag mapping (default: 0, maximum 3)
		-l     Log Filename (default: sourcefile.log)
		-st    Sort type of compare result (p.value 0|logFC 1, default:0)
		-rt    R template file (default: temple.R)
		-tl     Length of SAGE tags (default: 17)
		-cl     Clean the files created in the flow
		         0: Keep the files (default)
			 1: Clean the files
		-td    Target directory for save the result (default: sourcefile.dir)
		
		
    Usage

	exit(0);  
};


