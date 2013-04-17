Author: Dongliang Zhan, zhandongliang@genomics.org.cn
Version: 1.0, Date: 2012-2-15

This is kmer analysis for short reads to estimate genome size,repeat rate,heterozygosity.

Before you use this program,you have to install jellyfish to do simulation,and change the config file in the Bin Directory.

ok,let's begin to do the kmer analysis!

(1) use the jellyfish to count kmer frequency for solexa reads,by experience 30-60X data can get better results.
    you shoud reverse and completement the kmer.
    ps:
	you can use other software to count kmer frequency,but the software must contains information below:
    <1> total kmer number
    <2> total kmer type
    <3> a table contains information like:
	1	17141
	2	719
	3	161
	4	210
	5	250
	6	218
	7	205
	8	312
	9	357
	
	col.1 is the frequency of kmer, col.2 is the number of kmer type

(2) format the table file(format like <3>) using Bin/format.pl
    usage:
    perl Bin/format.pl <kmer.histo> <outFile> <number of kmer_kind(or node num):from *.stats District>
    parameters:
    <kmer.histo>:	the table in the format like <3>
    <outFile>:		the output file
    <kmer_kind>:	the number of kmer kind,you can get this value by using command like "grep Distinct *.stats" if you use jellyfish.

(3) do analysis using Bin/GenomeAnalysis.pl
    usage:
    perl ../Bin/GenomeAnalysis.pl <kmer.xls> <read_len> <k> <kmer_num> <expect_peak:default is the most high peak>
    parameters:
    <kmer.xls>: the format file generate by (2),the format is like:
    		1	0.393228927	617577847
		2	0.086437281	135752348
		3	0.038148596	59913517
		4	0.026984766	42380386
		5	0.022455646	35267267
		6	0.019265352	30256814
		7	0.016292055	25587162
		8	0.013584721	21335212
     col.1 is the frequency of kmer, col.2 is the percentage of kmer type, col.3 is the number of kmer type.

    <read_len>: the length of read
    <k>:	the length of kmer.
    <kmer_num>: the number of kmer from reads,you can get this value by using command like "grep Total *.stats" if you use jellyfish.
    <expect_peak>:this program may not find the peak if the peak is not outstanding,or the peak is no the main peak caused by error or snp,so you can set a peak to get more accurate result.

output:
***********************************************
peak:20
genomeSize:3501306188
repeat:0.920835298771084
hete:0.0509038265090368
read coverage:24.2195318843416

***********************************************
peak:40
genomeSize:1706215151
repeat:0.817742016882012
hete:0.0180414532931246
read coverage:49.7006469596802

This program uses each peak to do kmer analysis.

(4) simulation
    This module is used to simulate a genome with snp,repeat,and generate reads with errors without indels.And finally generate a kmer.xls.
    Note that the genome is very complex,the error module is also complex too.This simulation part can not match the real kmer distribution of the genome.

    usage:
    perl Bin/test.pl <size> <read_len> <error_rate> <kmer_size> <coverage> <repeat_rate> <hete_rate> <outDir>
    parameters:
    <size>:		the size of genome you want to simulate, by experience, 1000000  does well.
    <read_len>:		the read length
    <error_rate>:	from 0.001 to 0.01,by experience,the error rate is ~0.006 after filtering low quality data.
    <kmer_size>: 	the length of kmer
    <coverage>:		the bases of reads minus the genome size.for example,20 means we generate 20*G bp reads.
    <repeat_rate>:	the repeat rate of the genome.
    <hete_rate>:	we only simulate snp
    <outDir>:		the directory to put the simulate data into 

     The result derived by (3) may contains ~5% deviation,you may change the parameters to get a more likely kmer distribution picture.

(5) example
    u can look at the Example/test.sh to see the demo.
    

