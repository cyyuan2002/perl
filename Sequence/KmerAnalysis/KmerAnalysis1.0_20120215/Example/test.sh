perl ../Bin/format.pl 17mer.freq kmer.xls 1570530051
perl ../Bin/GenomeAnalysis.pl kmer.xls 99 17 71094946860
 perl ../Bin/test.pl 1000000 99 0.0062 17 54 0.82 0.018 Test
paste -d t kmer.xls Test/kmer.xls >test.xls
