# 3primeTag

3'mRNA-seq is a different approach to transcript analysis, sequencing only a short fragment at the 3' fragment extracted via polyA capturing technologies.
The main advantage is the reduced amount of sequence and therefore price and a clearer corelation between number of reads and number of transcript since each read represent one transcript.

An interesting pubpication to consult is "A comparison between whole transcript and 3’ RNA sequencing methods using Kapa and Lexogen library preparation methods' by Ma et al 2019 (https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-5393-3)




perl PolyACaptureNewV0.6.1.pl MappingV8/STAR_A4_01_2Aligned.out.bam A4-01 BED/STAR_A4_01_5-50-6-2.bed /biodata/ANDREAS/Cassava/ref/V8.1/annotation/Mesculenta_671_v8.1.gene_exons_plus.gff3 /biodata/ANDREAS/Cassava/ref/V8.1/assembly/SeqSplit LOG/STAR_A4_01_5 50-6-2.log 5 50 0.8 6 2



perl ParseBedMergeV0.3.pl BED/merged_A4_R1-D1.bed A4-01,A4-03,A4-05,A4-07,A4-09,A4-11 Results/merged_A4_R1-D1_pos.bed Results/merged_A4_R1-D1_id.txt Results/merged_A4_R1-D1.log 1 3

