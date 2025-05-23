# 3primeTag

3'mRNA-seq is a different approach to transcript analysis, sequencing only a short fragment at the 3' fragment extracted via polyA capturing technologies. The main advantage is the reduced amount of sequence and therefore price and a clearer corelation between number of reads and number of transcript since each read represent one transcript.

An interesting pubpication to consult is "A comparison between whole transcript and 3’ RNA sequencing methods using Kapa and Lexogen library preparation methods' by Ma et al 2019 (https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-5393-3).

In most cases, the 3'mRNA-seq data was analysed using the traditional approaches for the mRNA-seq data.

With the below scripts we developped an alternative approach to select clustering reads that represent a transcript by the presence of reads that contain the polyA tail.

We developed two scripts:
A) PolyACaptureNewV0.6.1.pl that used a coordinate sorted bam file to extract the read clusters representing each a transcrip, indicating the presence of reads with polyA tails and counting the reads in each cluster.

B) ParseBedMergeV0.3.pl that will parse the information of the merged data of the outputs of PolyACaptureNewV0.6.1.pl, more details below, and create a count matrix for the statistical downstream analysis.


perl PolyACaptureNewV0.6.1.pl MappingV8/STAR_A4_01_2Aligned.out.bam A4-01 BED/STAR_A4_01_5-50-6-2.bed /biodata/ANDREAS/Cassava/ref/V8.1/annotation/Mesculenta_671_v8.1.gene_exons_plus.gff3 /biodata/ANDREAS/Cassava/ref/V8.1/assembly/SeqSplit LOG/STAR_A4_01_5 50-6-2.log 5 50 0.8 6 2



perl ParseBedMergeV0.3.pl BED/merged_A4_R1-D1.bed A4-01,A4-03,A4-05,A4-07,A4-09,A4-11 Results/merged_A4_R1-D1_pos.bed Results/merged_A4_R1-D1_id.txt Results/merged_A4_R1-D1.log 1 3

