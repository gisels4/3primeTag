# 3primeTag

3'mRNA-seq is a different approach to transcript analysis, sequencing only a short fragment at the 3' fragment extracted via polyA capturing technologies. The main advantage is the reduced amount of sequence and therefore price and a clearer corelation between number of reads and number of transcript since each read represent one transcript.

An interesting pubpication to consult is "A comparison between whole transcript and 3’ RNA sequencing methods using Kapa and Lexogen library preparation methods' by Ma et al 2019 (https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-5393-3).

In most cases, the 3'mRNA-seq data was analysed using the traditional approaches for the mRNA-seq data.

With the below scripts we developped an alternative approach to select clustering reads that represent a transcript by the presence of reads that contain the polyA tail.

We developed two scripts:

A) PolyACaptureNewV0.6.1.pl that uses a coordinate sorted bam file to extract the read clusters representing each transcript, indicating the presence of reads with polyA tails and counting the reads in each cluster. **It is IMPORTANT that befor using the aligner reads were NOT clean of the polyAtails!!!!**

B) ParseBedMergeV0.3.pl that will parse the information of the merged data of the outputs of PolyACaptureNewV0.6.1.pl, more details below, and create a count matrix for the statistical downstream analysis.

We exlain with an example the workflow

In the example we sequenced samples of two conditions and each has three replicates:

A-1, A-2, A-3 and B-1, B-2, B-3

and therefore we have 6 coordinates sorted bam files

A-1.bam, A-2.bam, A-3.bam and B-1.bam, B-2.bam, B-3.bam

foreach bam file we run the script PolyACaptureNewV0.6.1.pl to extract the clusters.
PolyACaptureNewV0.6.1.pl takes as imput:
1) the bam file
2) the output file name
3) the reference gff3 file
4) the path to the reference sequence, split into the single chromosome and/or scaffold
5) the log output file name
6) treshold for the minimum size of the cluster to be accepted
7) treahold for minimal size of matched fragments
8) ration total clip length and number of A's, the polyA which is clipped by the aligner contains not only As


perl PolyACaptureNewV0.6.1.pl MappingV8/STAR_A4_01_2Aligned.out.bam A4-01 BED/STAR_A4_01_5-50-6-2.bed /biodata/ANDREAS/Cassava/ref/V8.1/annotation/Mesculenta_671_v8.1.gene_exons_plus.gff3 /biodata/ANDREAS/Cassava/ref/V8.1/assembly/SeqSplit LOG/STAR_A4_01_5 50-6-2.log 5 50 0.8 6 2



perl ParseBedMergeV0.3.pl BED/merged_A4_R1-D1.bed A4-01,A4-03,A4-05,A4-07,A4-09,A4-11 Results/merged_A4_R1-D1_pos.bed Results/merged_A4_R1-D1_id.txt Results/merged_A4_R1-D1.log 1 3

