# 3primeTag

3'mRNA-seq is a different approach to transcript analysis, sequencing only a short fragment at the 3' fragment extracted via polyA capturing technologies. The main advantage is the reduced amount of sequence and therefore price and a clearer corelation between number of reads and number of transcript since each read represent one transcript.

An interesting pubpication to consult is "A comparison between whole transcript and 3’ RNA sequencing methods using Kapa and Lexogen library preparation methods' by Ma et al 2019 (https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-5393-3).

In most cases, the 3'mRNA-seq data was analysed using the traditional approaches for the mRNA-seq data.

With the below scripts we developped an alternative approach to select clustering reads that represent a transcript by the presence of reads that contain the polyA tail.

We developed two scripts:

A) PolyACaptureNewV0.6.1.pl that uses a coordinate sorted bam file to extract the read clusters representing each transcript, indicating the presence of reads with polyA tails and counting the reads in each cluster. **It is IMPORTANT that befor using the aligner, reads were NOT clean of the polyAtails!!!!**

B) ParseBedMergeV0.3.pl that will parse the information of the merged data of the outputs of PolyACaptureNewV0.6.1.pl, more details below, and create a count matrix for the statistical downstream analysis.

We exlain with an example the workflow

In the example we sequenced samples of two conditions and each has three replicates:

A-1, A-2, A-3 and B-1, B-2, B-3

and therefore we have 6 coordinates sorted bam files (we use STAR and picard MergeSamFiles)

STAR-A-1.bam, STAR-A-2.bam, STAR-A-3.bam and STAR-B-1.bam, STAR-B-2.bam, STAR-B-3.bam

foreach bam file we run the script PolyACaptureNewV0.6.1.pl to extract the clusters.

PolyACaptureNewV0.6.1.pl takes as imput:
1) the bam file
2) sample name, **must not contain "_"!**
3) the output file name
4) the reference gff3 file
5) the path to the reference sequence, split into the single chromosome and/or scaffold
6) the log output file name
7) treshold for the minimum size of the cluster to be accepted (5)
8) treahold for minimal size of matched fragments (50)
9) treahold for the ration total clip length and number of A's, the polyA which is clipped by the aligner contains not only A's (0.8)
10) treshold for the minimal number of consecutive A in polyA tail (5)
11) treshold for min number of polyA sequences in the cluster to be accepted (2)
12) treshold for the number of reads to confirm the splicing site (intron) (2)

example comand:
perl PolyACaptureNewV0.6.1.pl MappingV8/Star-A-1.bam A-1 STAR_A-1_5-50-0.8-5-2-2.txt /biodata/Cassava/ref/V8.1/annotation/Mesculenta_671_v8.1.gene_exons.gff3 /biodata/Cassava/ref/V8.1/assembly/SeqSplit STAR_A-1_5-50-0.8-5-2-2.log 5 50 0.8 5 2 2

the main output file is STAR_A-1_5-50-0.8-5-2-2.txt and it is formated as a bed file format. However it contains at the end of each line (cluster) the position(s) of the polyA atttachement site.

To continue for a the expression analysis we have to remove this last data wit the following comand:

cut -f 1,2,3,4,5,6 STAR_A-1_5-50-0.8-5-2-2.txt > STAR_A-1_5-50-0.8-5-2-2.bed

for the differential expression analysis we need to combine the 2 samples with the 3 replicates and sort the data by the coordinates:

cat STAR_A-1_5-50-0.8-5-2-2.bed STAR_A-2_5-50-0.8-5-2-2.bed STAR_A-3_5-50-0.8-5-2-2.bed STAR_B-1_5-50-0.8-5-2-2.bed STAR_B-2_5-50-0.8-5-2-2.bed STAR_B-3_5-50-0.8-5-2-2.bed   | sort -k1,1 -k2,2n > Star-A-B.bed

As a next step we merge the clusters according to the position data using bedtools:

bedtools merge -s -c 4,5 -o collapse -i Star-A-B.bed > merged_Star-A-B.bed

As last step we parse the count data and create the count matrix for the statistical downstream analysis using ParseBedMergeV0.3.pl

PolyACaptureNewV0.6.1.pl takes as imput:
1) the merged bed file
2) a coma separated list of sample names
3) output file name for count matrix with coordinates
4) output file name for count matrix with gene names
5) log file name
6) treshold for minimum polyA containing clusters between samples
7) treshold for minimum clusters between the samples
8) print all isoforms (all) or only one (nr)

example comand:
perl ParseBedMergeV0.3.pl merged_Star-A-B.bed A-1,A-2,A-3,B-1,B-2,B-3 merged_Star-A-B_pos.bed merged_Star-A-B_id.txt merged_Star-A-B.log 1 3 iso

The file you are interesteds is merged_Star-A-B_id.txt looking like:

TranscriptID	ClusterStart	 A-1  A-2  A-3  B-1  B-2  B-3
Manes.01G000031.1.v8.1	19594	1398	1947	1531	1917	1579	1579
Manes.01G000031.1.v8.1	22474	17	22	7	15	13	13

for the statistical downstream analysis you have to parse out the ClusterStart column; we added that since we observed in some cases different cluster for the same gene. We select the clsuter with the highest count numbers.

cut -f 1,3,4,5,6,7,8 merged_Star-A-B_id.txt > matrix_Star-A-B_id.txt

**The workflow is quite clumsy and we are working on a more elegant version!!**

