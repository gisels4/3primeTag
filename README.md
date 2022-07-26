# 3primeTag
data analysis of 3printTag sequencing data


perl PolyACaptureNewV0.6.1.pl MappingV8/STAR_A4_01_2Aligned.out.bam A4-01 BED/STAR_A4_01_5-50-6-2.bed /biodata/ANDREAS/Cassava/ref/V8.1/annotation/Mesculenta_671_v8.1.gene_exons_plus.gff3 /biodata/ANDREAS/Cassava/ref/V8.1/assembly/SeqSplit LOG/STAR_A4_01_5 50-6-2.log 5 50 0.8 6 2



perl ParseBedMergeV0.3.pl BED/merged_A4_R1-D1.bed A4-01,A4-03,A4-05,A4-07,A4-09,A4-11 Results/merged_A4_R1-D1_pos.bed Results/merged_A4_R1-D1_id.txt Results/merged_A4_R1-D1.log 1 3

