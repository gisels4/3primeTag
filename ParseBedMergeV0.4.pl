#!/usr/local/bin/perl -w
############################################################################################################
#	ParseBedMergeV0.4.pl
#
#	the three replicates from PolyACaputreNewV0.6.pl get merged with cat and sorted by coordinates (sort -k 1,1 -k2,2n)
#	this file then gets merged with bedtools merge -s -c 4,5 -o collapse -i and the following output gets parsed here
#
#	Created by Andreas Gisel July 2022
#
##########################################################################################################


use strict;

my $file_in = shift;						# file from bedtools merge -s -c 4,5 -o collapse -i
my $samples = shift;						# sample names separated with coma A4-01,A4-03,A4-05,B4-01,B4-03,B4-05
my $file_out_count = shift;					# counts of all clusters fulfilling the criteria with position information
my $file_out_id = shift;					# counts of all clusters fulfilling the criteria with gene id information
my $file_log = shift;
my $polyA_count = shift;					# number pf polyA between samples required
my $clu_count = shift;						# number of clusters between samples required
my $iso = shift;							# print all isoforms (all) or only one (nr)

$| = 1;

my @samples = split/,/,$samples;

unless(open(LOG, ">$file_log"))
{
	print "cannot open $file_log!!\n";
	exit;
}


my @timeData = localtime(time);
print LOG "start reading file\t@timeData\n";
print "start reading file\t@timeData\n";

unless(open(IN, $file_in))
{
	print "cannot open $file_in!!\n";
	exit;
}

my %cluster_pos;								# hosts cluster information (according position) of all samples; chrom, cluster start, sampleID and count
my %cluster_id;									# hosts cluster information (gene id) of all samples; chrom, cluster start, sampleID and count
my %cluster_polyA;								# hosts cluster polyA tail info of all samples; chrom, cluster start, polyA (A/N)


while(<IN>)
{
	chomp;
	my @data = split/\t/;
	
		# Chromosome01	19594	20604	-	A4-15_R_7571_A_Manes.01G000031.1.v8.1&&Manes.01G000031.2.v8.1&&Manes.01G000031.3.v8.1_0,A4-13_R_7172_A_Manes.01G000031.1.v8.1&&Manes.01G000031.2.v8.1&&Manes.01G000031.3.v8.1_0,A4-17_R_6204_A_Manes.01G000031.1.v8.1&&Manes.01G000031.2.v8.1&&Manes.01G000031.3.v8.1_0,A4-01_R_8532_A_Manes.01G000031.1.v8.1&&Manes.01G000031.2.v8.1&&Manes.01G000031.3.v8.1_0,A4-03_R_10602_A_Manes.01G000031.1.v8.1&&Manes.01G000031.2.v8.1&&Manes.01G000031.3.v8.1_0,A4-05_R_8511_A_Manes.01G000031.1.v8.1&&Manes.01G000031.2.v8.1&&Manes.01G000031.3.v8.1_0,A4-15_R_7572_N_Manes.01G000031.1.v8.1&&Manes.01G000031.2.v8.1&&Manes.01G000031.3.v8.1_0,A4-17_R_6205_N_Manes.01G000031.1.v8.1&&Manes.01G000031.2.v8.1&&Manes.01G000031.3.v8.1_0	1675,1398,1523,1316,1405,1723,272,8
		# Chromosome01	22474	22835	-	A4-13_R_7175_A_Manes.01G000031.1.v8.1&&Manes.01G000031.2.v8.1&&Manes.01G000031.3.v8.1_0,A4-03_R_10604_N_Manes.01G000031.1.v8.1&&Manes.01G000031.2.v8.1&&Manes.01G000031.3.v8.1_0,A4-15_R_7574_N_Manes.01G000031.1.v8.1&&Manes.01G000031.2.v8.1&&Manes.01G000031.3.v8.1_0,A4-05_R_8512_N_Manes.01G000031.1.v8.1&&Manes.01G000031.2.v8.1&&Manes.01G000031.3.v8.1_0	17,13,22,6
	
	my @id = split/,/, $data[4];			# splits the cluster information of the different samples
	my @count = split/,/, $data[5];			# splits the counts corresponding to the clusters above (@id)
	my @info;								# hosts for each cluster in @id: sampleID, orientation(F/R), clusterID, geneID(s)
	my %polyA;								# counts polyAs for the cluster in the samples
	my %clus; 								# counts the cluster from the different samples
	
	foreach my $check (@id)
	{
		@info = split/_/,$check;			# separate the info within the sample ID
		
		$clus{$info[0]}++;					# counts samples

		if($info[3] eq "A")
		{
			$polyA{$info[0]}++;				# counts samples with poly A
		}
	}
	
	# if(keys %polyA >= $polyA_count and keys %clus >= $clu_count)
	if(keys %polyA >= $polyA_count)
	{
		for(my$i=0;$i<=$#id;$i++)
		{
			@info = split/_/,$id[$i];
			if($cluster_pos{$data[0]}{$data[1]}{$info[0]})
			{
				$cluster_pos{$data[0]}{$data[1]}{$info[0]} += $count[$i];
			}
			else
			{
				$cluster_pos{$data[0]}{$data[1]}{$info[0]} = $count[$i];
			}
			
			my @ids = split/&&/, $info[4];
			foreach my $id (@ids)
			{
				if($cluster_id{$id}{$data[1]}{$info[0]})						# key defined by gene ID, position and sample
				{
					$cluster_id{$id}{$data[1]}{$info[0]} += $count[$i];
				}
				else
				{
					$cluster_id{$id}{$data[1]}{$info[0]} = $count[$i];
				}
			}
		}
	}
}

unless(open(OUT, ">$file_out_count"))											# print cluster positions
{
	print "cannot open $file_out_count!!\n";
	exit;
}



for my $chrom ( sort keys %cluster_pos )
{
	for my $start ( sort { $a <=> $b } keys %{ $cluster_pos{$chrom}})
	{
		print OUT "$chrom\t$start";
		
		foreach my $sample (@samples)
		{
			if($cluster_pos{$chrom}{$start}{$sample})
			{
				print OUT "\t$cluster_pos{$chrom}{$start}{$sample}";
			}
			else
			{
				print OUT "\t0";
			}
		}
		print OUT "\n";
	}
}

close OUT;

unless(open(OUT, ">$file_out_id"))											# print cluster positions
{
	print "cannot open $file_out_id!!\n";
	exit;
}

print OUT "TranscriptID\tClusterStart";
foreach my $sample (@samples)
{
	print OUT "\t$sample"
}
print OUT "\n";

my %nrprint;

for my $id ( sort keys %cluster_id )
{
	for my $start ( sort keys %{ $cluster_id{$id}})
	{
		# print OUT "$id-$start";
		my $print = "A";
		my $pres = 0;
		
		foreach my $sample (@samples)
		{
			if($print eq "A")
			{
				if($cluster_id{$id}{$start}{$sample})
				{
					$print = "$cluster_id{$id}{$start}{$sample}";
					$pres++;

				}
				else
				{
					$print = "0";
				}
			}
			else
			{
				if($cluster_id{$id}{$start}{$sample})
				{
					$print = "$print\t$cluster_id{$id}{$start}{$sample}";
					$pres++;
				}
				else
				{
					$print = "$print\t0";
				}
			}
		}
		
		if($iso eq "all")
		{
			if($pres >= $clu_count)
			{
				print OUT "$id-$start\t$print\n";
			}
		}
		else
		{
			if(not $nrprint{$print} and $pres >= $clu_count)
			{
				$nrprint{$print} = "$id-$start";
			}
		}
	}
}

if($iso ne "all")
{
	for my $counts( sort { $nrprint{$a} cmp $nrprint{$b} } keys %nrprint)
	{
		print OUT "$nrprint{$counts}\t$counts\n";
	}
}

exit;
