#!/usr/local/bin/perl -w
############################################################################################################
#	PolyACaputreNewV0.6.pl
#	
#
#
#	Created by Andreas Gisel June 2022
#
##########################################################################################################


use strict;

use Bio::SeqIO;

my $file_in = shift;						# coordinate sorted BAM file
my $sample = shift;							# sample nname
my $file_out_all = shift;					# bed file of all clusters fulfilling the criteria
my $gff_file = shift;
my $path = shift;							# path to reference sequence (directory with ref seq for each chromosome)
my $file_log = shift;
my $cluster_size = shift;					# min cluster size								5
my $min_match = shift;						# minimal size of matched fragments 			55
my $ratio = shift;							# ration total clip length and number of A's	0.8
my $polyA_min = shift;						# min number of consecutive A in polyA tail		7
my $num_polyA = shift;						# min number of polyA sequences in the cluster	1
my $intron_tresh = shift;					# treshold for the number of reads to confirm the intron 2

$| = 1;

unless(open(LOG, ">$file_log"))
{
	print "cannot open $file_log!!\n";
	exit;
}


my @timeData = localtime(time);
print LOG "start reading BAM file\t@timeData\n";
print "start reading BAM file\t@timeData\n";

unless(open(IN, "samtools view $file_in |"))
{
	print "cannot open $file_in!!\n";
	exit;
}

my $Clu_start = 0;
my $Clu_end;
my $Clu_id;
my $Clu_counter = 1;
my %Clu_count;
my %Clu_pos;
my $chr = "A";
my %polyA;						# tells me whether the cluster contains a sequence with a polyA
my $chrom;						# whole chromosome sequence
my $counter = 0;
my $intron_ev = 0;					# counts the reads which gives the evedence for the intron

my $end_old = 0;					

while(<IN>)
{
	# NS500352:278:HHYCVBGX7:1:11208:17290:15192	0	Chromosome10	1314089	255	1S61M14S	*	0	0	GTCTGGCTGGGTGTATTCGTATGATTGTAATTATGTGTATGGTTGTGATTCCATTTCTGCGAAAAAAAAAAAAAAA	AAAAAEAEEAEEEE6EEAE/A<EEAEEAEAEEEEEEEAAEEEEE/EA/EAEAEEEE<E/E/E/EEEEA/<EEEEEE	NH:i:1	HI:i:1	AS:i:58	nM:i:1
	
	next if(/@/);										# remove SAM header
	next if (/NH:i:[2-9]/);								# remove multiple hits
	
	# $counter++;
	# if(($counter % 100000) == 0 )
	# {
	# 	print "$counter\n";
	# }
	
	my @data = split/\t/;
	
	if($chr ne $data[2])
	{
		print "Chr: $chr - $data[2]\n";
		$chr = $data[2];
		my $file;
				
		if($chr =~ /ome(\d\d)/)
		{
			$file = "$path/chromosome$1.fasta";
		}
		elsif($chr =~ /fold(\d\d)/)
		{
			$file = "$path/scaffold$1.fasta";
		}
		else
		{
			print "problem with ref files $file!\n";
			exit;
		}

		print "$file\n";
		
		my $inseq = Bio::SeqIO->new(-file   => "$file");
		my $seq = $inseq->next_seq;
		$chrom = $seq->seq();
	}
	
	my @end = findEnd($data[3], $data[5]);				# input hit position, CIGAR. Calculates according to the CIGAR the end position of the read
									# return $end1, $end2, $intron[0], $flag, @clip
	my $end;
	
	if($end[2] > 0)
	{
		if($end[2] == $end_old)
		{
			$intron_ev++;
		}
		else
		{
			$end_old = $end[2];
			$intron_ev = 1;		
		}
	}
	
	if($intron_ev > $intron_tresh)				# if $intron_ev larger than the treshold we will accept the intron und use the corresponding end position
	{
		$end = $end[0];
	}
	else
	{
		$end = $end[1];
	}
	


	# print "$end[2] - $intron_ev - $end - $end[0] - $end[1]\n";
	# if($end[2]>0)
	# {
	# 	sleep(1);
	# }

	
	next if($end[3] == 4);								# hit fragment shorter than min_match
	
	if($data[3] >= $Clu_start)						# if new hit position is bigger than last hit we are still on the same chromosome
	{
		if($Clu_start == 0)
		{
			$Clu_start = $data[3];
			$Clu_end = $data[3] + 100;
		}
		
		if($data[3] <= $Clu_end)					# if the new hit position is smaller than the last end position the read belongs to the cluster
		{
			if($end[0] > $Clu_end)
			{
				$Clu_end = $end;
			}
			
			if($data[1] == 0)
			{
				$Clu_id = $sample ."_F_$Clu_counter";
				
				if($end[3] == 2 or $end[3] == 3)
				{
					polyA($end[5], $data[9], $Clu_id, $end, $data[2]);		# check for polyA (polyT) - send clip data from above, sequence
				}
				else
				{
					print LOG "$Clu_id\t$data[9]\t$data[3]\t$data[5]\n";
				}

			}
			elsif($data[1] == 16)
			{
				$Clu_id = $sample ."_R_$Clu_counter";
				
				if($end[3] == 1 or $end[3] == 3)
				{
					polyT($end[4], $data[9], $Clu_id, $data[3], $data[2]);		# check for polyA (polyT) - send clip data from above, sequence
				}
				else
				{
					print LOG "$Clu_id\t$data[9]\t$data[3]\t$data[5]\n";				
				}
			}
			
			$Clu_count{$Clu_id}++;
		}
		else										# new cluster
		{
			if($data[1] == 0)
			{
				$Clu_pos{$data[2]}{$Clu_id}[0] = $Clu_start;
				$Clu_pos{$data[2]}{$Clu_id}[1] = $Clu_end;
			}
			elsif($data[1] == 16)
			{
				$Clu_pos{$data[2]}{$Clu_id}[0] = $Clu_start;
				$Clu_pos{$data[2]}{$Clu_id}[1] = $Clu_end;
			}
			
			$Clu_counter++;
			$Clu_end = $end;
			$Clu_start = $data[3];
			
			if($data[1] == 0)
			{
				$Clu_id = $sample ."_F_$Clu_counter";
				
				if($end[3] == 2 or $end[3] == 3)
				{
					polyA($end[5], $data[9], $Clu_id, $end, $data[2]);		# check for polyA (polyT) - send clip data from above, sequence
				}
				else
				{
					print LOG "$Clu_id\t$data[9]\t$data[3]\t$data[5]\n";				
				}

			}
			elsif($data[1] == 16)
			{
				$Clu_id = $sample ."_R_$Clu_counter";
				
				if($end[3] == 1 or $end[3] == 3)
				{
					polyT($end[4], $data[9], $Clu_id, $data[3], $data[2]);		# check for polyA (polyT) - send clip data from above, sequence
				}
				else
				{
					print LOG "$Clu_id\t$data[9]\t$data[3]\t$data[5]\n";				
				}
			}

			$Clu_count{$Clu_id}++;
		}
	}
	else											# new chromosome
	{
		if($data[1] == 0)
		{
			$Clu_pos{$chr}{$Clu_id}[0] = $Clu_start;
			$Clu_pos{$chr}{$Clu_id}[1] = $Clu_end;
		}
		elsif($data[1] == 16)
		{
			$Clu_pos{$chr}{$Clu_id}[0] = $Clu_start;
			$Clu_pos{$chr}{$Clu_id}[1] = $Clu_end;
		}
		
		$Clu_counter++;
		$Clu_end = $end;
		$Clu_start = $data[3];
		
		if($data[1] == 0)
		{
			$Clu_id = $sample ."_F_$Clu_counter";
			
			if($end[3] == 2 or $end[3] == 3)
			{
				polyA($end[5], $data[9], $Clu_id, $end);			# check for polyA (polyT) - send clip data from above, sequence
			}
			else
			{
					print LOG "$Clu_id\t$data[9]\t$data[3]\t$data[5]\n";			
			}
		}
		elsif($data[1] == 16)
		{
			$Clu_id = $sample ."_R_$Clu_counter";
			if($end[3] == 1 or $end[3] == 3)
			{
				polyT($end[4], $data[9], $Clu_id, $data[3]);		# check for polyA (polyT) - send clip data from above, sequence
			}
			else
			{
					print LOG "$Clu_id\t$data[9]\t$data[3]\t$data[5]\n";			
			}
		}

		$Clu_count{$Clu_id}++;
	}
}

@timeData = localtime(time);
print LOG "done reading BAM file, starting GFF reading\t@timeData\n";
print "done reading BAM file, starting GFF reading\t@timeData\n";

unless(open(IN, $gff_file))
{
	print "cannot open $gff_file!!\n";
	exit;
}


my %gff_id;							# positions of mRNA and ids						

while(<IN>)
{
	next if(/^##/);
	
	# Chromosome01	phytozomev13	mRNA	19589	28327	.	-	.	ID=Manes.01G000031.1.v8.1;Name=Manes.01G000031.1;pacid=50282548;longest=1;ancestorIdentifier=Manes.12G099301.1.v7.1;Parent=Manes.01G000031.v8.1
	# Chromosome01	phytozomev13	prom_200	28328	28527	.	-	.	ID=Manes.01G000031.1.v8.1.prom200;Parent=Manes.01G000031.1.v8.1

	
	my @data = split/\t/;
	my $id;
	
	if($data[2] =~ /mRNA/)	 	# get feature ID
	{
		my @dataID = split/;/,$data[8];
		if($dataID[0] =~ /ID=(.*)/)
		{
			$id = $1;
		}
		else
		{
			print "problem with ID extraction!\n";
		}
		
		# for my $i($data[3]..$data[4])
		for(my$i=$data[3];$i<=$data[4];$i=$i+1)							# register only each third position of this feature to speed up ad safe space
		{
			if($gff_id{$data[0]}{$data[6]}{$i})
			{
				$gff_id{$data[0]}{$data[6]}{$i} = "$gff_id{$data[0]}{$data[6]}{$i}&&$id";
			}
			#elsif($gff_id{$data[0]}{$data[6]}{$i+1})
			#{
			#	$gff_id{$data[0]}{$data[6]}{$i} = "$gff_id{$data[0]}{$data[6]}{$i+1}&&$id";
			#}
			#elsif($gff_id{$data[0]}{$data[6]}{$i+2})
			#{
			#	$gff_id{$data[0]}{$data[6]}{$i} = "$gff_id{$data[0]}{$data[6]}{$i+2}&&$id";
			#}
			else
			{
				$gff_id{$data[0]}{$data[6]}{$i} = "$id";
			}
		}
	}
}

my $size = keys %gff_id;

@timeData = localtime(time);
print LOG "done with GFF reading $size - print to file \t@timeData\n";
print "done with GFF reading $size  - print to file \t@timeData\n";


unless(open(OUT, ">$file_out_all"))											# print cluster positions
{
	print "cannot open $file_out_all!!\n";
	exit;
}


for my $chrom ( keys %Clu_pos )
{
	for my $clus ( keys %{ $Clu_pos{$chrom}})
	{
		if($Clu_count{$clus} >= $cluster_size)
		{
			my $size = $Clu_pos{$chrom}{$clus}[1] - $Clu_pos{$chrom}{$clus}[0];
			
			if($clus =~ /F/)
			{
				my $strand = "+";
				my $hsize = 0;
				
				for my $clus1 ( keys %{ $polyA{$clus})
				{
					$hsize += $polyA{$clus}{$clus1};
				}
				# my $hsize = keys %{ $polyA{$clus}};
				
				print LOG "hsizeF $clus : $hsize\n";
				my @id = findID($chrom,$strand,$Clu_pos{$chrom}{$clus}[0],$Clu_pos{$chrom}{$clus}[1]);


				if($hsize >= $num_polyA)
				{
					my $clusA = $clus."_A";
					print OUT "$chrom\t$Clu_pos{$chrom}{$clus}[0]\t$Clu_pos{$chrom}{$clus}[1]\t". $clusA ."_$id[0]_$id[1]\t$Clu_count{$clus}\t$strand\t$size";
				}
				else
				{
					my $clusA = $clus."_N";
					print OUT "$chrom\t$Clu_pos{$chrom}{$clus}[0]\t$Clu_pos{$chrom}{$clus}[1]\t". $clusA ."_$id[0]_$id[1]\t$Clu_count{$clus}\t$strand\t$size";
				}
								
				for my $pos (keys %{ $polyA{$clus}})
				{
					print OUT "\t$pos-$polyA{$clus}{$pos}";
				}
				
				print OUT "\n";
			}
			elsif($clus =~ /R/)
			{
				my $strand = "-";
				my $hsize = 0;
				
				for my $clus1 ( keys %{ $polyA{$clus})
				{
					$hsize += $polyA{$clus}{$clus1};
				}
				# my $hsize = keys %{ $polyA{$clus}};
				
				print LOG "hsizeR $clus : $hsize\n";
				
				my @id = findID($chrom,$strand,$Clu_pos{$chrom}{$clus}[0],$Clu_pos{$chrom}{$clus}[1]);

				
				if($hsize >= $num_polyA)
				{
					my $clusA = $clus."_A";
					print OUT "$chrom\t$Clu_pos{$chrom}{$clus}[0]\t$Clu_pos{$chrom}{$clus}[1]\t". $clusA ."_$id[0]_$id[1]\t$Clu_count{$clus}\t$strand\t$size";
				}
				else
				{
					my $clusA = $clus."_N";
					print OUT "$chrom\t$Clu_pos{$chrom}{$clus}[0]\t$Clu_pos{$chrom}{$clus}[1]\t". $clusA ."_$id[0]_$id[1]\t$Clu_count{$clus}\t$strand\t$size";
				}
				
#				my @id = findID($chrom,$strand,$Clu_pos{$chrom}{$clus}[0],$Clu_pos{$chrom}{$clus}[0]);
#
#				print OUT "\t$id[0]-$id[1]";


				for my $pos (keys %{ $polyA{$clus}})
				{
					print OUT "\t$pos-$polyA{$clus}{$pos}";
				}
				
				print OUT "\n";
			}
			else
			{
				print "sh...... did not work!\n";
			}
		}
	}
}

@timeData = localtime(time);
print LOG "done\t@timeData  - cluster number $Clu_counter\n";
print "done\t@timeData - cluster number $Clu_counter\n";


sub findEnd
{
	my @list = @_;				#input: $data[3] (position), $data[5] (CIGAR)
	
	my @clip =(0,0);			# array with clip sizes
	my @match = (0,0,0);			# array with match sizes
	my @intron = (0,0,0);			# array with intron size
	my $flag = 0;				# flag for 5'S (1), 3'S (2), both (3), M too short (4) or intron (5)
	
	
	if($list[1] =~ /^(\d+)S/)
	{
		$clip[0] = $1;
		$flag = 1;
	}
	
	if($list[1] =~ /(\d+)S$/)
	{
		$clip[1] = $1;
		
		if($flag == 1)
		{
			$flag = 3;
		}
		else
		{
			$flag = 2;
		}
	}

	my $i = 0;
	while($list[1] =~ /(\d+)M/g)
	{
		$match[$i] = $1;
		$i++;
	}
	
	
	$i = 0;
	while($list[1] =~ /(\d+)N/g)
	{
		$intron[$i] = $1;
		$i++;
	}
	
	my $match_size1 = 0;			# match size with evtl intron
	my $match_size2 = 0;			# match size without intron
	
	foreach my $n (@match)
	{
		$match_size1 += $n;
		$match_size2 += $n;
	}
	
	if($match_size2 < $min_match)
	{
		$flag = 4;
	}

	foreach my $n (@intron)
	{
		$match_size1 += $n;
	}
	
	my $end1 = $list[0] + $match_size1;
	my $end2 = $list[0] + $match_size2;
	
	return ($end1, $end2, $intron[0], $flag, @clip);
}


sub polyA 
{
	my @list = @_;									#input  $end[1] (right clip length), $data[9] (sequence), $Clu_id, $end[0] ( end of read), $data[2] (Chromosome)
	
	my $tail = substr $list[1] ,-$list[0];				# get the sequence of the soft clipped fragment
	my $len = length($tail);
		
	my $query = check_query($list[4], $list[3], $list[0], "A");							#send chrom, hit start, skip length, orientation
	
	my @polyA = split //,$tail;
	my $count = 0;
	for(my$j=1;$j<=$#polyA;$j++)
	#foreach my $j (@polyA)
	{
		if($polyA[$j] eq "A")
		{
			$count++;
		}
		else
		{
			last;
		}
	}
	
	#my $len = length($tail);
	
	#if($count > 0 and $len/$count > $ratio)
	if($count >= $polyA_min)
	{
		if($polyA{$list[2]})						# check whether a polyA(T) is already assigned to the cluster
		{
			my $pos = "A-$list[3]";
			$polyA{$list[2]}{$pos}++;
	
			# $polyA{$list[2]} = "$polyA{$list[2]}&&A-$list[3]";
			print LOG "$list[2]\t$list[1]\t$tail\t$query\t$len\t$list[3]\t$polyA{$list[2]}{$pos}\t$count\n";
		}
		else
		{
			my $pos = "A-$list[3]";
			$polyA{$list[2]}{$pos}++;

			#$polyA{$list[2]} = "A-$list[3]";
			print LOG "$list[2]\t$list[1]\t$tail\t$query\t$list[3]\t$polyA{$list[2]}{$pos}\t$count\n";
		}
	}
}

sub polyT
{
	my @list = @_;									#input  $end[1] (right clip length), $data[9] (sequence), $Clu_id, $data3 ( start of the read)

	my $tail = substr $list[1], 0, $list[0];	    # get the sequence of the soft clipped fragment
	my $len = length($tail);
	my $query = check_query($list[4], $list[3], $list[0], "T");							#send chrom, hit start, skip length, orientation


	my @polyT = split //,$tail;
	my $count = 0;
	for(my$j=$#polyT-1;$j>=0;$j--)
	{
		if($polyT[$j] eq "T")
		{
			$count++;
		}
		else
		{
			last;
		}
	}
	
	# my $len = length($tail);
	
	# if($count > 0 and $len/$count > $ratio)
	if($count >= $polyA_min)
	{
		if($polyA{$list[2]})						# check whether a polyA(T) is already assigned to the cluster
		{
			my $pos = "T-$list[3]";
			$polyA{$list[2]}{$pos}++;

			# $polyA{$list[2]} = "$polyA{$list[2]}&&T-$list[3]";
			print LOG "$list[2]\t$list[1]\t$tail\t$query\t$len\t$list[3]\t$polyA{$list[2]}{$pos}\t$count\n";
		}
		else
		{
			my $pos = "T-$list[3]";
			$polyA{$list[2]}{$pos}++;

			# $polyA{$list[2]} = "T-$list[3]";
			print LOG "$list[2]\t$list[1]\t$tail\t$query\t$list[3]\t$polyA{$list[2]}{$pos}\t$count\n";
		}

	}
}

sub findID
{
	my @list = @_;							# input: $chrom (chromosome), $strand, $Clu_pos{$chrom}{$clus}[0] (cluster start), $Clu_pos{$chrom}{$clus}[1] (cluster end)
	
	my $rangeS1 = $list[2]-1;							 # shift to be find the feature since we have only every third position
	my $rangeE1 = $list[3]-1;
	my $rangeS2 = $list[2]-2;
	my $rangeE2 = $list[3]-2;		

	
	my $off_setF1 = $list[3]+200;							# off-sets to find downstream clusters which might belong to the gene
	my $off_setF2 = $list[3]+500;
	my $off_setF3 = $list[3]+1000;
	my $off_setF4 = $list[3]+1500;
	
	
	my $off_setR1 = $list[2]-200;							# off-sets to find downstream clusters which might belong to the gene
	my $off_setR2 = $list[2]-500;
	my $off_setR3 = $list[2]-1000;
	my $off_setR4 = $list[2]-1500;
	
	
	if($gff_id{$list[0]}{$list[1]}{$list[2]})
	{
		return ($gff_id{$list[0]}{$list[1]}{$list[2]},0);
	}
	elsif($gff_id{$list[0]}{$list[1]}{$list[3]})
	{
		return ($gff_id{$list[0]}{$list[1]}{$list[3]},0);
	}
	#elsif($gff_id{$list[0]}{$list[1]}{$rangeS1})
	#{
	#	return ($gff_id{$list[0]}{$list[1]}{$rangeS1},0);
	#}
	#elsif($gff_id{$list[0]}{$list[1]}{$rangeE1})
	#{
	#	return ($gff_id{$list[0]}{$list[1]}{$rangeE1},0);
	#}
	#elsif($gff_id{$list[0]}{$list[1]}{$rangeS2})
	#{
	#	return ($gff_id{$list[0]}{$list[1]}{$rangeS2},0);
	#}
	#elsif($gff_id{$list[0]}{$list[1]}{$rangeE2})
	#{
	#	return ($gff_id{$list[0]}{$list[1]}{$rangeE2},0);
	#}

	elsif($gff_id{$list[0]}{$list[1]}{$off_setF1} and $list[1] eq "+")
	{
		return ($gff_id{$list[0]}{$list[1]}{$off_setF1},200);
	}
	elsif($gff_id{$list[0]}{$list[1]}{$off_setF2} and $list[1] eq "+")
	{
		return ($gff_id{$list[0]}{$list[1]}{$off_setF2},500);
	}
	elsif($gff_id{$list[0]}{$list[1]}{$off_setF3} and $list[1] eq "+")
	{
		return ($gff_id{$list[0]}{$list[1]}{$off_setF3},1000);
	}
	elsif($gff_id{$list[0]}{$list[1]}{$off_setF4} and $list[1] eq "+")
	{
		return ($gff_id{$list[0]}{$list[1]}{$off_setF4},1500);
	}


	elsif($gff_id{$list[0]}{$list[1]}{$off_setR1} and $list[1] eq "-")
	{
		return ($gff_id{$list[0]}{$list[1]}{$off_setR1},200);
	}
	elsif($gff_id{$list[0]}{$list[1]}{$off_setR2} and $list[1] eq "-")
	{
		return ($gff_id{$list[0]}{$list[1]}{$off_setR2},500);
	}
	elsif($gff_id{$list[0]}{$list[1]}{$off_setR3} and $list[1] eq "-")
	{
		return ($gff_id{$list[0]}{$list[1]}{$off_setR3},1000);
	}
	elsif($gff_id{$list[0]}{$list[1]}{$off_setR4} and $list[1] eq "-")
	{
		return ($gff_id{$list[0]}{$list[1]}{$off_setR4},1500);
	}

	else
	{
		return ("NoAnnotation",0);
	}
}

sub check_query	
{
	my @list = @_;	
	my $start;											# input $list[4] (chromosome), $list[3] (hit end), $list[0] (skip length)), orientation
	
	if($list[3] eq "T")
	{
		$start = $list[1]-$list[2];
	}
	else
	{
		$start = $list[1]
	}
	
	my $query = substr $chrom, $start, $list[2];
	
	return $query;
}
