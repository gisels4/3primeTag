#!/usr/bin/env perl
###################################################################################################################
# 	
#	analyze_RNAhybrid_V5.pl
#
#	takes the compact output of RNAhybrid and filters for p-value and recalculate the position on the genome of 
#	the potential target
#	Further it adds the annotations to the folding hits
#
# 	Andreas Gisel January 2012
#			modified Andreas Gisel September 2017
#
###################################################################################################################

use strict;
# use Bio::SeqIO;

my $rnahybrid_data = shift;		#RNAhybrid output
my $query_seq = shift;	# fasta input for RNAhybrid
my $tresh_p_value = shift;	# folding p-value
my $tresh_score = shift;	# Carrington score
my $abund = shift;		# abundance for the smallRNA
my $file_out = shift;
my $file_anno = shift;	# annotation file such as 
my $file_gff = shift;	# GFF for the annotation of the RNAhybrid hits

my $file_log = "log.txt";
unless(open(LOG, ">$file_log"))
{
	print "cannot open $file_out!!\n";
	exit;
}

my $time0 = time;
my $time1 = time;
my $time = $time1-$time0;
print LOG "Start SeqIO - $time\n";

# test the changes


# my $seq_in  = Bio::SeqIO->new( -format => 'fasta',
#                              -file => $query_seq);
                               
my %query_seq;

#while( my $seq = $seq_in->next_seq() ) 
#{
#	$query_seq{$seq->id} = $seq->seq;
#}

unless(open(SEQ, $query_seq))
{
	print "cannot open $query_seq!!\n";
	exit;
}

my $id_seq;
my $seq;

while(<SEQ>)
{
	if(/>(.*)/)
	{
		$id_seq = $1;
		chomp $id_seq;
		$seq = <SEQ>;
		chomp $seq;
		$query_seq{$id_seq} = $seq;
	}

}

$time1 = time;
$time = $time1-$time0;
print LOG "end SeqIO - $time\n";

#########################
# read gff file 

print LOG "start gff read \n";

unless(open(IN, $file_gff))
{
	print LOG "cannot open $file_gff!!\n";
	exit;
}

my $id;		# gene ID
my @chr;	# chromosome
my @feat;	# feature
my @start;	# feature start
my @end;	# feature end
my @strand;	# strand orientation
my @xyz;		# gene ID

while(<IN>)	
{	
	next if(/##/);
	# Chromosome01	phytozomev10	gene	4	803	.	+	.	ID=Manes.01G000100.v6.1;Name=Manes.01G000100
	my @data = split /\t/;

	# ID=Manes.08G002000.v6.1;Name=Manes.08G002000
	# Chr1	TAIR10	gene	3631	5899	.	+	.	ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010
	# Chr1	TAIR10	exon	3996	4276	.	+	.	Parent=AT1G01010.1
	
	if($data[2] =~ /gene/)		# calculate promoter and terminal regions relative to the gene position
	{
		$data[8] =~ /ID=(.*)\;Note/;
		$id = $1;
		$id = "$id.1";
		
		#chomp $data[8];
		#print "A $data[0] - $data[2] - $data[3] - $data[8] = $id\n";
		
		if($data[6] eq "+")
		{
			push @chr,$data[0];				# promoter regions from 0 to 1500
			push @feat, "prom-200";
			push @start, $data[3]-201;
			push @end, $data[3]-1;
			push @strand, $data[6];
			push @xyz, $id;
			
			push @chr,$data[0];
			push @feat, "prom-500";
			push @start, $data[3]-501;
			push @end, $data[3]-201;
			push @strand, $data[6];
			push @xyz, $id;
			
			push @chr,$data[0];
			push @feat, "prom-1000";
			push @start, $data[3]-1001;
			push @end, $data[3]-501;
			push @strand, $data[6];
			push @xyz, $id;
			
			push @chr,$data[0];
			push @feat, "prom-1500";
			push @start, $data[3]-1501;
			push @end, $data[3]-1001;
			push @strand, $data[6];
			push @xyz, $id;
			
			push @chr,$data[0];				# terminal regions from 0 to 1500
			push @feat, "term-200";
			push @start, $data[4]+1;
			push @end, $data[4]+201;
			push @strand, $data[6];
			push @xyz, $id;
			
			push @chr,$data[0];
			push @feat, "term-500";
			push @start, $data[4]+201;
			push @end, $data[4]+501;
			push @strand, $data[6];
			push @xyz, $id;
			
			push @chr,$data[0];
			push @feat, "term-1000";
			push @start, $data[4]+501;
			push @end, $data[4]+1001;
			push @strand, $data[6];
			push @xyz, $id;
			
			push @chr,$data[0];
			push @feat, "term-1500";
			push @start, $data[4]+1001;
			push @end, $data[4]+1501;
			push @strand, $data[6];
			push @xyz, $id;
		}
		elsif ($data[6] eq "-")
		{
			push @chr,$data[0];
			push @feat, "prom-200";
			push @start, $data[4]+201;
			push @end, $data[4]+1;
			push @strand, $data[6];
			push @xyz, $id;
			
			push @chr,$data[0];
			push @feat, "prom-500";
			push @start, $data[4]+501;
			push @end, $data[4]+201;
			push @strand, $data[6];
			push @xyz, $id;
			
			push @chr,$data[0];
			push @feat, "prom-1000";
			push @start, $data[4]+1001;
			push @end, $data[4]+501;
			push @strand, $data[6];
			push @xyz, $id;
			
			push @chr,$data[0];
			push @feat, "prom-1500";
			push @start, $data[4]+1501;
			push @end, $data[4]+1001;
			push @strand, $data[6];
			push @xyz, $id;
			
			push @chr,$data[0];			# terminal regions from 0 to 1500
			push @feat, "term-200";
			push @start, $data[3]-1;
			push @end, $data[3]-201;
			push @strand, $data[6];
			push @xyz, $id;
			
			push @chr,$data[0];
			push @feat, "term-500";
			push @start, $data[3]-201;
			push @end, $data[3]-501;
			push @strand, $data[6];
			push @xyz, $id;
			
			push @chr,$data[0];
			push @feat, "term-1000";
			push @start, $data[3]-501;
			push @end, $data[3]-1001;
			push @strand, $data[6];
			push @xyz, $id;
			
			push @chr,$data[0];
			push @feat, "term-1500";
			push @start, $data[3]-1001;
			push @end, $data[3]-1501;
			push @strand, $data[6];
			push @xyz, $id;
		}
		else
		{
			print "HUCH! \n";
			exit;
		}
	}
	
	if($data[2] !~ /chromosome/)
	{
		my @id = split/;/, $data[8];
		
		if($id[0] =~ /=(.*)\.\d,/)
		{
			$id = $1;
			$id = "$id.1";			
		}
		elsif($id[0] =~ /=(.*)\.\d/)
		{
			$id = $1;
			$id = "$id.1";		
		}
		elsif($id[0] =~ /=(.*)/)
		{
			$id = $1;
			$id = "$id.1";			
		}
		
		if ($data[8] =~ /Alias=(.*)/)
		{
			$id = "$id\t$1";
		}
		
		
		#chomp $data[8];
		#print "B $data[0] - $data[2] - $data[3] - $data[8] = $id\n";

		

		#print "B:$id\n";
		#print "BB:$data[0] - $data[2] - $data[3] - $data[4] - $data[6] - $id\n";
		
		push @chr,$data[0];
		push @feat, $data[2];
		push @start, $data[3];
		push @end, $data[4];
		push @strand, $data[6];
		push @xyz, $id;
	}
}

close IN;

$time1 = time;
$time = $time1-$time0;
print LOG "end gff - $time\n";

my %anno;
if ($file_anno)
{
	unless(open(ANNO, $file_anno))
	{
		print "cannot open $file_anno!!\n";
		exit;
	}
	
	while(<ANNO>)
	{
		my @data = split /\t/;
		chomp $_;
		# ppa000001m	PF07728,PF07726	PTHR22908	KOG1808			AT1G67120.1		ATP binding / ATPase/ nucleoside-triphosphatase/ nucleotide binding / transcription factor binding
		# 32359217	Manes.01G000300	Manes.01G000300.1	Manes.01G000300.1.p	PF07851	PTHR21433	KOG4758			GO:0016021	AT4G10430.3		TMPIT-like protein
		# AT1G01010.1	protein_coding	NAC domain containing protein 1		NAC domain containing protein 1 (NAC001); FUNCTIONS IN: sequence-specific DNA binding transcription factor activity; INVOLVED IN: multicellular organismal development, regulation of transcription; LOCATED IN: cellular_component unknown; EXPRESSED IN: 7 plant structures; EXPRESSED DURING: 4 anthesis, C globular stage, petal differentiation and expansion stage; CONTAINS InterPro DOMAIN/s: No apical meristem (NAM) protein (InterPro:IPR003441); BEST Arabidopsis thaliana protein match is: NAC domain containing protein 69 (TAIR:AT4G01550.1); Has 2503 Blast hits to 2496 proteins in 69 species: Archae - 0; Bacteria - 0; Metazoa - 0; Fungi - 0; Plants - 2502; Viruses - 0; Other Eukaryotes - 1 (source: NCBI BLink).

		$anno{$data[0]} = $_;
	}
}


$time1 = time;
$time = $time1-$time0;
print LOG "end annotation - $time\n";

unless(open(RNA, $rnahybrid_data))
{
	print LOG "cannot open $rnahybrid_data!!\n";
	exit;
}

my $flag;

unless(open(OUT, ">$file_out"))
{
	print LOG "cannot open $file_out!!\n";
	exit;
}


print LOG "start RNAhybrid read \n";

while(<RNA>)
{
	my @data = split /:/;
	# Chromosome04_24813622-24814620:999:D3-318137_19_x39:19:-36.0:0.043152:11:U  C                C: CU CCAACUGAGCUAUCCC : GA GGUUGACUCGAUAGGG :   U  
	# 4_5562090-5563088:999:G3-5313925_24_35:24:-48.8:0.001474:60:U                        A: AACUUCGAUCGUAGUAGUGGGCAU : UUGAAGCUAGCAUCAUCACCCGUA :   
	
	my @pos = split /_/, $data[0];
	
	
	my $pos_seq;
	
	if($pos[1] =~ /(\d+)-/)		# get coordinate of the subject sequence fragment
	{
		$pos_seq = $1;
	}
	else
	{
		print "Wrong coordination data\n";
		exit;
	}
	
	########
	# add Chr to the chromosome information and fix Mt and Pt
	
	if($pos[0] =~ /Mt/)
	{
		$pos[0] = "ChrM";
	}
	elsif($pos[0] =~ /Pt/)
	{
		$pos[0] = "ChrC";
	}
	else
	{
		$pos[0] = "Chr$pos[0]";
	}
	

	my $pos_targ = $pos_seq+$data[6]-1;		# calculate the position of the hit
	
	my $freq;
	if($data[2] =~ /_\d+_x(\d+)/)		# extract the sequence abundance
	{
		$freq = $1;
	}
	else
	{
			print "Errore in sequence frequency\n";
	}

	
	$data[7] =~ tr/ /*/;
	$data[8] =~ tr/ /*/;
	$data[9] =~ tr/ /*/;
	$data[10] =~ tr/ /*/;
	
	my $fold = "$data[7]\t$data[8]\t$data[9]\t$data[10]";
	my $fold1 = "$data[7]\n$data[8]\n$data[9]\n$data[10]";

#	print "$fold1";
	chomp $fold;
	
	my @fold1 = split //, $data[7];
	my @fold2 = split //, $data[8];
	my @fold3 = split //, $data[9];
	my @fold4 = split //, $data[10];
	chomp $data[10];
	my $text;
	
	if($query_seq{$data[2]})
	{
		$text = "$pos[0]\t$pos[1]\t$pos_targ\t$data[2]\t$query_seq{$data[2]}\t$data[7]\t$data[8]\t$data[9]\t$data[10]\t";
	}
	else
	{
		$text = "$pos[0]\t$pos[1]\t$pos_targ\t$data[2]\t - \t$data[7]\t$data[8]\t$data[9]\t$data[10]\t";
		print LOG "missing seq?? $data[2]\n";
	}
	
	
	#############
	# Calculate the Carrington score
	#
	
	# TTTGAGGAACAAAAGTTCTTT
	
	# A    U                A
	#  GAAG GCUUUUGUUCCUCAAA 
	#  UUUC UGAAAACAAGGAGUUU 
	# U                      
	
	# load hybrid pairs
	my $score = 0;
	my $flag = 0;
	my %pair1 = ( 'A', 'U', 'U', 'A', 'G', 'C' , 'C', 'G');
	my %pair2 = ( 'G', 'U', 'U', 'G');
	
	my $count = 0;
	for (my $i=$#fold1;$i>=0;$i--)
	{
		if($fold3[$i] eq "*" and $fold4[$i] eq "*" and ($count == 0 or $count == length($query_seq{$data[2]})))
		{
			next;
		}
		elsif($fold3[$i] eq "*" and $fold4[$i] eq "*")
		{
			if($count >= 2 and $count <= 13)
			{
				$score = $score + 2;
				$text = $text ."-2";
			}
			else
			{
				$score = $score + 1;
				$text = $text ."-1";
			}
		}
		elsif($fold4[$i] ne "*")
		{
			$count++;
			if($count >= 2 and $count <= 13)
			{
				$score = $score + 2;
				$text = $text ."-2";
			}
			else
			{
				$score = $score + 1;
				$text = $text ."-1";
			}
		}
		else
		{
			$count++;
			if ($pair1{$fold3[$i]} eq $fold2[$i])
			{
				$text = $text ."-$fold3[$i]$fold2[$i]0";
			}
			elsif($pair2{$fold3[$i]} eq $fold2[$i])
			{
				if($count >= 2 and $count <= 13)
				{
					$score = $score + 1;
					$text = $text ."-$fold3[$i]$fold2[$i]1";
				}
				else
				{
					$score = $score + 0.5;
					$text = $text ."-$fold3[$i]$fold2[$i]0.5";
				}
			}
			else
			{
				if($count >= 2 and $count <= 13)
				{
					$score = $score + 2;
					$text = $text ."-$fold3[$i]$fold2[$i]2";
				}
				else
				{
					$score = $score + 1;
					$text = $text ."-$fold3[$i]$fold2[$i]1";
				}
			}
		}
	}
	
	$text = $text ."\t$score\n";
	# print LOG "$text\n";
	
	if($query_seq{$data[2]})
	{
		if($score < $tresh_score and $data[5] < $tresh_p_value and $freq > $abund)
		{
			$flag = 0;
			for(my$i=0;$i<=$#chr;$i++)
			{	
				if($pos[0] eq $chr[$i] and $pos_targ >= $start[$i] and $pos_targ <= $end[$i])
				{
					if(!$anno{$xyz[$i]})
					{
						print OUT "$query_seq{$data[2]}\t$freq\t$data[2]\t$pos[0]-$pos[1]\t$strand[$i]\t$pos_targ\t$data[4]\t$data[5]\t$fold\t$score\t$feat[$i]\t$xyz[$i]\t no annotation\n";
					}
					else
					{
						print OUT "$query_seq{$data[2]}\t$freq\t$data[2]\t$pos[0]-$pos[1]\t$strand[$i]\t$pos_targ\t$data[4]\t$data[5]\t$fold\t$score\t$feat[$i]\t$xyz[$i]\t$anno{$xyz[$i]}\n";
					}
					$flag = 1;
				}
			}
			
			if($flag == 0)
			{
				print OUT "$query_seq{$data[2]}\t$freq\t$data[2]\t$pos[0]-$pos[1]\t---\t$pos_targ\t$data[4]\t$data[5]\t$fold\t$score\tno annotation\n";
			}
		}
	}
	else
	{
		print OUT " --- \t --- \t$data[2]\t$pos[0]-$pos[1]\t---\t$pos_targ\t$data[4]\t$data[5]\t$fold\t$score\tno annotation\n";
	}
	
}

$time1 = time;
$time = $time1-$time0;
print LOG "end RNAhybrid - $time\n";

exit;


# perl analyze_RNAhybridV5.pl file_9_out.txt all_exp-n_semi_nonzero_reg_Unique.fasta 0.05 5 10 newtest.txt Mesculenta_305_v6.1.annotation_info.txt Mesculenta_305_v6.1.gene_exons.gff3


