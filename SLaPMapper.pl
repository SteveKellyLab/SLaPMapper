#! /usr/bin/perl -w

use Getopt::Std;

&run_program;

exit;

sub run_program	{
				&get_options;				
				&paths;
				&define_splice;
				&get_genome;
				&prepare_mapping_reference;		
				&prepare_reads;
				&find_splice_reads;
				&map_reads;		
				&read_splice_sam_file;
				&print_splice_results;
				&read_polyA_sam_file;
				&print_polyA_results;
				unless ($gffbit eq 0)
					{
					&read_gff;
					&read_trans_splice;
					&read_polyA_sites;
					}
				&tidy_up;
				}
				
sub paths	{
			$bowtie_build = '/usr/local/bowtie/bowtie2-2.1.0/bowtie2-build'; 
			$bowtie = '/usr/local/bowtie/bowtie2-2.1.0/bowtie2'; 
			$fastq_join = '/usr/local/ea-utils/ea-utils.1.1.2-537/fastq-join';
			$trim = 'java -jar /usr/local/trimmomatic/Trimmomatic-0.30/trimmomatic-0.30.jar';
			$adapt = '/usr/local/trimmomatic/Trimmomatic-0.30/adapters/all_adaptors.fasta';
			$threads = 10;
			$cut = 20;
			$filterA = 1;
			}				
				
sub known_splice	{        												          
					my $p_serpens =                'CTATTCTAGATACAGTTTCTGTACTTTATTG';
					my $p_davidii =            'ACGCTAAAAATTGTTACAGTTTCTGTACATTATTG';
					my $t_borrelli =              'GCTATAAAAGTCACAGTTTCTGTACTTTATTG';
					my $t_carassii =   'TTGTTGTTGTTGTTATTATTGATACAGTTTCTGTACTATATTG';
					my $t_theileri =                   'TATTGATACAGTTTCTGTACTATATTG';
					my $e_gracilis = 			      		    'AGTGTCTATTTTTTTTCG';
					my $t_brucei = 			   'AACGCTATTATTAGAACAGTTTCTGTACTATATTG';
					my $l_mexicana =       'AACTAACGCTATATAAGTATCAGTTTCTGTACTTTATTG';							                                          			                     
					my $unknown = 			                         'GGCCTTGTCTGTT';
					my $splice_leader_c =					    'CGTTTCTGTACTATATTG';
					my $splice_leader_a =						'AGTTTCTGTACTATATTG';
					my $l_pyrrhocoris =    'AACTAACGCTATATAAGTATCAGTTTCTGTACTTTATTG';					                                                   
					my $perkinsela =        'AATTTCTGCTAAAATAGTTCAGTTTCTGTACTTAATTG';
					my $blechomonas_1 =    'AACTAACGCTATTATTGATACAGTTTCTGTACTTTATTG';  
					my $b_ayalai = 		   'AACTAACGCTATTATTGATACAGTTTCTGTACTATATTG';                              
					}
				
sub get_genome	{
				my %tmp = &SlurpFasta($genome_file);
				foreach my $a (keys(%tmp))
					{
					#$a =~ m/(.+?)\s/;
					#my $n = $1;
					my $n = $a;
					$tmp{$a} =~ tr/a-z/A-Z/;
					$seqs{$n} = $tmp{$a};
					}				
				}
				
sub prepare_mapping_reference	{
								system "$bowtie_build $genome_file mapping_reference";
								}
								
sub map_reads	{
				system "$bowtie -x mapping_reference -p $threads -U spliced_leader_reads.fq -S spliced_leader.sam";
				system "$bowtie -x mapping_reference -p $threads -U polyA_reads.fq -S polyA.sam";				
				}
								
sub prepare_reads	{					
					my %count = ();
					my $i = 0;					
					open FILE, $reads1;
					while (<FILE>)
						{
						my $n = $_;
						my $s = <FILE>;
						my $b = <FILE>;
						my $q = <FILE>;
						chomp $q;
						my @t = split (//, $q);
						foreach my $a (@t)
							{
							++$count{ord($a)};		
							}
						if ($i == 10000)
							{
							last;
							}
						++$i;
						}
					close FILE;

					my $encoding = '-phred64';
					
					$print_encoding = "phred64";
					
					my @sort = sort {$a <=> $b} keys %count;
					if ($sort[0] < 60)
						{
						$encoding = '-phred33';
						$print_encoding = "phred33";
						}
						
					system "$trim PE -threads $threads $encoding $reads1 $reads2 trimmed\_$reads1 unpaired_$reads1 trimmed_$reads2 unpaired_$reads2 ILLUMINACLIP\:$adapt\:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:5:15 MINLEN:50";
					system "$fastq_join trimmed\_$reads1 trimmed_$reads2 -o prepared_input_reads\.";
								
					my $j = 1;
					my @files = ("prepared_input_reads.join", "prepared_input_reads.un1", "prepared_input_reads.un2", "unpaired\_$reads1", "unpaired\_$reads2");
					open FILE, ">filtered_mapping_read_set.fastq";
					foreach my $f (@files)
						{
						open IN, $f;
						while (<IN>)
							{
							my $n = $_;
							chomp $n;
							my $s = <IN>;
							chomp $s;
							my $b = <IN>;
							my $q = <IN>;
							chomp $q;
							print FILE "\@read_$j\/1\n$s\n\+\n$q\n";
							++$j;
							}
						close IN;
						}
					close FILE;
					}
					
sub find_splice_reads	{
						open FILE, 'filtered_mapping_read_set.fastq';					
						open OUT1, ">spliced_leader_reads.fq";
						open OUT2, ">polyA_reads.fq";
						open OUT3, ">splice_and_polyA_cleaned_reads.fastq";					
						while (<FILE>)
							{
							my $n = $_;
							chomp $n;
							my $s = <FILE>;
							chomp $s;
							my $b = <FILE>;
							my $q = <FILE>;
							chomp $q;
							my @res1 = &splice_check($s, $q);
							my $p = 0;
							if ($res1[0] ne 0)
								{
								print OUT1 "$n\n$res1[0]\n\+\n$res1[1]\n";
								print OUT3 "$n\n$res1[0]\n\+\n$res1[1]\n";
								$p = 1;
								}															
							my @res3 = &polyAT_check($s, $q);
							if ($res3[0] ne 0)
								{
								my $n1 = substr($n, 0, -2);
								my $n2 = substr($n, -2, 2);
								print OUT2 "$n1\_$res3[2]$n2\n$res3[0]\n\+\n$res3[1]\n";
								unless ($p == 1)
									{
									print OUT3 "$n\n$res3[0]\n\+\n$res3[1]\n";
									$p = 1;
									}
								}
							if ($p == 0)
								{
								print OUT3 "$n\n$s\n\+\n$q\n";
								}							
							}
						close FILE;
						close OUT1;
						close OUT2;						
						}
						
sub polyAT_check	{
					my $seq = $_[0];
					my $qual = $_[1];
					my $block = 0;
					my $fail = 0;
					if ($seq =~ m/^(.+?)(A{5,})$/)
						{						
						my $frag = $1;
						my $aend = length ($2);
						my $l = length $frag;
						if ($l > $cut)
							{
							my $q = substr($qual, -$l);
							return ($frag, $q, $aend);
							}
						$block = 1;
						}
					unless ($block == 1)
						{
						$seq = reverse $seq;
						$seq =~ tr/ACGT/TGCA/;
						if ($seq =~ m/^(.+?)(A{5,})$/)
							{
							my $frag = $1;
							my $aend = length ($2);
							my $l = length $frag;
							if ($l > $cut)
								{
								my $q = reverse $qual;
								$q = substr($q, -$l);
								return ($frag, $q, $aend);
								}
							}
						}
					return ($fail, $fail, $fail);
					}						
						
sub splice_check	{
					my $seq = $_[0];
					my $qual = $_[1];
					my $block = 0;
					my $fail = 0;					
					if ($seq =~ m/$splice_frag(.+)/)
						{
						my $frag = $1;
						unless ($frag =~ m/^GTATGAGAAGCTCCC/)
							{								
							my $l = length $frag;
							if ($l > $cut)
								{
								my $q = substr($qual, -$l);
								return ($frag, $q);
								}
							}
						$block = 1;
						}
					unless ($block == 1)
						{
						$seq = reverse $seq;
						$seq =~ tr/ACGT/TGCA/;
						if ($seq =~ m/$splice_frag(.+)/)
							{
							my $frag = $1;
							unless ($frag =~ m/^GTATGAGAAGCTCCC/)
								{
								my $l = length $frag;
								if ($l > $cut)
									{
									my $q = reverse $qual;
									$q = substr($q, -$l);
									return ($frag, $q);
									}
								}
							}
						}
					return ($fail, $fail);
					}				
				
sub read_splice_sam_file	{
							my @ref = split('', $splice_frag);													
							%all_splice_data = ();
							%splice_sites = ();									
							open FILE, 'spliced_leader.sam';
							while (<FILE>)
								{
								my $l = $_;
								chomp $l;
								my @t = split(/\t/, $l);
								if ((exists $t[10])&&($t[2] ne '*'))
									{
									my $len = length $t[9];
									my $cig = $t[5];
									my $mapq = $t[4];
									$cig =~ m/(\d+)M/;
									my $ml = $1;
									if (($ml == $len)&&($mapq > 0))
										{
										if ($t[1] == 0)
											{									
											my $s2 = substr ($seqs{$t[2]}, $t[3] - 1 - $alen, $alen);	
											my $s3 = substr ($seqs{$t[2]}, $t[3] - 3, 2);							
											my @ts2 = split ('', $s2);
											my $i = 0;
											my $s = 0;
											if ($#ref == $#ts2)
												{
												while ($i <= $#ref)
													{
													if ($ref[$i] eq $ts2[$i])
														{
														++$s;
														}
													++$i;
													}
												unless ($s >= $alen - 2)
													{
													++$all_splice_data{$t[2]}{$t[3] - 1}{'forward'};
													$splice_sites{$t[2]}{$t[3] - 1} = $s3;
													}
												}
											}
										elsif ($t[1] == 16)
											{									
											my $s2 = substr ($seqs{$t[2]}, $t[3] - 1 + $len, $alen);									
											my $s3 = substr ($seqs{$t[2]}, $t[3] - 1 + $len, 2);
											$s3 = reverse $s3;											
											$s3 =~ tr/ACGT/TGCA/;
											$s2 = reverse $s2;
											$s2 =~ tr/ACGT/TGCA/;
											my @ts2 = split ('', $s2);									
											my $i = 0;
											my $s = 0;
											if ($#ref == $#ts2)
												{
												while ($i <= $#ref)
													{
													if ($ref[$i] eq $ts2[$i])
														{
														++$s;
														}
													++$i;
													}
												unless ($s >= $alen - 2)
													{
													++$all_splice_data{$t[2]}{$t[3] - 1 + $len}{'reverse'};
													$splice_sites{$t[2]}{$t[3] - 1 + $len} = $s3;
													}
												}									
											}
										}
									}
								}
							close FILE;					
							}

sub print_splice_results	{
							open FILE, ">prelim_trans_splice_sites.xls";
							print FILE "Chromosome\tposition\tCount\tStrand\tdinucleotide\n";
							$fc = 0;
							$rc = 0;
							my @s1 = sort {$a cmp $b} keys %all_splice_data;
							foreach my $c (@s1)
								{
								my @temp = keys %{$all_splice_data{$c}};
								my @s2 = sort {$a <=> $b} @temp;
								foreach my $p (@s2)
									{
									if (exists $all_splice_data{$c}{$p}{'forward'})
										{
										print FILE "$c\t$p\t$all_splice_data{$c}{$p}{'forward'}\t\+\t$splice_sites{$c}{$p}\n";
										++$fc;
										}
									if (exists $all_splice_data{$c}{$p}{'reverse'})
										{
										print FILE "$c\t$p\t$all_splice_data{$c}{$p}{'reverse'}\t\-\t$splice_sites{$c}{$p}\n";
										++$rc;
										}
									}
								}
							close FILE;	
							print "$fc forward sites $rc reverse sites\n";				
							}
							
sub read_polyA_sam_file	{
						%all_polya_data = ();									
						open FILE, 'polyA.sam';
						open OUT, ">accepted_polyA_reads.sam";						
						$real_polyA_site = 0;
						while (<FILE>)
							{
							my $l = $_;
							chomp $l;
							my @t = split(/\t/, $l);
							if (exists $t[10])
								{
								if ($t[2] ne '*')
									{
									my $len = length $t[9];
									my $cig = $t[5];
									my $mapq = $t[4];
									$cig =~ m/(\d+)M/;
									my $ml = $1;
									$t[0] =~ m/\_(\d+)\//;
									my $allen = $1;
									if (($ml == $len)&&($mapq > 0))
										{
										if ($t[1] == 0)
											{
											my $s2 = substr ($seqs{$t[2]}, $t[3] - 1 + $len, $allen);									
											my $ds = ($s2 =~ tr/A//);
											if ($filterA == 0)
												{
												$ds = 0;
												}
											if ($ds < $allen/2)
												{
												++$all_polya_data{$t[2]}{$t[3] - 1 + $len}{'forward'};
												++$real_polyA_site;
												print OUT "$l\n";
												} 
											}
										elsif ($t[1] == 16)
											{										
											my $s2 = substr ($seqs{$t[2]}, $t[3] - 1 - $allen, $allen);									
											my $ds = ($s2 =~ tr/T//);
											if ($filterA == 0)
												{
												$ds = 0;
												}
											if ($ds < $allen/2)
												{
												++$all_polya_data{$t[2]}{$t[3] - 1}{'reverse'};
												++$real_polyA_site;
												print OUT "$l\n";
												} 
											}
										}
									}
								}
							else
								{
								print OUT "$l\n";
								}
							}
						close FILE;	
						close OUT;				
						}

sub print_polyA_results	{
						open FILE, ">polyA_sites.xls";				
						print FILE "Chromosome\tPosition\tCount\tStrand\n";
						my @s1 = sort {$a cmp $b} keys %all_polya_data;
						foreach my $c (@s1)
							{
							my @temp = keys %{$all_polya_data{$c}};
							my @s2 = sort {$a <=> $b} @temp;
							foreach my $p (@s2)
								{
								if (exists $all_polya_data{$c}{$p}{'forward'})
									{
									print FILE "$c\t$p\t$all_polya_data{$c}{$p}{'forward'}\t\+\n";
									}
								if (exists $all_polya_data{$c}{$p}{'reverse'})
									{
									print FILE "$c\t$p\t$all_polya_data{$c}{$p}{'reverse'}\t\-\n";
									}
								}
							}
						close FILE;					
						}
						
sub read_gff	{
				open FILE, $gff;				
				while (<FILE>)
					{
					my $l = $_;
					chomp $l;
					my @t = split(/\t/, $l);
					if (exists $t[8])
						{
						my $c = $t[0];
						$c =~ s/v5\.1/v5/gmis;
						my $s = $t[3];
						my $e = $t[4];
						if ($t[2] eq 'CDS')
							{
							my $a = ();
							if ($t[8] =~ m/\;/)
								{
								$t[8] =~ m/ID\=(.+?)\;/;
								$a = $1;
								}
							else
								{
								$t[8] =~ m/ID\=(.+)/;
								$a = $1;
								}							
							if ($a =~ s/cds\_//)
								{
								$a =~ s/\-\d+//;
								unless (exists $all_data{$c}{$a})
									{									
									$all_data{$c}{$a}{'start'} = $s;
									$all_data{$c}{$a}{'end'} = $e;
									$all_data{$c}{$a}{'strand'} = $t[6];
									}
								}
							}
						}
					}
				close FILE;
				foreach my $c (keys(%all_data))
					{
					@$c = sort {$all_data{$c}{$a}{'start'} <=> $all_data{$c}{$b}{'start'}} keys %{$all_data{$c}};
					}
				}
				
sub read_trans_splice	{
						open FILE, "prelim_trans_splice_sites.xls";
						open OUT, ">trans_splice_sites.xls";
						open OUT2, ">trans_splice_sites.bed";
						open OUT3, ">trans_splice_sites.gff";
						print OUT "Chromosome\tPosition\tDinucleotide\tSite-strand\tClosest CDS\tCDS Strand\tDistance to CDS\tStatus\tCount\n";
						print OUT2 "track name=\"Trans-splice sites\" visibility=2\n";
						my $i = 0;						
						my $j = 1;
						while (<FILE>)
							{
							my $l = $_;
							chomp $l;
							my @t = split(/\t/, $l);
							unless($i == 0)
								{
								if (exists $t[2])
									{
									my $c = $t[0];
									my $ps = $t[1] - 2;
									my $pe = $ps + 2;
									
									my $psgff = $t[1] - 1;
									my $pegff = $psgff + 1;
									
									my $count = $t[2];
									my $strand = $t[3];
									if ($strand eq '-')
										{
										$ps = $t[1];
										$pe = $ps + 2;
										$psgff = $t[1] + 1;
										$pegff = $psgff + 1;
										}
									my $min = 1000000000000;
									my $rem = 'undefined';
									my $block = 0;
									foreach my $a (@$c)
										{
										if (($ps > $all_data{$c}{$a}{'start'})&&($ps < $all_data{$c}{$a}{'end'}))
											{
											$block = 1;
											$rem = $a;
											last;
											}
										else
											{
											if (($t[3] eq '+')&&($all_data{$c}{$a}{'strand'} eq '+'))
												{
												my $d = sqrt(($ps - $all_data{$c}{$a}{'start'})**2);
												if (($d < $min)&&($ps < $all_data{$c}{$a}{'start'}))
													{
													$min = $d;
													$rem = $a;
													$rems = 1;
													}
												}
											elsif (($t[3] eq '-')&&($all_data{$c}{$a}{'strand'} eq '-'))
												{
												my $d = sqrt(($ps - $all_data{$c}{$a}{'end'})**2);
												if (($d < $min)&&($ps > $all_data{$c}{$a}{'end'}))
													{
													$min = $d;
													$rem = $a;
													$rems = 1;
													}
												}
											elsif (($t[3] eq '-')&&($all_data{$c}{$a}{'strand'} eq '+'))
												{
												my $d = sqrt(($ps - $all_data{$c}{$a}{'start'})**2);
												if ($d < $min)
													{
													$min = $d;
													$rem = $a;
													$rems = 0;
													}
												}
											elsif (($t[3] eq '+')&&($all_data{$c}{$a}{'strand'} eq '-'))
												{
												my $d = sqrt(($ps - $all_data{$c}{$a}{'end'})**2);
												if ($d < $min)
													{
													$min = $d;
													$rem = $a;
													$rems = 0;
													}
												}
											}
										}
									if ($rem eq 'undefined')
										{
										print OUT "$t[0]\t$psgff\t$t[4]\t$t[3]\t$rem\tNA\tNA\tNA\t$count";
										}
									elsif ($block == 1)
										{								
										print OUT "$t[0]\t$psgff\t$t[4]\t$t[3]\t$rem\t$all_data{$c}{$rem}{'strand'}\t0\tInternal\t$count";
										}
									elsif ($rems == 1)
										{
										print OUT "$t[0]\t$psgff\t$t[4]\t$t[3]\t$rem\t$all_data{$c}{$rem}{'strand'}\t$min\tOK\t$count";
										}	
									elsif ($rems == 0)
										{
										print OUT "$t[0]\t$psgff\t$t[4]\t$t[3]\t$rem\t$all_data{$c}{$rem}{'strand'}\tNA\tOpposite\t$count";
										}									
									print OUT2 "$t[0]\t$ps\t$pe\tSAS\t$count\t$t[3]\n";
									print OUT3 "$t[0]\tSLaP\tSAS\t$psgff\t$pegff\t.\t$t[3]\t\.\tID\=SAS\_$j\;count\=$count\;\n";
									print OUT "\n";									
									++$j;
									}
								}
							$i = 1;
							}
						close FILE;
						close OUT;
						close OUT2;
						close OUT3;
						}
						
sub read_polyA_sites	{
						open FILE, "polyA_sites.xls";
						open OUT, ">gene_linked_polyA_sites.xls";
						open OUT3, ">polyA_sites.gff";
						print OUT "Chromosome\tPosition\tSite strand\tClosest CDS\tCDS Strand\tDistance to CDS\tStatus\tCount\n";
						my $i = 0;
						my $j = 1;
						while (<FILE>)
							{
							my $l = $_;
							chomp $l;
							my @t = split(/\t/, $l);
							unless ($i == 0)
								{
								if (exists $t[2])
									{
									my $c = $t[0];
									
									my $ps = $t[1] - 1;
									my $pe = $ps + 1;
									
									my $psgff = $t[1];
									my $pegff = $psgff;
																		
									my $strand = $t[3];
									if ($strand eq '-')
										{
										$ps = $t[1];
										$pe = $ps + 1;
										$psgff = $t[1] + 1;
										$pegff = $psgff;
										}
									my $count = $t[2];									
									my $min = 1000000000000;
									my $rem = 'undefined';
									my $block = 0;
									foreach my $a (@$c)
										{
										if (($ps > $all_data{$c}{$a}{'start'})&&($ps < $all_data{$c}{$a}{'end'}))
											{
											$block = 1;
											$rem = $a;
											last;
											}
										else
											{
											if (($strand eq '+')&&($all_data{$c}{$a}{'strand'} eq '+'))
												{
												my $d = sqrt(($ps - $all_data{$c}{$a}{'end'})**2);
												if (($d < $min)&&($ps > $all_data{$c}{$a}{'end'}))
													{
													$min = $d;
													$rem = $a;
													$rems = 1;
													}
												}
											elsif (($strand eq '-')&&($all_data{$c}{$a}{'strand'} eq '-'))
												{
												my $d = sqrt(($ps - $all_data{$c}{$a}{'start'})**2);
												if (($d < $min)&&($ps < $all_data{$c}{$a}{'start'}))
													{
													$min = $d;
													$rem = $a;
													$rems = 1;
													}
												}
											elsif (($strand eq '-')&&($all_data{$c}{$a}{'strand'} eq '+'))
												{
												my $d = sqrt(($ps - $all_data{$c}{$a}{'end'})**2);
												if ($d < $min)
													{
													$min = $d;
													$rem = $a;
													$rems = 0;
													}
												}
											elsif (($strand eq '+')&&($all_data{$c}{$a}{'strand'} eq '-'))
												{
												my $d = sqrt(($ps - $all_data{$c}{$a}{'start'})**2);
												if ($d < $min)
													{
													$min = $d;
													$rem = $a;
													$rems = 0;
													}
												}
											}
										}
									if ($rem eq 'undefined')
										{
										print OUT "$t[0]\t$psgff\t$t[3]\t$rem\tNA\tNA\tNA\t$count";
										}
									elsif ($block == 1)
										{								
										print OUT "$t[0]\t$psgff\t$t[3]\t$rem\t$all_data{$c}{$rem}{'strand'}\t0\tInternal\t$count";
										}
									elsif ($rems == 1)
										{
										print OUT "$t[0]\t$psgff\t$t[3]\t$rem\t$all_data{$c}{$rem}{'strand'}\t$min\tOK\t$count";
										}	
									elsif ($rems == 0)
										{
										print OUT "$t[0]\t$psgff\t$t[3]\t$rem\t$all_data{$c}{$rem}{'strand'}\tNA\tOpposite\t$count";
										}	
									print OUT3 "$t[0]\tSLaP\tPAS\t$psgff\t$pegff\t.\t$t[3]\t\.\tID\=PAS\_$j\;count\=$count\;\n";
									print OUT "\n";	
									}									
								++$j;
								}
							$i = 1;
							}
						close FILE;
						close OUT;
						close OUT3;
						}				

					
sub define_splice	{
					$splice_frag = substr($spliced_leader_seq, -$alen);					
					}
					
sub SlurpFasta	{
				my %seqs = ();
				my @tmp = ();
				my $i = 0;
				my $acc = ();
				my $next = ();
				open FILE, $_[0];
				while (<FILE>)
					{
					my $l = $_;
					if ($l =~ m/^>(\S+)/)
						{
						$next = $1;
						if ($i == 0)
							{
							$acc = $next;
							}
						else
							{
							$seqs{$acc} = &process_seq(@tmp);
							$acc = $next;
							@tmp = ();
							}						
						$i = 1;						
						}
					else
						{
						push (@tmp, $l);
						}
					}
				close FILE;
				if (defined $acc)
					{
					$seqs{$acc} = &process_seq(@tmp);
					}
				return %seqs;
				}
					
sub process_seq	{
				my $j = join '', @_;
				$j =~ s/[\s\n]//g;
				return $j;	
				}
			
sub tidy_up	{
			my @files = qw/prepared_input_reads.join prepared_input_reads.un1 prepared_input_reads.un2 filtered_mapping_read_set.fastq spliced_leader_reads.fq polyA_reads.fq mapping_reference.1.bt2 mapping_reference.2.bt2 mapping_reference.3.bt2 mapping_reference.4.bt2 mapping_reference.rev.1.bt2 mapping_reference.rev.2.bt2/;
			push (@files, "unpaired\_$reads1");
			push (@files, "unpaired\_$reads2");
			push (@files, "trimmed\_$reads1");
			push (@files, "trimmed\_$reads2");
			foreach my $f (@files)
				{
				if (-e $f)
					{
					system "rm $f";
					}					
				}			
			}
			
sub exit_screen	{
				my $line = "_______________________________________________\n";
				system "clear";
				print $line;
				print "Steve's spliced leader and polyA finder program\n";
				print $line;
				print "MANDATORY INPUT FILES\n\n";
				print "-g Genome fasta file\n";
				print "-l Paired end read file 1\n";
				print "-r Paired end read file 2\n";
				print "-a GFF file for genome FASTA\n";
				print $line;
				print "OTHER OPTIONS\n\n";
				print "-i Spliced leader sequence\n";
				print "-s Minimum spliced leader length\n";
				print $line;
				print "\n";
				exit;
				}
			
sub get_options	{
				&getopts('g:l:r:t:s:i:a:',\%parameters);
				if (exists $parameters{'g'})
					{
					$genome_file = $parameters{'g'};
					chomp $genome_file;
					}
				else
					{
					&exit_screen;
					}
				if (exists $parameters{'l'})
					{
					$reads1 = $parameters{'l'};
					chomp $reads1;
					}
				else
					{
					&exit_screen;
					}
				if (exists $parameters{'a'})
					{
					$gff = $parameters{'a'};
					chomp $gff;
					$gffbit = 1;
					}
				else
					{
					$gffbit = 0;
					}
				if (exists $parameters{'r'})
					{
					$reads2 = $parameters{'r'};
					chomp $reads2;
					}
				else
					{
					&exit_screen;
					}
				if (exists $parameters{'s'})
					{
					$alen = $parameters{'s'};
					chomp $alen;
					}
				else
					{
					$alen = 12;
					}
				if (exists $parameters{'i'})
					{
					$spliced_leader_seq = $parameters{'i'};
					chomp $spliced_leader_seq;
					}
				else
					{					                        
					print "you must specify a spliced leader sequence eg AACGCTATTATTAGAACAGTTTCTGTACTATATTG\n";
					exit;
					}				
				}
				


