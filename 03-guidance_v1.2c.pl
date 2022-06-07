## 20 February 2015
#!/usr/bin/perl

##########################################################################
# version 1.2c
# Adapted to Katz GUIDANCE outputs
# Within a folder there is one folder for each FASTA GUIDANCE RUN output
#
# The FASTA file will be unaligned 
#  and it will not have masked residues nor columns removed, only low scored sequences excluded
#
# The PHYLIP file will be aligned 
#  and filtered using the specified thresholds for sequence, column and residue scores.  
##########################################################################


use strict;
use warnings;
use Getopt::Long;


# Script input options

my ($filename,$dir,$outdir,$maskcutoff,$columncutoff,$datatype,$prefix, $seqcutoff);
my $getoptResult = GetOptions 	("filename=s"   => \$filename,  # by default the aln is named MSA.MAFFT.aln.With_Names
								 "inDir=s"      => \$dir, # the directory with all subdirectories of each GUIDANCE RUN
								 "outDir=s"     => \$outdir, # a new folder were to save the filtered files (FASTA and PHYLIP) 
								 "siteCutoff=s" => \$maskcutoff, # residue guidance score cutoff below which residue will be masked (only for PHYLIP)
								 "colCutoff=s"  => \$columncutoff, # column guidance score cutoff below which they will be excluded from final aln (only for PHYLIP)
								 "seqCutoff=s"  => \$seqcutoff,	# sequence guidance score cutoff below which they will be excluded from final files ((FASTA and PHYLIP)						 
								 "dataType=s"   => \$datatype, # datatype could be [DNA|RNA] or AA
								 "prefix=s"     => \$prefix, # prefix to be added for the output files
								 );
								 
								die "ERROR: siteCutoff must be a number\n" if ($maskcutoff !~/^\d\.?(\d+)?$/);
								die "ERROR: colCutoff must be a number\n" if ($columncutoff !~/^\d\.?(\d+)?$/);
								die "ERROR: seqCutoff must be a number\n" if ($seqcutoff !~/^\d\.?(\d+)?$/); 
								if (($datatype ne "nuc") and ($datatype ne "codon") and ($datatype ne "AA")){
									die "ERROR: dataType should be nuc, codon or AA (case sensitive)\n";
								}
								die "ERROR: prefix must not be an empty string\n" if ($prefix eq "");
								 
#$filename =~ s/\.aln//;

# Variables that will store the numbers of initial and final sequences, columns in aln and non masked residues
# this information is printed to STDOUT
	
	print "here is the name $filename";
	my ($InitialSeq,$FinalSeq);
	my ($InitialCol,$FinalCol);
	my ($InitialRes,$FinalRes);


sub MakeTempAlignment (){

	print "here we are before Make";
	## This function will reformat the aligned FASTA file (*.MAFFT).aln.With_Names
	## creating a temporary aln in FASTA format
	## each line will have one sequence
	
	system ("cp ${dir}/${filename}.aln.With_Names ${dir}/${filename}.temp");
#	system ("grep '^>' ${dir}/Seqs.Orig*fas | perl -p -e 's/^>//g' > ${dir}/seq_names"); 
	system ("grep '^>' ${dir}/Seqs.Orig*fas | perl -p -e 's/^>//g' > ${dir}/seq_names");   # MACR corrected here, "i" argument requires input file
	system("perl -pi -e 's/(>.+)\n\$/!\$1!\n/g' ${dir}/${filename}.temp"); # add a "!" at the beginning and end of sequence name
#	system("perl -pi -e 's/\n\$//g' ${dir}/${filename}.temp"); # remove end of lines "\n"	
	system("perl -pi -e 's/\n//g' ${dir}/${filename}.temp"); # remove end of lines "\n"      # MACR --- It wasn't working properly. So, I modified here
	system("perl -pi -e 's/!>/\n>/g' ${dir}/${filename}.temp");
	system("perl -pi -e 's/!/\n/g' ${dir}/${filename}.temp");

	# keeping only the sequences
	system("grep -v '^>' ${dir}/${filename}.temp > ${dir}/${filename}.temp2"); 
#	system("awk 'NR > 1' ${dir}/${filename}.temp | awk '{split(\$0,a,\"!\") ; print a[2]}' > ${dir}/${filename}.temp2"); 
	system("grep -v '^\$' ${dir}/${filename}.temp2 > ${dir}/${filename}.temp; rm ${dir}/${filename}.temp2");

#	system("perl -pi -e 's/^NC_\\d+//g' ${dir}/${filename}.temp");
#	system("perl -pi -e 's/^\n\$//g' ${dir}/${filename}.temp");
#	system("perl -pi -e 's/\\.//g' ${dir}/${filename}.temp");	
}


sub MaskResidues_RemoveColumns (){
	
	## this function will mask residues in a given alignment and remove the aln columns with score below the defined threshold
	## a cutoff value should be defined from 0 (no cutoff) to 1
	## the score for each residue is found in the file *.MAFFT.Guidance2_res_pair_res.scr
	## the file contaning the column score values is *.MAFFT.Guidance2_res_pair_col.scr

	# opening temporary sequence file - created using the previous function MakeTempAlignment
	# this file does not have the sequence names, only the sequences
	open   (INFILE, "<${dir}/${filename}.temp") or die ("Cannot open temp.aligned file $!\n");
	
	my @align; # array to store each sequence (as a string)
	my @align2; # same as @align, but to be used for the new FASTA file (no columns removal nor residue masking, only low scored sequence exclusion
	my $seqnumber = 0; # variable to store total number of sequences in the temporary file
	my $align_length; # variable to store the total lenght of the alignment, i.e. total number of columns
	
	while(<INFILE>){		
		my($line) = $_;	
		chomp($line);
		my @seq = split(//,$line);
		my @seq2 = split(//,$line);
		$align_length = length($line);
		$align[$seqnumber] = \@seq;
		$align2[$seqnumber] = \@seq2;
		$seqnumber++;
	}


	##################################################
	##												##
	##			Part 1 - Masking Residues			##
	##		(using threshold from $maskcutoff)		##
	##												##
	##################################################

	# creating temporary res_score file
	system("cp ${dir}/${filename}.Guidance2_res_pair_res.scr ${dir}/${filename}.res_score"); 
	
	# removing header from temporary res_score file -> lines starting with "#"
	system("grep -v '#' ${dir}/${filename}.res_score > ${dir}/${filename}.res_score2"); 
	system("mv ${dir}/${filename}.res_score2 ${dir}/${filename}.res_score");

	# reformat ${filename}.res_score
	# depending on the dataset type [nuc|codon|aa] the res_score table has a different delimiter (guidance definitions!)
	system("perl -pi -e 's/^ +//g' ${dir}/${filename}.res_score"); # removing whitespace at the beginning of each line
	system("perl -pi -e 's/ +/\t/g' ${dir}/${filename}.res_score");# replacing whitespaces with "\t"

	# opening temporary res_score file
	open   (RESIDUESCORE, "<${dir}/${filename}.res_score") or die ("Cannot open ${dir}/${filename}.res_score $!\n");
	
	## mask residues
	## reading residues score file created by the script:
	## each column is separated by the symbol "!"
	while(<RESIDUESCORE>){
		my($line) = $_;
		chomp($line);
		my @row;
		@row = split ("\t",$line);

		my $aligncol     = $row[0];
		my $alignrow     = $row[1];
		my $residuescore = $row[2];
#		print $aligncol,"\t",$alignrow,"\t",$residuescore,"\n";
		
		# replace site with X or n, if residue score is below maskcutoff
		# it will not go through the positions with no score (e.g. gaps)
		if ($residuescore < $maskcutoff){ 
#			print $aligncol,"\t",$alignrow,"\t",$residuescore,"\n";
			if ($datatype eq "AA"){
#				print $align[$alignrow-1][$aligncol-1],"\n";
				$align[$alignrow-1][$aligncol-1] = "X";
#				print $align[$alignrow-1][$aligncol-1],"\n";
			}
			else {
				$align[$alignrow-1][$aligncol-1] = "n"
			}
		}
	}


	##################################################
	##												##
	##			Part 2 - Removing Columns			##
	##		(using threshold from $columncutoff)	##
	##												##
	##################################################
	
	# Remove unreliable columns below confidence score
	# the file contaning these score values is *.MAFFT.Guidance2_res_pair_col.scr
	system("cp ${dir}/${filename}.Guidance2_res_pair_col.scr ${dir}/${filename}.col_score");

	# reformatting ${filename}.col_score
	# removing whitespace at the beginning of each line
	system("perl -pi -e 's/^ +//g' ${dir}/${filename}.col_score"); 
	# replacing whitespaces with "\t"
	system("perl -pi -e 's/ +/\t/g' ${dir}/${filename}.col_score");
	
	# Opening reformated ${filename}.col_score
	open (COLUMNSCORE, "<${dir}/${filename}.col_score") or die ("Cannot open ${dir}/${filename}.col_score $!\n");
	
	## mask columns
	my $columns = " "; ## string with all present columns numbers separated by space

	while(<COLUMNSCORE>){
		my($line) = $_;

		chomp($line);
		next if $line =~ /^#/; # skipping header; starts with "#"
#		print $line,"\n";
		my @row = split ("\t",$line);

		my $col      = $row[0];
		my $colscore = $row[1];	
		$columns .= $col." ";
		
		# substitute all column with the symbol "!" if the score is below threshold
		if ($colscore < $columncutoff){			
#			print "Substitute column ", $col, " with score: ${colscore} \n";
			for my $seq (0..$seqnumber-1){
#				print $align[$seq][$col-1],"\n";
				$align[$seq][$col-1] = "!";
#				print $align[$seq][$col-1],"\n";
			}
		} 
	}


	## remove columns that are not present in the string $columns
	$InitialCol = $align_length;
	for my $colnumber (1..$align_length+1){
		my $substring = " ".$colnumber." ";
		if ($columns !~ /$substring/) {	
			#substitute column with "!"
			for my $seq (0..$seqnumber-1){
				$align[$seq][$colnumber-1] = "!";
			}			
		}
	}
	

	my @seqs;      # for PHYLIP
	my @seqs2; 	   # for FASTA
	my $sequence;  # for PHYLIP
	my $sequence2; # for FASTA
	
	for my $seq (0..$#align){
#		print $seq,"\n";
		$sequence = "";
		$sequence2 = "";
		for my $residue (0..$align_length - 1){
			my $aa = $align[$seq][$residue];
			my $aa2 = $align2[$seq][$residue]; # for FASTA
#			print $aa;			
			$sequence2 .= $aa2;

			next if $aa eq "!";
			$sequence .= $aa;
#			print $align[$seq][$residue];
		} 
#		print $sequence,"\n";
#		print $sequence2,"\n";
		push (@seqs,$sequence);
		push (@seqs2,$sequence2);
	}	

	

	## including sequence names in the aligment
	my $seqnames = `cat ${dir}/seq_names`;
#	print $seqnames,"\n";
	my @seqnames = split ("\n",$seqnames);

	## creating final files
	open   (OUTFASTA, ">${outdir}/${prefix}.${filename}.final.fasta");
	open   (OUTPHYLIP, ">${outdir}/${prefix}.${filename}.final.phy");
	
	my $finalseqnumber = $#seqnames + 1;
	my $finalseqlength = length($sequence);
	
	## for the PHYLIP outfile, the lenght from ^ to the sequence itself must be the same
	## let's define that that length must be the length of the longest sequence plus two withspaces
	my $seqNameLength = 0;
	for my $seqname (0..$#seqnames){
		my $sequence_name = $seqnames[$seqname];
		if (length($sequence_name) > $seqNameLength){
			$seqNameLength = length($sequence_name);
#			print $seqNameLength,"\n";
		}
	}
	
	
	print OUTPHYLIP $finalseqnumber," ",$finalseqlength,"\n"; 
	for my $seqname (0..$#seqnames){
		my $sequence_name = $seqnames[$seqname];
		my $whitespaces   = ($seqNameLength - length($sequence_name) + 2);
		$sequence_name    = $sequence_name." "x$whitespaces;
#		print $sequence_name,"!\n"; #ok

		my $sequence  = $seqs[$seqname]; # for phylip
		my $sequence2 = $seqs2[$seqname]; # for FASTA
		
		# final PHYLIP file (low scored sequences were not removed yet!, see next function)
		print OUTPHYLIP $sequence_name,$sequence,"\n";

		# final FASTA file (low scored sequences were not removed yet!, see next function)
		$sequence_name =~ s/ //g;
		print OUTFASTA ">",$sequence_name,"\n";
		
		$sequence2 =~ s/-//g;
		print OUTFASTA $sequence2,"\n";
	}

	# removing *temp files
	system("rm ${dir}/${filename}.temp ${dir}/seq_names ${dir}/${filename}.res_score ${dir}/${filename}.col_score");
}


sub RemoveLowScoredSeq (){
	## this function will remove sequences in a given alignment 
	## a cutoff value should be defined from 0 (no cutoff) to 1
	## the GUIDANCE score for the sequence is found in the file MSA.MAFFT.Guidance2_res_pair_seq.scr
	

	# Opening GUIDANCE score file MSA.MAFFT.Guidance2_res_pair_seq.scr
	open (SEQSCORE, "<${dir}/${filename}.Guidance2_res_pair_seq.scr");
#	print "${dir}/${filename}.Guidance2_res_pair_seq.scr\n";
	my @seq_to_remove;  # for Phylip
	my @seq_to_remove2; # for FASTA

	while(<SEQSCORE>){
		my($line) = $_;
		next if $line =~ /^#/; # skipping header; starts with "#"

		chomp($line);
		$line =~ s/^\s+//g;
		my @row;
		@row = split (" ",$line);
	
		my $seqnumber = $row[0];
		my $seqscore  = $row[1];
#		print $seqnumber,"-\t-",$seqscore,"\n";
		if ($seqscore < $seqcutoff){
#			print $seqnumber,"-\t-",$seqscore,"\n";
			push (@seq_to_remove,$seqnumber+1); # for PHYLIP
			push (@seq_to_remove2,$seqnumber*2-1); # for FASTA; removing the line with sequence name
			push (@seq_to_remove2,$seqnumber*2); # for FASTA; removing the line with the sequence itself
#			print "Sequence to remove: $seqnumber \n";
		}
	}
	

		# for PHYLIP
		my $lines_to_remove = join("d;",@seq_to_remove);
#		print $lines_to_remove,"\n";
		$lines_to_remove = $lines_to_remove."d";
		
		# for FASTA
		my $lines_to_remove2 = join("d;",@seq_to_remove2);
		$lines_to_remove2 = $lines_to_remove2."d";
#		print $#seq_to_remove,"\n";
 
		if ($#seq_to_remove >= 0){
			system ("sed '${lines_to_remove}' ${outdir}/${prefix}.${filename}.final.phy > ${outdir}/${prefix}.${filename}.notfinal.phy2");
			system("rm ${outdir}/${prefix}.${filename}.final.phy");
#			print ("sed '${lines_to_remove}' ${outdir}/${prefix}.${filename}.final.phy > ${outdir}/${prefix}.${filename}.final.phy2\n");
			system ("sed '${lines_to_remove2}' ${outdir}/${prefix}.${filename}.final.fasta > ${outdir}/${prefix}.${filename}.notfinal.fasta2");
			system("rm ${outdir}/${prefix}.${filename}.final.fasta");
#			print ("sed '${lines_to_remove2}' ${outdir}/${prefix}.${filename}.final.fasta > ${outdir}/${prefix}.${filename}.final.fasta2\n");
		}
		else{
			system ("cp ${outdir}/${prefix}.${filename}.final.phy ${outdir}/${prefix}.${filename}.final.phy2");
			system ("cp ${outdir}/${prefix}.${filename}.final.fasta ${outdir}/${prefix}.${filename}.final.fasta2");
		}


		my $phylip_line1 = `head -1 ${outdir}/${prefix}.${filename}.*inal.phy2`;
#		print $phylip_line1;
		chomp($phylip_line1);
		my @line = split (" ",$phylip_line1);
		my $NumberofSeq     = $line[0];
		my $AlnLength       = $line[1];
		my $FinalNumberSeqs = $NumberofSeq;
		
		$FinalCol = $AlnLength;

		# change first line of PHYLIP file if sequences are removed 
		# (total number of sequences should be smaller than that the original file)		
		if ($#seq_to_remove >= 0){
			$FinalNumberSeqs = $NumberofSeq - ($#seq_to_remove + 1);	
			system ("perl -pi -e s/\"$phylip_line1\"/'$FinalNumberSeqs $AlnLength'/ ${outdir}/${prefix}.${filename}.notfinal.phy2");
		}

		$InitialSeq = $NumberofSeq;
		$FinalSeq   = $FinalNumberSeqs;
}


sub Print_Info (){
	print "\n\n##########################################################\n"; 
	print "########### Running script 03-guidance_v1.2c.pl ##########\n";
	print "##########################################################\n\n";
	print my $date_time = `date`;
	print "\nFILE = ${dir}/${filename}.aln.With_Names\n\n";
	print "DIR=${dir}\n";
	print "FILENAME=${filename}\n\n";
	print "TYPE\tCUTOFF\tBEFORE\tFINAL\tPERC_KEPT(%)\n";

	# sequences
	print "SEQ\t$seqcutoff\t$InitialSeq\t$FinalSeq\t";
	printf "%5.2f\n",$FinalSeq/$InitialSeq*100;  # MACR -- I modified all printf sentences (4). There was an error, more arguments than needed, ",'\n'" at the end

	# columns	
	print "COL\t$columncutoff\t$InitialCol\t$FinalCol\t";
	printf "%5.2f\n",$FinalCol/$InitialCol*100;

#	$InitialRes = `sed '1d' ${outdir}/${prefix}.${filename}.*.phy2 | awk '{print \$2}'| tr -dc 'A-Z'| wc -c | perl -pi -e 's/\\s//g'`;
	$InitialRes = `sed '1d' ${outdir}/${prefix}.${filename}.*.phy2 | awk '{print \$2}'| tr -dc 'A-Z'| wc -c | perl -p -e 's/\\s//g'`;	# MACR -- corrected here, "i" requires input file
	
	chomp($InitialRes);
	$InitialRes =~ s/\s//g;

	# final res means unmasked ones
#	my $MaskedRes = `sed '1d' ${outdir}/${prefix}.${filename}.*.phy2 | awk '{print \$2}'| tr -dc 'X'| wc -c | perl -pi -e 's/\\s//g'`;	# MACR -- corrected here, "i" requires input file
	my $MaskedRes = `sed '1d' ${outdir}/${prefix}.${filename}.*.phy2 | awk '{print \$2}'| tr -dc 'X'| wc -c | perl -p -e 's/\\s//g'`;
	chomp($MaskedRes);
	$MaskedRes =~ s/\s//g;	
	$FinalRes = $InitialRes - $MaskedRes;
	
	print "RES\t$maskcutoff\t$InitialRes\t$FinalRes\t";
	
	if ($InitialRes > 0){
		printf "%5.2f\n",$FinalRes/$InitialRes*100;
	}
	else{
		printf "%5.2f\n","0";
	}
	
	if ($maskcutoff == 0){
		print "\n#If Residue cutoff is 0, but if in the final aln there are some masked residues,\nit means that the original FASTA file already contained these masked residues (\"X\")\n";
	}

	# if $FinalSeq == 0, then no sequences were retained -> no final alignment should be kept
	# if $FinalSeq <= 3, there is no point to keep the alignment either for further phylogenetics analysis
	if ($FinalSeq <= 3){
		#system("rm ${outdir}/${prefix}.${filename}.*final.*");
		print "\nEXCLUDED ALN with 3 or less sequences:\n${dir}${filename}.aln.With_Names\n";
	}
	else{
		if ($FinalSeq == $InitialSeq){
			system("mv ${outdir}/${prefix}.${filename}.final.phy2 ${outdir}/${prefix}.${filename}.final.phy");
			system("mv ${outdir}/${prefix}.${filename}.final.fasta2 ${outdir}/${prefix}.${filename}.final.fasta");
		}
		else{
			system("mv ${outdir}/${prefix}.${filename}.notfinal.phy2 ${outdir}/${prefix}.${filename}.notfinal.phy");
			system("mv ${outdir}/${prefix}.${filename}.notfinal.fasta2 ${outdir}/${prefix}.${filename}.notfinal.fasta");			
		}
	}
	print "---\n";
}

MakeTempAlignment ();
MaskResidues_RemoveColumns ();
RemoveLowScoredSeq ();
Print_Info ();

exit;
