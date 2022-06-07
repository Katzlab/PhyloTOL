
echo Node being used is ${HOSTNAME//+([[:alpha:]-.])}
unset SGE_ROOT

mypwd=$(echo `pwd`)
echo $mypwd

## version 2.2 

## command example: ./exeGuidance_v2.sh -i /home/mmfonseca/2015-Katz/2015-September/00-FASTA/OG5_129999.fas -o /state/partition1/mmfonseca/2015-Katz/ -t 2 -c /home/mmfonseca/2015-Katz/2015-September/01-Guidance_out -s 0.4 -l 0.5 -r 0.4 -m nr

# PARSING THE  ARGUMENTS
usage() { echo "Usage: $0 [-i FASTA] [-o SCRATCHDIR] [-t NTHREADS] [-c OUTHOMEDIR] [-g GUIDITER] [-s SEQCUT] [-l COLCUT] [-r RESCUT] [-m MODE]" ; exit 1; }

while getopts ":h:i:o:t:c:g:s:l:r:m:" opt; do
        case $opt in
                i)
                        FASTA=$OPTARG
                        ;;
                o)
                        SCRATCHDIR=$OPTARG
                        ;;
                t)
                        NTHREADS=$OPTARG
                        ;;
                c)
                        OUTHOMEDIR=$OPTARG
                        ;;
                g)
                        GUIDITER=$OPTARG
                        ;;
                s)
                        SEQCUT=$OPTARG
                        ;;
                l)
                        COLCUT=$OPTARG
                        ;;

                r)
                        RESCUT=$OPTARG
                        ;;

                m)
                        MODE=$OPTARG
                        ;;                        
                        
                h)
                        usage
                        ;;
               \?)
                        echo "Invalid option: -$OPTARG"
                        usage
                        exit 1
                        ;;
                :)
                        if [ "$OPTARG" != "h" ]; then
                        	echo "Option -$OPTARG requires an argument."
                        fi
                        echo
                        echo "------------------------------------------------------------------------------"
                        echo "  exeGuidance_v2.2.sh (06/06/2016)"
                        echo "      for Tree of Eukaryotes project"
                        echo
                        echo "      script developed by Miguel M Fonseca"
                        echo "      mig.m.fonseca@gmail.com"
                        echo "------------------------------------------------------------------------------"
                        echo
                        echo "To execute "exeGuidance_v2.2.sh" you need the following arguments:"
                        echo
                        echo "-i"
                        echo "    absolute path to the unaligned fasta file."
                        echo
                        echo "-o"
                        echo "    absolute path to temporary directory where guidance run will be executed and files will be stored."
                        echo "    Inside this directory, the script will make a directory with the name of the gene"
                        echo "    (it might change a bit in your case, I might need to adapt the script)"
                        echo
                        echo "-c"
                        echo "    absolute path to the output directory of the guidance execution."
                        echo "    Inside this directory, a new directory will be created, specifically for the fasta file that you are running."
                        echo "    (again we might need to adapt it for your case)"
                        echo 
                        echo "-t"
                        echo "    Number of threads to used in the guidance analysis."
                        echo "    This option is not very efficient for guidance, hence just use \"-t 1\" or \"-t 2\""
                        echo
                        echo "-g"
                        echo "    maximum number of Guidance iterations. Choose a number or n"
                        echo
                        echo "-s"
                        echo "    sequence cutoff. Choose a value between 0 and 1"
                        echo
                        echo "-l"
                        echo "    columns cutoff. Choose a value between 0 and 1"
                        echo
                        echo "-r"
                        echo "    residues cutoff. Choose a value between 0 and 1"
                        echo
                        echo "-m"
                        echo "    mode. If 'nr', the RAxML is not executed. Any other value will be the opposite
                        echo
                        echo " # NOTE: to change the threshold values: column, residue or sequence cutoff you need to change lines 166, 167 and 235 of the script."
                        usage
                        exit 1
                        ;;
          esac
done
shift $(( OPTIND - 1 ))

if [ -z "${FASTA}" ] || [ -z "${SCRATCHDIR}" ] || [ -z "${NTHREADS}" ] || [ -z "${OUTHOMEDIR}" ] || [ -z "${GUIDITER}" ] || [ -z "${SEQCUT}" ] || [ -z "${COLCUT}" ] || [ -z "${RESCUT}" ] || [ -z "${MODE}" ]; then usage; fi


	##########################
	## Loading Modules
	##########################
		
#		module load openmpi/1.8.4
#		module load openmpi/1.8.5
#		module load perl/5.20.1 
#		module load bio/mafft/7.212
#		module load bio/guidance/1.5

	###########################

	
	##########################
	## Getting Gene family name
	##########################
	
		declare -a myarr=(`echo "$FASTA" | sed 's/\// /g'`)
		lastpos=$(echo ${#myarr[*]})
		genename=$(echo ${myarr[$lastpos-1]})
	
	##########################
	
	
	
	####################################################
	##												  ## 
	##				EXECUTING GUIDANCE 1.5			  ##
	##												  ##
	####################################################	

		run=1
		exe=1

	while [ "$exe" -eq 1 ]
	do
		echo "[$run] Guidance iteration" 
		touch $SCRATCHDIR/RUN$run    # MACR
		date > $SCRATCHDIR/RUN$run   # MACR
		
		## for the first iteration
		## cp FASTA to scratch directory
		if [ "$run" == 1 ]; then
			## create scratch folder
			mkdir -p $SCRATCHDIR/$genename
			## cp FASTA to scratch directory
			cp $FASTA $SCRATCHDIR/$genename
		fi	
		
		# MACR ---- In the two lines commented MACR above and this one I corrected a bug. It was creating a file RUN1 before creating its containing folder.
		mv $SCRATCHDIR/RUN$run $SCRATCHDIR/$genename/RUN$run 
														  
		
		## Execute Guidance 2.0? actually version 2.0 is very slow, we will use version 1.5
		## MACR -- We are using guidance 2.02

			##########################
			## calculate number of  
			## sequences in FASTA	
			##########################
			
			## identify if FASTA file has less or more than 200 sequences
			## if less, then gensi algorithm will be used with MAFFT
			## if more, then auto option will be used

				nlines=$(grep '^>' $SCRATCHDIR/$genename/$genename | wc -l  | awk '{print $1}')
				if [ "$nlines" -le 200 ]; then
					mafftalg="genafpair"
				else
					mafftalg="auto"
				fi
	
				echo "[$run] number of sequences in $SCRATCHDIR/$genename/$genename = $nlines"
				echo "[$run] mafftalg = $mafftalg"	
			
			
			if [ "$nlines" -ge 4 ]; then
			
				##########################
				## execute guidance 2.02
				##########################
			
				scratchFASTA="$SCRATCHDIR/$genename/$genename"
				scratchDIR="$SCRATCHDIR/$genename/"
				guidancepl="guidance.v2.02/www/Guidance/guidance.pl" # new guidance 2.02
				seqcutoff="$SEQCUT"
				colcutoff="$COLCUT"
				rescutoff="$RESCUT"
				
				echo ""
				echo "[$run] lauching guidance.pl ..."
				#echo "perl $guidancepl --seqFile $scratchFASTA --msaProgram MAFFT --seqType aa --outDir $scratchDIR --seqCutoff $seqcutoff --colCutoff $colcutoff --proc_num $NTHREADS --outOrder as_input --bootstraps 100 --MSA_Param \\\-'-${mafftalg} --maxiterate 1000 --bl 62 --anysymbol --thread 1'" | sh
				echo "perl $guidancepl --seqFile $scratchFASTA --msaProgram MAFFT --seqType aa --outDir $scratchDIR --seqCutoff $seqcutoff --colCutoff $colcutoff --outOrder as_input --bootstraps 10 --MSA_Param \\\-'-${mafftalg} --maxiterate 1000 --thread $NTHREADS --bl 62 --anysymbol'" | sh
				echo "[$run] ... closing guidance.pl"
				echo ""
				
				## Here I create the variable $done, which will help to specify number of guidance iterations tu run
				if [ "$GUIDITER" = "n" ]; then
					DONE=0
				else
					DONE=$GUIDITER		
				fi
				
				## get number of sequences below seq threshold
				NLOWSEQ=$(grep -v '#' $SCRATCHDIR/$genename/MSA.MAFFT.Guidance2_res_pair_seq.scr | awk '$2 < "'$seqcutoff'" {print}' | wc -l)
				echo "[$run] Number of sequences below threshold: $NLOWSEQ"
				## get number of sequences above seq threshold
				NABOVESEQ=$(grep -v '#' $SCRATCHDIR/$genename/MSA.MAFFT.Guidance2_res_pair_seq.scr | awk '$2 >= "'$seqcutoff'" {print}' | wc -l)
				echo "[$run] Number of sequences above threshold: $NABOVESEQ"
				echo ""
				
				## if number of final sequences is below 4 -> finish guidance
				if [ "$NABOVESEQ" -le 3 ]; then
					exe=0
					echo "[$run] Final number of sequences is below 4 -> Quitting Guidance"
					touch $SCRATCHDIR/$genename/NAboveCutoffBelow4
					touch $SCRATCHDIR/$genename/NSeqBelow4
					
					#mv $SCRATCHDIR/$genename/Seqs.Orig.fas.FIXED.Without_low_SP_Seq $SCRATCHDIR/$genename/$genename.less3
				fi
				
				## if there are no sequences to be removed -> finish guidance
				if [ "$NLOWSEQ" -eq 0 ]; then
					exe=0
					touch $SCRATCHDIR/$genename/NLowCutoffEquals0
					echo "[$run] There are no sequences to be removed -> Quitting Guidance"
					
					if [ "$NABOVESEQ" -ge 4 ]; then
						echo "[$run] There are no sequences to be removed and number of sequences above threshold >= 4 -> Quitting Guidance"
					fi
				fi
				
				## if maximum of guidance iterations -> finish guidance MACR
				if [ "$DONE" -eq "$run" ]; then
					exe=0
					echo "[$run] maximum guidance iterations: $DONE <-- Finished"
					
					if [ "$NABOVESEQ" -ge 4 ]; then
						echo "[$run] $DONE iterations already done and number of sequences above threshold >= 4 -> Quitting Guidance"
					fi
				fi
	

				if [ "$NLOWSEQ" -ge 1 ]; then
					if [ "$DONE" != "$run" ]; then
						echo "[$run] There are sequences to be removed -> Go to next Guidance iteration"
						echo "[$run] copying file Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names as the new input file:"
						echo "cp $SCRATCHDIR/$genename/Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names $SCRATCHDIR/$genename/$genename"
						cp $SCRATCHDIR/$genename/Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names $SCRATCHDIR/$genename/$genename
						echo "[$run] Removing Guidance output files (except log) that will not be needed for the next iteration"
						rm $SCRATCHDIR/$genename/Seqs* $SCRATCHDIR/$genename/MSA.MAFFT* $SCRATCHDIR/$genename/ENDS_OK
						echo "[$run] rm $SCRATCHDIR/$genename/Seqs* $SCRATCHDIR/$genename/MSA.MAFFT* $SCRATCHDIR/$genename/ENDS_OK"
						echo "[$run] Actualizing iteration number..."
						run=$(( run+1 ))
						echo "[$run] ...Iteration number updated"
						echo "---"	
					fi				
				fi
			else
				echo "[$run] GUIDANCE ABORTED [2]: the initial sequence file contained 3 or less sequences in total"
				echo "[$run] GUIDANCE ABORTED [3]: number of sequences in $SCRATCHDIR/$genename/$genename = $nlines"
				exe=0
				touch $SCRATCHDIR/$genename/NInitialSeqBelow4
				touch $SCRATCHDIR/$genename/NSeqBelow4
			fi
	done
	
	##########################
	## execute 
	## 03-guidance_v1.2c.pl
	##########################
	
#	if [ ! -f $SCRATCHDIR/$genename/NInitialSeqBelow4 ] | [ ! -f $SCRATCHDIR/$genename/NAboveCutoffBelow4 ]; then  # MACR modified to continue pipeline after guidance crash 
	if [ ! -f $SCRATCHDIR/$genename/NSeqBelow4 ]; then  # MACR modified to continue pipeline after guidance crash 
#	if [ -f $]
				
		## After Guidance execution let's see if there were sequences removed
		## to do it, we have to execute the script 03-guidance_v1.2c.pl
		echo ""
		echo "[$run] Launching perl script 03-guidance_v1.2c.pl..."
	
		myscript="03-guidance_v1.2c.pl"
		Guidance_out_DIR="$scratchDIR"
	#	SITE_CUTOFF=0.4  # MAC - I replaced this by a rescutoff
		COL_CUTOFF="$colcutoff"
		SEQ_CUTOFF="$seqcutoff"
		RES_CUTOFF="$rescutoff"
		DATATYPE="AA"
		PREFIX="${genename}.run${run}"

#		echo "perl $myscript -filename MSA.MAFFT.aln.With_Names -inDir $Guidance_out_DIR -outDir $Guidance_out_DIR -siteCutoff $RES_CUTOFF -colCutoff $COL_CUTOFF -seqCutoff $SEQ_CUTOFF -dataType $DATATYPE -prefix $PREFIX" | sh > $Guidance_out_DIR/${PREFIX}.log
#		echo "perl $myscript --filename MSA.MAFFT --inDir $Guidance_out_DIR --outDir $Guidance_out_DIR --siteCutoff $RES_CUTOFF --colCutoff $COL_CUTOFF --seqCutoff $SEQ_CUTOFF --dataType $DATATYPE --prefix $PREFIX" #| sh > $Guidance_out_DIR/${PREFIX}.log
		perl $myscript --filename MSA.MAFFT --inDir $Guidance_out_DIR --outDir $Guidance_out_DIR --siteCutoff $RES_CUTOFF --colCutoff $COL_CUTOFF --seqCutoff $SEQ_CUTOFF --dataType $DATATYPE --prefix $PREFIX   # MACR modified here, it wasn't working properly
		echo "[$run] ...perl script 03-guidance_v1.2c.pl finished"

	
		echo "[$run] Removing temporary files MSA.MAFFT* Seqs* *fasta"
		rm $SCRATCHDIR/$genename/MSA.MAFFT*
		rm $SCRATCHDIR/$genename/Seqs* 
		rm $SCRATCHDIR/$genename/*fasta
		echo ""

		##########################
		## saving Post-guidance ##    MACR
		##########################

		echo "renaming post-guidance alignment" 
		cp $SCRATCHDIR/$genename/*phy $SCRATCHDIR/$genename/$genename.Post
	
		##########################
		## execute trimal
		##########################
		#module load bio/trimal/1.3
		# MACR -- Changes here to read trimAl from local folder
		
		echo "[$run] Executing trimal to remove columns with 95% or more gaps..."
		echo "[$run] `head -1 $SCRATCHDIR/$genename/*phy` Head -1 of alignment BEFORE filtering gaps"
#		trimal -in $SCRATCHDIR/$genename/*phy -out $SCRATCHDIR/$genename/$genename.95gapTrimmed.nex -gapthreshold 0.05 -nexus
#		trimal -in $SCRATCHDIR/$genename/*phy -out $SCRATCHDIR/$genename/$genename.95gapTrimmed.fas -gapthreshold 0.05 -fasta
#		trimal -in $SCRATCHDIR/$genename/*phy -out $SCRATCHDIR/$genename/$genename.95gapTrimmed.phy -gapthreshold 0.05 -phylip
		trimal-trimAl/source/trimal -in $SCRATCHDIR/$genename/*phy -out $SCRATCHDIR/$genename/$genename.95gapTrimmed.nex -gapthreshold 0.05 -nexus
		trimal-trimAl/source/trimal -in $SCRATCHDIR/$genename/*phy -out $SCRATCHDIR/$genename/$genename.95gapTrimmed.fas -gapthreshold 0.05 -fasta
		trimal-trimAl/source/trimal -in $SCRATCHDIR/$genename/*phy -out $SCRATCHDIR/$genename/$genename.95gapTrimmed.phy -gapthreshold 0.05 -phylip
		echo "[$run] `head -1 $SCRATCHDIR/$genename/$genename.95gapTrimmed.phy` Head -1 of alignment AFTER filtering gaps"
		echo "[$run] ...trimal finished"
		echo ""


		if [ "$MODE" != "nr" ]; then
		
			##########################
			## execute RAxML
			##########################

			## Now, that we have the final alignment, 
			## we can inferr the respective gene tree
			## using RAxML

			if [ "$NABOVESEQ" -ge 4 ]; then			
				echo "[$run] Launching RAxML..."
			#	module load bio/raxml/8.1.7
				
				## Define RAxML input file
					cp $SCRATCHDIR/$genename/$genename.95gapTrimmed.phy $SCRATCHDIR/$genename/$genename.RAxML.phy
					PHYLIP="$SCRATCHDIR/$genename/$genename.RAxML.phy"
		
				## RAxML commands
#				raxmlHPC -s $PHYLIP -m PROTGAMMALG -f d -p 12345 -# 1 -n $genename # -T 3
				raxmlHPC-PTHREADS-AVX2 -s $PHYLIP -m PROTGAMMALG -f d -p 12345 -# 1 -n $genename -T $NTHREADS
				echo "[$run] ...RAxML finished"
				echo ""
			else
				echo "[$run] RAxML ABORTED: RAxML did not run because the number of sequences in final file is below 4"
				echo ""
			fi
		fi
	fi   ## MACR -- changes here in order to continue even if guidance fails with numb of seqs < 4. I took next lines out of the loop

	## move output to home directory
	echo "[$run] moving scratch output to home output directory"
	echo "[$run] mkdir $OUTHOMEDIR/$genename.output/"
#	echo "[$run] cp $scratchDIR/* $OUTHOMEDIR/$genename.output/"
	echo "[$run] cp $SCRATCHDIR/$genename/* $OUTHOMEDIR/$genename.output/"
	echo "[$run] rm -r $scratchDIR/"
	mkdir $OUTHOMEDIR/$genename.output/
#	cp $scratchDIR/* $OUTHOMEDIR/$genename.output/
	cp $SCRATCHDIR/$genename/* $OUTHOMEDIR/$genename.output/
#	rm -r $scratchDIR/
	rm -r $SCRATCHDIR/$genename/
	
	if [ "$MODE" != "nr" ]; then
		## move other RAxML outputs that are saved in the directory where we lauched the bash file
		echo "[$run] home RAxML files to home output directory"
		echo "[$run] mv $mypwd/RAxML_*$genename* $OUTHOMEDIR/$genename.output/"
		mv $mypwd/RAxML_*$genename* $OUTHOMEDIR/$genename.output/ 
	fi
#	fi
