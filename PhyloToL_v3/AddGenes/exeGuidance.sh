echo Node being used is ${HOSTNAME//+([[:alpha:]-.])}
unset SGE_ROOT

mypwd=$(echo `pwd`)
echo $mypwd

## version 1 --- Adapted from exeGuidance v2.2 witten by Miguel -- MACR

## command example: ./exeGuidance_v2.sh -i /home/mmfonseca/2015-Katz/2015-September/00-FASTA/OG5_129999.fas -o /state/partition1/mmfonseca/2015-Katz/ -t 2 -c /home/mmfonseca/2015-Katz/2015-September/01-Guidance_out -s 0.4 -l 0.5 -r 0.4 -m nr

# PARSING THE  ARGUMENTS
usage() { echo "Usage: $0 [-i FASTA] [-o SCRATCHDIR] [-t NTHREADS] [-c OUTHOMEDIR] [-s SEQCUT]" ; exit 1; }

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

		s)
				SEQCUT=$OPTARG
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
				echo "  exeGuidance.sh"
				echo 
				echo "      Mario A. Cerón-Romero"
				echo "      mig.m.fonseca@gmail.com"
				echo "------------------------------------------------------------------------------"
				echo
				echo "To execute "exeGuidance.sh" you need the following arguments:"
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
				echo "-s"
				echo "    sequence cutoff. Choose a value between 0 and 1"
				echo
				usage
				exit 1
				;;
	esac
done
shift $(( OPTIND - 1 ))

if [ -z "${FASTA}" ] || [ -z "${SCRATCHDIR}" ] || [ -z "${NTHREADS}" ] || [ -z "${OUTHOMEDIR}" ] || [ -z "${SEQCUT}" ] ; then usage; fi

	

# Getting Gene family name
declare -a myarr=(`echo "$FASTA" | sed 's/\// /g'`)
lastpos=$(echo ${#myarr[*]})
genename=$(echo ${myarr[$lastpos-1]})

		
## create scratch folder
mkdir -p $SCRATCHDIR/$genename
## cp FASTA to scratch directory
cp $FASTA $SCRATCHDIR/$genename

## **** Execute Guidance 2.02 *****

## identify if FASTA file has less or more than 200 sequences
## if less, then gensi algorithm will be used with MAFFT
## if more, then auto option will be used

nlines=$(grep '^>' $SCRATCHDIR/$genename/$genename | wc -l  | awk '{print $1}')
if [ "$nlines" -le 200 ]; then
	mafftalg="genafpair"
else
	mafftalg="auto"
fi

echo "number of sequences in $SCRATCHDIR/$genename/$genename = $nlines"
echo "mafftalg = $mafftalg"	
			
			
if [ "$nlines" -ge 4 ]; then
	outD="$SCRATCHDIR/OGs_filteredNoName/"
	scratchFASTA="$SCRATCHDIR/$genename/$genename"
	scratchDIR="$SCRATCHDIR/$genename/"
#	guidancepl="./guidance.v2.02/www/Guidance/guidance.pl" 
	guidancepl="~/PhyloTOL_paral1/Scripts/guidance.v2.02/www/Guidance/guidance.pl"
	seqcutoff="$SEQCUT"
	colcutoff="0"
	
	echo ""
	echo "lauching guidance.pl ..."
	echo "perl $guidancepl --seqFile $scratchFASTA --msaProgram MAFFT --seqType aa --outDir $scratchDIR --seqCutoff $seqcutoff --colCutoff $colcutoff --outOrder as_input --bootstraps 10 --MSA_Param \\\-'-${mafftalg} --maxiterate 1000 --thread $NTHREADS --bl 62 --anysymbol'" | sh
	echo "... closing guidance.pl"
	echo ""
	
	## get number of sequences below seq threshold
	NLOWSEQ=$(grep -v '#' $SCRATCHDIR/$genename/MSA.MAFFT.Guidance2_res_pair_seq.scr | awk '$2 < "'$seqcutoff'" {print}' | wc -l)
	echo "Number of sequences below threshold: $NLOWSEQ"
	## get number of sequences above seq threshold
	NABOVESEQ=$(grep -v '#' $SCRATCHDIR/$genename/MSA.MAFFT.Guidance2_res_pair_seq.scr | awk '$2 >= "'$seqcutoff'" {print}' | wc -l)
	echo "Number of sequences above threshold: $NABOVESEQ"
	echo ""
	
	
	## if there are no sequences to be removed -> finish guidance
	if [ "$NLOWSEQ" -eq 0 ]; then
		echo "There are no sequences to be removed -> Original file will be kept"
		cp $FASTA $outD/$genename
		
	fi	

	if [ "$NLOWSEQ" -ge 1 ]; then
			echo "There are sequences to be removed"
			echo "Saving file Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names"
			echo "cp $SCRATCHDIR/$genename/Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names $SCRATCHDIR/$genename/$genename"
			cp $SCRATCHDIR/$genename/Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names $outD/$genename
	fi
else
	echo "GUIDANCE ABORTED [2]: the initial sequence file contained 3 or less sequences in total"
	echo "GUIDANCE ABORTED [3]: number of sequences in $SCRATCHDIR/$genename/$genename = $nlines"
fi

echo "Removing temporal Guidance files"
rm -r $SCRATCHDIR/$genename
echo "rm -r $SCRATCHDIR/$genename"

