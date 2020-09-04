#!/usr/bin/env python3

##__Updated__: 2020-08-24
##__Author__: Xyrus Maurer-Alcala; maurerax@gmail.com; xyrus.maurer-alcala@izb.unibe.ch
##__Usage__: python 5_GCodeTranslate.py --help


##########################################################################################
## This script is intended to aid in identifying the genetic code of the data given		##
##                                                                                      ##
## Prior to running this script, ensure the following:									##
##                                                                                      ##
## 1. You have assembled your transcriptome and COPIED the 'assembly' file              ##
##    (contigs.fasta, or scaffolds.fasta) to the PostAssembly Folder                    ##
## 2. Removed small sequences (usually sequences < 300bp) with 1_ContigFiltStats.py     ##
## 3. Removed SSU/LSU sequences from your Fasta File                                    ##
## 4. Classified your sequences as Strongly Prokaryotic/Eukaryotic or Undetermined      ##
## 5. Classified the Non-Strongly Prokaryotic sequences into OGs                        ##
## 6. You either know (or have inferred) the genetic code of the organism               ##
##                                                                                      ##
##             E-mail Xyrus (author) for help if needed: maurerax@gmail.com             ##
##                                                                                      ##
##                              Next Script(s) to Run:                                  ##
##                 6_FilterPatials.py (in FinalizeTranscripts Folder)                   ##
##                                                                                      ##
##########################################################################################

import argparse, os, re, sys
from argparse import RawTextHelpFormatter,SUPPRESS

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data.CodonTable import CodonTable


#-------------------------- Set-up Codon Tables (Genetic Codes) --------------------------#

blepharisma_table = CodonTable(forward_table={
	'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
	'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
	'TAT': 'Y', 'TAC': 'Y',
	'TGT': 'C', 'TGC': 'C', 'TGA': 'W', 'TGG': 'W',
	'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
	'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
	'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
	'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
	'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
	'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
	'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
	'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
	'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
	'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
	'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
	'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'},
	start_codons = [ 'ATG'],
	stop_codons = ['TAA','TAG'])

condylostoma_table = CodonTable(forward_table={
	'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
	'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
	'TAT': 'Y', 'TAC': 'Y', 'TAA': 'Q', 'TAG': 'Q',
	'TGT': 'C', 'TGC': 'C', 'TGA': 'W', 'TGG': 'W',
	'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
	'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
	'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
	'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
	'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
	'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
	'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
	'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
	'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
	'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
	'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
	'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'},
	start_codons = [ 'ATG'],
	stop_codons = [''])

c_uncinata_table = CodonTable(forward_table={
	'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
	'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
	'TAT': 'Y', 'TAC': 'Y',             'TAG': 'Q',
	'TGT': 'C', 'TGC': 'C', 'TGA': 'Q', 'TGG': 'W',
	'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
	'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
	'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
	'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
	'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
	'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
	'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
	'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
	'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
	'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
	'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
	'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'},
	start_codons = [ 'ATG'],
	stop_codons = ['TAA'])

euplotes_table = CodonTable(forward_table={
	'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
	'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
	'TAT': 'Y', 'TAC': 'Y',
	'TGT': 'C', 'TGC': 'C', 'TGA': 'C', 'TGG': 'W',
	'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
	'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
	'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
	'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
	'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
	'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
	'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
	'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
	'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
	'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
	'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
	'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'},
	start_codons = [ 'ATG'],
	stop_codons = ['TAA','TAG'])

myrionecta_table = CodonTable(forward_table={
	'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
	'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
	'TAT': 'Y', 'TAC': 'Y', 'TAA': 'Y', 'TAG': 'Y',
	'TGT': 'C', 'TGC': 'C',             'TGG': 'W',
	'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
	'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
	'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
	'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
	'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
	'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
	'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
	'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
	'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
	'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
	'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
	'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'},
	start_codons = [ 'ATG'],
	stop_codons = ['TGA'])

no_stop_table = CodonTable(forward_table={
	'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
	'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
	'TAT': 'Y', 'TAC': 'Y', 'TAA': 'X', 'TAG': 'X',
	'TGT': 'C', 'TGC': 'C', 'TGA': 'X', 'TGG': 'W',
	'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
	'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
	'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
	'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
	'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
	'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
	'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
	'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
	'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
	'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
	'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
	'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'},
	start_codons = [ 'ATG'],
	stop_codons = [''])

peritrich_table = CodonTable(forward_table={
	'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
	'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
	'TAT': 'Y', 'TAC': 'Y', 'TAA': 'E', 'TAG': 'E',
	'TGT': 'C', 'TGC': 'C',             'TGG': 'W',
	'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
	'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
	'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
	'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
	'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
	'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
	'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
	'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
	'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
	'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
	'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
	'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'},
	start_codons = [ 'ATG'],
	stop_codons = ['TGA'])

tag_table = CodonTable(forward_table={
	'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
	'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
	'TAT': 'Y', 'TAC': 'Y', 'TAA': 'Q',
	'TGT': 'C', 'TGC': 'C', 'TGA': 'Q', 'TGG': 'W',
	'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
	'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
	'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
	'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
	'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
	'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
	'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
	'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
	'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
	'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
	'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
	'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'},
	start_codons = [ 'ATG'],
	stop_codons = ['TAG'])


#------------------------------ Colors For Print Statements ------------------------------#
class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   ORANGE = '\033[38;5;214m'
   PURPLE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'


#------------------------------- Main Functions of Script --------------------------------#

###########################################################################################
###------------------------- Checks the Command Line Arguments -------------------------###
###########################################################################################

def check_args():

	parser = argparse.ArgumentParser(description=
	color.BOLD + '\n\nThis script will '+color.RED+'Translate '+color.END+color.BOLD+'a '\
	'given Fasta file of CDS\nsequences using a given'+color.PURPLE+' Genetic Code.'+color.END+\
	color.BOLD+usage_msg(), usage=SUPPRESS, formatter_class=RawTextHelpFormatter)


	required_arg_group = parser.add_argument_group(color.ORANGE+color.BOLD+'Required Options'+color.END)

	required_arg_group.add_argument('--input_file','-in', action='store',
	help=color.BOLD+color.GREEN+' Fasta file with CDSs\n'+color.END)

	optional_arg_group = parser.add_argument_group(color.ORANGE+color.BOLD+'Options'+color.END)

	optional_arg_group.add_argument('--genetic_code','-g', action='store', default='universal',
	help=color.BOLD+color.GREEN+' Genetic code to use for translation\n (default = '\
	'"universal")\n'+color.END)

	optional_arg_group.add_argument('--no_RP','-no_RP', action='store_true',
	help=color.BOLD+color.GREEN+' Allows files to "skip" the removal\n of Partial Transcripts\n'\
	+color.END)

	optional_arg_group.add_argument('--list_codes','-codes', action='store_true',
	help=color.BOLD+color.GREEN+' Lists supported genetic codes\n'+color.END)

	optional_arg_group.add_argument('-author', action='store_true',
	help=color.BOLD+color.GREEN+' Prints author contact information\n'+color.END)


	if len(sys.argv[1:]) == 0:
		print (parser.description)
		print ('\n')
		sys.exit()

	args = parser.parse_args()

	quit_eval = return_more_info(args)
	if quit_eval > 0:
		sys.exit()

	args.tsv_file = args.input_file.replace('.fasta','_allOGCleanresults.tsv')

	### Adding in names to 'arg' class for more easy use throughout the script
	args.ntd_out = args.input_file.split('.fas')[0]+'_'+args.genetic_code.title()+'_NTD.ORF.fasta'
	args.aa_out = args.input_file.split('.fas')[0]+'_'+args.genetic_code.title()+'_AA.ORF.fasta'
	args.tsv_out = args.input_file.split('.fas')[0]+'_'+args.genetic_code.title()+'_allOGCleanresults.tsv'

	args.home_folder = '../'+args.input_file.split('/')[1]
	args.Usearch_Folder = '../'+args.input_file.split('/')[1]+'/UsearchOG'
	args.StopFreq = '../'+args.input_file.split('/')[1]+'/StopCodonFreq'
	args.ForPartials = '../FinalizeTranscripts/'

	return args


###########################################################################################
###------------------------------- Script Usage Message --------------------------------###
###########################################################################################

def usage_msg():
	return (color.BOLD+color.RED+'\n\nExample usage:'+color.CYAN+' python 3g_GCodeTranslate.py'\
	' --input_file ../Stentor_coeruleus.WGS.CDS.Prep/Stentor_coeruleus.WGS.CDS.Renamed.fasta'\
	' --genetic_code Universal'+color.END)


##########################################################################################
###-------- Storage for LARGE (Annoying) Print Statements for Flagged Options ---------###
##########################################################################################

def return_more_info(args):

	valid_arg = 0

	supported_gcodes_names = ['bleph','blepharisma','chilo','chilodonella','condy',\
	'condylostoma','none','eup','euplotes','peritrich','vorticella','ciliate','universal',\
	'taa','tag','tga', 'mesodinium', 'myrionecta']

	supported_gcodes_list = ['Blepharisma\t(TGA = W)','Chilodonella\t(TAG/TGA = Q)','Ciliate\t\t(TAR = Q)',\
	'Conylostoma\t(TAR = Q, TGA = W)','Euplotes\t(TGA = C)','Peritrich\t(TAR = E)','Mesodinium\t(TAR = Y)',\
	'Myrionecta\t(TAR = Y)','None\t\t(TGA/TAG/TAA = X)','Universal\t(TGA/TAG/TAA = STOP)',\
	'TAA\t\t(TAG/TGA = Q)', 'TAG\t\t(TRA = Q)', 'TGA\t\t(TAR = Q)']

	author = (color.BOLD+color.ORANGE+'\n\n\tQuestions/Comments? Email Xyrus (author) at'\
	' maurerax@gmail.com\n\n'+color.END)


	if args.genetic_code != None and args.genetic_code.lower() not in supported_gcodes_names:
		print (color.BOLD+color.RED+'\nProvided genetic code is currently unsupported.\n\n'\
		'If you have a new genetic code, please contact the author (with some evidence).\n\n'\
		'Otherwise, use one of the currently supported genetic codes.\n'+color.END)
		print (color.BOLD+color.ORANGE+'\n'.join(supported_gcodes_list)+'\n\n'+color.END)
		print (author)
		valid_arg += 1

	if args.list_codes == True:
		print (color.BOLD+color.RED+'\nThese are the currently supported genetic codes.\n'+color.END)
		print (color.BOLD+color.ORANGE+'\n'.join(supported_gcodes_list)+'\n\n'+color.END)
		valid_arg += 1
		sys.exit()

	if args.author == True:
		print (author)
		valid_arg += 1
		sys.exit()

	if args.input_file != None:
		if os.path.isfile(args.input_file) != False:
			if args.input_file.split('/')[-1] not in os.listdir('/'.join(args.input_file.split('/')[:-1])):
				print (color.BOLD+color.RED+'\nError:'+color.END+color.BOLD+' The provided Fasta file '\
				'('+color.DARKCYAN+args.input_file.split('/')[-1]+color.END+color.BOLD+')\ndoes not'\
				' exist or is incorrectly formatted.\n\nDouble-check then try again!\n\n'+color.END)
				valid_arg += 1
			elif args.input_file.endswith('WTA_NBU.Renamed.fasta') != True:
				print (color.BOLD+'\n\nInvalid Fasta File! Only Fasta Files that were processed'\
				' with '+color.GREEN+'3_CountOGsUsearch.py '+color.END+color.BOLD+'are valid\n\n'\
				'However, to bypass that issue, Fasta Files MUST end with '+color.CYAN+\
				'"WTA_NBU.Renamed.fasta"\n\n'+color.END)
				valid_arg += 1
		else:
			print (color.BOLD+color.RED+'\nError:'+color.END+color.BOLD+' The provided Fasta file '\
			'('+color.DARKCYAN+args.input_file.split('/')[-1]+color.END+color.BOLD+')\ndoes not'\
			' exist or is incorrectly formatted.\n\nDouble-check then try again!\n\n'+color.END)
			valid_arg += 1

	return valid_arg


###########################################################################################
###--------------------------- Does the Inital Folder Prep -----------------------------###
###########################################################################################

def prep_folders(args):

	OG_folder = '../'+args.input_file.split('/')[1]+'/UsearchOG/'

	if os.path.isdir(OG_folder) != True:
		os.system('mkdir '+OG_folder)

	if os.path.isdir(OG_folder+'OG_CompData') != True:
		os.system('mkdir '+OG_folder+'OG_CompData')

	if os.path.isdir('../TranslatedTranscriptomes') != True:
		os.system('mkdir ../TranslatedTranscriptomes')



##########################################################################################
###---------------- Scans 5-Prime End of Transcript for In-Frame "ATG" ----------------###
##########################################################################################

def check_new_start_new(some_seq, low_lim, upper_lim, old_start, codon_table):

	## Looks for in-frame STOP codons in the UTR of the transcript
	prime5 = str(Seq(some_seq[low_lim:upper_lim]).translate(table=codon_table)).replace('*','x')
	in_frame_stops = [stops.start() for stops in re.finditer('x',prime5)]

	## Looks for in-frame START codons in the UTR of the transcript
	in_frame_starts = [starts.start() for starts in re.finditer('M',prime5)]

	## Checks that there are NO in-frame STOP codons between the possible "new" START codon
	## and the aligned portion of the transcript -- THIS is double checked!
	if len(in_frame_starts) != 0:
		if len(in_frame_stops) != 0:
			if in_frame_stops[-1] < in_frame_starts[-1]:
				new_start = low_lim+in_frame_starts[-1]*3
			else:
				new_start = old_start
		else:
			new_start = low_lim+in_frame_starts[-1]*3
	else:
		new_start = old_start

	## Skips the double-checking if there are no GOOD potential START codons
	if new_start == old_start:
		updated_start = old_start

	else:
	## Double checks that there are NO IN-FRAME stop codons between the NEW-SUGGESTED Start
	## position and the OLD-SUPPORTED stop position!
		between_new_old_start = str(Seq(some_seq[new_start:old_start]).translate(table=1)).replace('*','x')
		in_frame_stops_check = [stops.start() for stops in re.finditer('x',between_new_old_start)]
		in_frame_starts_check = [starts.start() for starts in re.finditer('M',between_new_old_start)]
		if len(in_frame_starts_check) != 0:
			if len(in_frame_stops_check) != 0:
				if in_frame_stops_check[-1] < in_frame_starts_check[-1]:
					updated_start = new_start+in_frame_starts_check[-1]*3
				else:
					updated_start = old_start
			else:
				updated_start = new_start
		else:
			updated_start = new_start

	return updated_start


##########################################################################################
###--------------- Extracts the ORF from the Fasta File and SpreadSheet ---------------###
##########################################################################################

def extract_ORF(prot_dict, codon_table, args):

	print (color.BOLD+'\n\nExtracting '+color.PURPLE+'ORFs'+color.END+color.BOLD+' from'\
	' the transcriptomic data-set\n\n'+color.END)

	for k, v in prot_dict.items():

	## Attempting to find the most-likely START (ATG) position in the transcript (tricky)
	## Skips this if the initial Methionine (ATG) is likely present
	## (e.g. the alignment position of the protein = '1')
		prot_start = int(v[3].split('..')[0])
		old_start = v[1]
		if prot_start != 1:
			min_dist, max_dist = round_down_three(prot_start)
			min_start = old_start-min_dist
			max_start = old_start-max_dist
			if min_start < 0:
				min_start = old_start
			if max_start < 0:
				max_start = min_start%3
#			print k+'\tOld_start\t'+str(old_start)+'\tMin_Dist/Start\t'+str(min_dist)+'/'+str(min_start)+'\tMax_Dist/Start\t'+str(max_dist)+'/'+str(max_start)+'\n'
			updated_start = check_new_start_new(v[-1], max_start, min_start, old_start, codon_table)
		else:
			updated_start = old_start
		temp = prot_dict[k][-1][updated_start:]

	## Uses the given genetic code to identify the stop position of the ORF
		temp_prot = str(Seq(temp).translate(table=codon_table))
		if '*' in temp_prot:
			stop_pos = (temp_prot.index('*')+1)*3
			prot_dict[k].append(temp[:stop_pos])
		else:
			stop_pos = prot_dict[k][2] - prot_dict[k][1]
			prot_dict[k].append(temp[:stop_pos])

	## Awkward_list is populated with unexpectedly SHORT ORFs!
	## Reasons for being short include:
		# An error Xyrus introduced
		# Not as great genetic code decision (in-frame stop)
		# Crummy sequence/assembly quality (false in-frame stop codons)

	awkward_list = []
	look_good = []

	for k, v in prot_dict.items():
		expected_min = len(v[-2][v[1]:v[2]])-1
		if len(v[-1]) < expected_min:
			awkward_list.append(k)
		else:
			look_good.append(k)

	if len(awkward_list) != 0:
		with open('UnexpexctedShortStuffBlameXyrus.txt','w+') as x:
			for entry in awkward_list:
				x.write(entry+'\n')
	else:
		pass

	print (color.BOLD+'\n\nTranslating '+color.PURPLE+'ORFs'+color.END+color.BOLD+' from'\
	' using the '+color.DARKCYAN+args.genetic_code.title()+' genetic code'+color.END)

	for k, v in prot_dict.items():
		prot_dict[k].append(str(Seq(v[-1]).translate(table=codon_table)).rstrip('*'))

	return prot_dict

##########################################################################################
###------------ Grabs the Coding Coordinates from the OG-BLAST SpreadSheet ------------###
##########################################################################################

def prep_translations(args):

	print (color.BOLD+'\n\nGrabbing useful info from the '+color.ORANGE+args.input_file\
	.split('/')[-1]+color.END+color.BOLD+' Fasta File\nand from the '+color.ORANGE+args.tsv_file.\
	split('/')[-1]+color.END+color.BOLD+' OG-Assignment Spreadsheet'+color.END)

	inTSV = [i.rstrip('\n') for i in open(args.tsv_file).readlines() if i != '\n']
	inFasta = [i for i in SeqIO.parse(args.input_file,'fasta')]

	# ORF identification step here, uses the 'allOGCleanresults.tsv file to identify the ORF
	prot_dict = {}

	# Special scenario! Only for when the genetic code is not particularly useful ...
	if args.genetic_code.lower() == 'none' or args.genetic_code.lower() == 'condylostoma'\
	or args.genetic_code.lower() == 'condy':
		for i in inTSV:
			prot_dict.setdefault(i.split('\t')[0],[])
			if int(i.split('\t')[6]) < int(i.split('\t')[7]):
	## Saves the Transcript Orientation (Coding vs. Template Strand)
				prot_dict[i.split('\t')[0]].append('F')
	## Collects initial Start and Stop positions from the BLAST alignment
				prot_dict[i.split('\t')[0]].append(int(i.split('\t')[6])-1)
				prot_dict[i.split('\t')[0]].append(int(i.split('\t')[7])+3)
	## Implied Amino Acid alignment positions (e.g. does the alignment start at the 1st Methionine?)
				prot_dict[i.split('\t')[0]].append('..'.join(i.split('\t')[-4:-2]))

			if int(i.split('\t')[7]) < int(i.split('\t')[6]):
	## Saves the Transcript Orientation (Coding vs. Template Strand)
				prot_dict[i.split('\t')[0]].append('RC')
	## Collects initial Start and Stop positions from the BLAST alignment
				prot_dict[i.split('\t')[0]].append(int(i.split('_Len')[1].split('_')[0])-int(i.split('\t')[6]))
				prot_dict[i.split('\t')[0]].append(int(i.split('_Len')[1].split('_')[0])-int(i.split('\t')[7])+1)
	## Implied Amino Acid alignment positions (e.g. does the alignment start at the 1st Methionine?)
				prot_dict[i.split('\t')[0]].append('..'.join(i.split('\t')[-4:-2]))

	## Makes sure that the dictionary has the transcript in the correct orientation
		for i in inFasta:
			if i.description in prot_dict.keys():
				if 'RC' == prot_dict[i.description][0]:
					prot_dict[i.description].append(str(i.seq.reverse_complement()))
				else:
					prot_dict[i.description].append(str(i.seq))

	else:
		for i in inTSV:
			prot_dict.setdefault(i.split('\t')[0],[])
			if int(i.split('\t')[6]) < int(i.split('\t')[7]):
	## Saves the Transcript Orientation (Coding vs. Template Strand)
				prot_dict[i.split('\t')[0]].append('F')
				prot_dict[i.split('\t')[0]].append(int(i.split('\t')[6])-1)
				prot_dict[i.split('\t')[0]].append(int(i.split('\t')[7])+3)
	## Implied Amino Acid alignment positions (e.g. does the alignment start at the 1st Methionine?)
				prot_dict[i.split('\t')[0]].append('..'.join(i.split('\t')[-4:-2]))
			if int(i.split('\t')[7]) < int(i.split('\t')[6]):
	## Saves the Transcript Orientation (Coding vs. Template Strand)
				prot_dict[i.split('\t')[0]].append('RC')
	## Collects initial Start and Stop positions from the BLAST alignment (but in the "correct" orientation)
				prot_dict[i.split('\t')[0]].append(int(i.split('_Len')[1].split('_')[0])-int(i.split('\t')[6]))
				prot_dict[i.split('\t')[0]].append(int(i.split('_Len')[1].split('_')[0])-int(i.split('\t')[7])+1)
	## Implied Amino Acid alignment positions (e.g. does the alignment start at the 1st Methionine?)
				prot_dict[i.split('\t')[0]].append('..'.join(i.split('\t')[-4:-2]))

	## Makes sure that the dictionary has the transcript in the correct orientation
		for i in inFasta:
			if i.description in prot_dict.keys():
				if 'RC' == prot_dict[i.description][0]:
					prot_dict[i.description].append(str(i.seq.reverse_complement()))
				else:
					prot_dict[i.description].append(str(i.seq))

	return prot_dict


##########################################################################################
###------------------------ Rounds Down Values to Nearest "3" -------------------------###
##########################################################################################

def round_down_three(num):
	min_val = int(num*3*.5)-int(num*3*.5)%3
	max_val = int(num*6)-int(num*6)%3
	return min_val, max_val


##########################################################################################
###--------------------- Makes Translation Steps (Later) Easier -----------------------###
##########################################################################################

def standardize_gcode(given_code):
	if given_code == 'ciliate' or given_code == 'tga':
		codon_table = 6
	elif given_code == 'chilodonella' or given_code == 'chilo' or given_code == 'taa':
		codon_table = c_uncinata_table
	elif given_code == 'blepharisma' or given_code == 'bleph':
		codon_table = blepharisma_table
	elif given_code == 'euplotes' or given_code == 'eup':
		codon_table = euplotes_table
	elif given_code == 'myrionecta' or given_code == 'mesodinium':
		codon_table = myrionecta_table
	elif given_code == 'peritrich' or given_code == 'vorticella':
		codon_table = peritrich_table
	elif given_code == 'none':
		codon_table = no_stop_table
	elif given_code == 'condylostoma' or given_code == 'condy':
		codon_table = condylostoma_table
	elif given_code == 'tag':
		codon_table = tag_table
	elif given_code == 'universal':
		codon_table = 1
	else:
		print (color.BOLD+color.RED+'\n\nNo valid genetic code provided!\n\n'+color.END+\
		color.BOLD+'Using the "Universal" genetic code (by default)\n\nPlease check that the'\
		' code you wish to use is supported:'+color.CYAN+'\n\npython 5_GCodeTranslate.py'\
		' -list_codes\n\n'+color.END)
		codon_table = 1

	return codon_table


###########################################################################################
###------------------ Updates Spreadsheet with Updated Contig Names --------------------###
###########################################################################################

def update_spreadsheet(args, updated_spreadsheet_dict):
	if os.path.isdir('../'+args.input_file.split('/')[1]+'/UsearchOG') != True:
		os.system('mkdir ../'+args.input_file.split('/')[1]+'/UsearchOG/')
	else:
		pass

# 	os.system('cp '+args.tsv_file+' ../'+args.tsv_file.split('/')[1]+'/UsearchOG/'+args.tsv_file\
# 	.split('/')[-1])

	inTSV = [line.rstrip('\n') for line in open(args.tsv_file).readlines() if line != '\n'\
	 and line.split('\t')[0] in updated_spreadsheet_dict.keys()]

	updatedTSV = [updated_spreadsheet_dict[line.split('\t')[0]]+'\t'+'\t'.join(line.split('\t')[1:]) for line in inTSV]

	with open(args.tsv_out,'w+') as w:
		w.write('\n'.join(updatedTSV))


###########################################################################################
###-------------------- Updates Log With OG Assignment Information ---------------------###
###########################################################################################

def update_log(filename, codon_table):

	if os.path.isdir('../PostAssembly_Logs/') != True:
		os.system('mkdir ../PostAssembly_Logs/')
	else:
		pass

	ntd_ORF = [i for i in SeqIO.parse(filename.split('.fas')[0]+'_'+gcode.title()+'_ORF.fasta','fasta')]
	aa_ORF = [i for i in SeqIO.parse(filename.split('.fas')[0]+'_'+gcode.title()+'_ORF.aa.fasta','fasta')]

	min_ntd_ORF = str(min([len(i.seq) for i in ntd_ORF]))
	max_ntd_ORF = str(max([len(i.seq) for i in ntd_ORF]))
	avg_ntd_ORF = '%.2f' % (sum([len(i.seq) for i in ntd_ORF])/float(len(ntd_ORF)))

	min_aa_ORF = str(min([len(i.seq) for i in aa_ORF]))
	max_aa_ORF = str(max([len(i.seq) for i in aa_ORF]))
	avg_aa_ORF = '%.2f' % (sum([len(i.seq) for i in aa_ORF])/float(len(aa_ORF)))

	for Logname in os.listdir(os.curdir+'./PostAssembly_Logs/'):
		if Logname.startswith(filename.split('/')[2].split('_WTA')[0]) and Logname.endswith('Log.txt'):
			with open('../PostAssembly_Logs/'+Logname,'a') as LogFile:
				LogFile.write('Nucleotide ORFs\t'+str(len(ntd_ORF))+'\tn/a\tn/a\n')
				LogFile.write('Nucleotide ORF Lengths\t'+avg_ntd_ORF+'\t'+min_ntd_ORF+'\t'+max_ntd_ORF+'\n')
				LogFile.write('Protein ORFs\t'+str(len(aa_ORF))+'\tn/a\tn/a\n')
				LogFile.write('Protein ORF Lengths\t'+avg_aa_ORF+'\t'+min_aa_ORF+'\t'+max_aa_ORF+'\n')


##########################################################################################
###----------------------- Write File with Provided Genetic Code ----------------------###
##########################################################################################

def write_data_out(prot_dict, codon_table, args):

	update_spreadsheet_dict = {}

	for k, v in prot_dict.items():
		if 'Cov' in k:
			new_name = k.split('_Len')[0]+'_Len'+str(len(v[-2]))+'_'+'_'.join(k.split('_')[-3:])
			update_spreadsheet_dict[k] = new_name
		else:
			new_name = k.split('_Len')[0]+'_Len'+str(len(v[-2]))+'_'+'_'.join(k.split('_')[-2:])
			update_spreadsheet_dict[k] = new_name


	with open(args.ntd_out,'w+') as w:
		print (color.BOLD+'\n\nWriting FASTA file with '+color.PURPLE+'ORF'+color.END+color.BOLD\
		+' sequences using the '+color.DARKCYAN+args.genetic_code.title()+' genetic code'+color.END)

		for k, v in prot_dict.items():
			w.write('>'+update_spreadsheet_dict[k]+'\n'+str(v[-2])+'\n')

	with open(args.aa_out, 'w+') as w:
		print (color.BOLD+'\n\nWriting FASTA file with '+color.PURPLE+'Translated ORF'+color.END+color.BOLD\
		+' sequences using the '+color.DARKCYAN+args.genetic_code.title()+' genetic code'+color.END)

		for k, v in prot_dict.items():
			w.write('>'+update_spreadsheet_dict[k]+'\n'+str(v[-1])+'\n')

	return update_spreadsheet_dict


##########################################################################################
###--------------------- Cleans up the Folder and Moves Final Files -------------------###
##########################################################################################

def clean_up(args):

	if args.input_file.split('.fas')[0].split('/')[-1]+'_StopCodonStats.tsv' in os.listdir(args.home_folder):
		os.system('mv '+args.input_file.split('.fas')[0]+'_StopCodonStats.tsv '+args.StopFreq)

	os.system('mv '+args.tsv_file+' '+args.Usearch_Folder)
	os.system('mv '+args.input_file+' '+args.Usearch_Folder)

	if args.no_RP == True:
		if os.path.isdir(args.ForPartials+'ToRename/') != True:
			os.system('mkdir '+args.ForPartials+'ToRename/')

		os.system('cp '+args.ntd_out+' '+args.ForPartials+'ToRename/')
		os.system('cp '+args.aa_out+' '+args.ForPartials+'ToRename/')
		os.system('cp '+args.tsv_out+' '+args.ForPartials+'ToRename/')

	else:
		os.system('cp '+args.tsv_out+' '+args.ForPartials)
		os.system('cp '+args.ntd_out+' '+args.ForPartials)
		os.system('cp '+args.aa_out+' '+args.ForPartials)

	os.system('mv '+args.home_folder+' ../TranslatedTranscriptomes')


###########################################################################################
###-------------------------------- Next Script Message --------------------------------###
###########################################################################################

def next_script(args):

	print (color.BOLD+'\n\nLook for '+color.DARKCYAN+args.ntd_out.split('/')[-1]+color.END+\
	color.BOLD+',\n'+color.DARKCYAN+args.aa_out.split('/')[-1]+color.END+color.BOLD+', and\n'\
	+color.DARKCYAN+args.tsv_out.split('/')[-1]+color.END+color.BOLD+',\nwhich are in the '+\
	color.ORANGE+args.home_folder.split('/')[-1]+' Folder'+color.END)

	if args.no_RP == True:
		print(color.BOLD+'\n\nNext Script is: '+color.GREEN+'7_FinalRename.py'+color.END+color.BOLD+\
		' in the '+color.PURPLE+'RemovePartials Folder'+color.END+color.BOLD+'\nwith a copy of'\
		' the outputs of this script!'+color.END)
		print(color.BOLD+'\n\nRemember that you have chosen '+color.RED+'NOT '+color.END+color.BOLD+\
		'to remove partials\nand are skipping to the renaming step!\n\n'+color.END)

	else:
		print(color.BOLD+'\n\nNext Script is: '+color.GREEN+'6_FilterPartials.py'+color.END+color.BOLD+\
		' in the '+color.PURPLE+'FinalizeTranscripts Folder'+color.END+color.BOLD+'\nwith a copy of'\
		' the outputs of this script!\n\n'+color.END)



##########################################################################################
###--------------- Checks Command Line Arguments and Calls on Functions ---------------###
##########################################################################################

def main():

	args = check_args()

	prep_folders(args)

	codon_table = standardize_gcode(args.genetic_code.lower())

	prot_dict_Prepped = prep_translations(args)

	prot_dict_Final = extract_ORF(prot_dict_Prepped, codon_table, args)

	new_spreadsheet_names = write_data_out(prot_dict_Final, codon_table, args)

	update_spreadsheet(args, new_spreadsheet_names)

#	update_log(fasta_file, gcode)

	clean_up(args)

	next_script(args)

main()
