#!/usr/bin/env python3.5

##__Updated__: 18_08_2017
##__Author__: Xyrus Maurer-Alcala; maurerax@gmail.com
##__Usage__: python 2_Auto_rRNA_BvE.py --help

##########################################################################################
## This script is intended to identify and isolate SSU/LSU sequences and then classify	##
## the output as STRONGLY Bacterial or STRONGLY Eukaryotic or UNDETERMINED/UNKNOWN		##
## Prior to running this script, ensure the following:									##
##																						##
## 1. You have assembled your transcriptome and COPIED the 'assembly' file 				##
##    (contigs.fasta, or scaffolds.fasta) to the PostAssembly Folder					##
## 2. Removed small sequences (usually sequences < 300bp) with ContigFilterPlusStats.py	##
##																						##
## 								COMMAND Example Below									##
##																						##
## 			E-mail Xyrus (author) for help if needed: maurerax@gmail.com				##
##																						##
##							Next Script(s) to Run: 										##
##							  3_CountOGsUsearch.py										##
##																						##
##########################################################################################


import argparse, os, sys
from argparse import RawTextHelpFormatter,SUPPRESS


#------------------------------ Colors For Print Statements ------------------------------#

class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   ORANGE = '\033[38;5;214m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'



#------------------------------- Main Functions of Script --------------------------------#

###########################################################################################
###--------------------- Parses and Checks Command-Line Arguments ----------------------###
###########################################################################################
	
def check_args():

	parser = argparse.ArgumentParser(description=
	color.BOLD+'\nThis script will remove '+color.RED+'rDNA contigs (both SSU and LSU)'+color.END\
	+color.BOLD+'\nfrom your Assembly using a set of '+color.RED+'SSU/LSU rDNAs'+color.END\
	+color.BOLD+' from diverse\n'+color.ORANGE+'Eukaryotes, Bacteria and Archaea'\
	+color.END+color.BOLD+'.\n\nThen, this script will categorize non-rRNA'\
	'sequences as'+color.ORANGE+'\nSTRONGLY '+color.END+color.BOLD+color.RED+'Eukaryotic'\
	' OR Prokaryotic'+color.END+color.BOLD+' using a set of Proteins\nfrom diverse'+\
	color.ORANGE+'Eukaryotes, Bacteria and Archaea'+color.END+color.BOLD+'.'+color.END\
	+usage_msg(), usage=SUPPRESS,formatter_class=RawTextHelpFormatter)
	
	optional_arg_group = parser.add_argument_group(color.ORANGE+color.BOLD+'Options'+color.END)

	optional_arg_group.add_argument('--input_file','-in', action='store',
	help=color.BOLD+color.GREEN+'Fasta file of Nucleotide sequences\n'+color.END)
	optional_arg_group.add_argument('--threads','-t', default='2',
	help=color.BOLD+color.GREEN+' Number of threads to use for BLAST\n (default = 2)\n'+color.END)
	optional_arg_group.add_argument('-author', action='store_true',
	help=color.BOLD+color.GREEN+' Print author contact information\n'+color.END)

	infile_list = []
	
	if len(sys.argv[1:]) == 0:
		if len([i for i in os.listdir(os.curdir+'./') if i.endswith('bp.fasta')]) > 0:
			infile_list += ['../'+i for i in os.listdir(os.curdir+'./') if i.endswith('bp.fasta')]
		else:		
			print (parser.description)
			print ('\n')
			sys.exit()

	elif '-author' in sys.argv[1:]:
		print (parser.description)
		print (color.BOLD+color.ORANGE+'\n\n\tQuestions/Comments? Email Xyrus (author) at maurerax@gmail.com\n\n'+color.END)
		sys.exit()

	elif '-h' in sys.argv[1:] or '--help' in sys.argv[1:]:
			pass
			args = parser.parse_args()
	
	args = parser.parse_args()
	
	if args.input_file != None:
		if args.input_file.endswith('bp.fasta') == True and args.input_file.strip('../') in os.listdir(os.curdir+'./'):
			infile_list += [args.input_file]
		else:
			print (color.BOLD + '\nCheck that you are giving an appropriately Named/Processed'\
			'Fasta file(s) to this script\n\nNOTE that this script CURRENTLY expects your'\
			' Fasta files to contain '+color.RED+ '"rna"'+color.END+color.BOLD+' in \nthe Fasta File'\
			' Name and must end with ' + color.RED + '"bp.fasta"' + color.END)
			print (usage_msg()+'\n')
			sys.exit()
	else:
		if infile_list[0] == '':
			print (color.BOLD+'\nNo valid input files were provided/exist\n\nCheck that '\
			'the outputs from 1_ContigFiltStats.py are in the main '+color.RED+'Transcriptome'\
			+color.END+color.BOLD+' folder before trying again\n\n'+color.END)
			sys.exit()

	return infile_list, args
	
	
###########################################################################################
###------------------------------- Script Usage Message --------------------------------###
###########################################################################################

def usage_msg():
	return color.BOLD+color.RED+'\n\n\nExample usage (For Single File):'+color.CYAN+' python'\
	' 2_Auto_rRNA_BvE.py --input_file ../Op_me_Xxma_rna.200bp.fasta'+color.END+color.BOLD+\
	color.RED+'\n\nExample usage (For Numerous Files):'+color.CYAN+' python 2_Auto_rRNA_BvE.py'\
	+color.END


###########################################################################################
###---------------------- Removes rRNA Sequences from Fasta Files ----------------------###
###########################################################################################

def remove_rRNA(SizeFilt_Fasta, args):
	os.system('python 2a_remove_rRNA.py --input_file '+SizeFilt_Fasta+' --threads '+args.threads)
	

###########################################################################################
###-------------- Compares Contigs Against Euk and Prok Protein Databases --------------###
###########################################################################################

def remove_Bact(rRNA_Free_Fasta, args):
	fasta_withBact = rRNA_Free_Fasta.split('/')[-1].split('_rna')[0]+'_NorRNAseqs.fasta'
	fasta_folder = rRNA_Free_Fasta.split('_rna')[0]
	os.system('python 2b_remove_Bact.py --input_file '+fasta_folder+'/'+fasta_withBact)

##########################################################################################
###--------------- Checks Command Line Arguments and Calls on Functions ---------------###
##########################################################################################

def main():

	Fasta_FileList, args = check_args()

	for Fasta_File in Fasta_FileList:
		remove_rRNA(Fasta_File, args)
		remove_Bact(Fasta_File, args)	

main()