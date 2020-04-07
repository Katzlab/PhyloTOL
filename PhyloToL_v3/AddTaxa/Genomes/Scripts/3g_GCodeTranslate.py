#!/usr/bin/env python3.5

##__Updated__: 19_09_2017
##__Author__: Xyrus Maurer-Alcala; maurerax@gmail.com
##__Usage__: python 3g_GCodeTranslate.py --help


##############################################################################
##																			##
## Translates CDSs sequences using the Provided Genetic Code. 				##
##																			##
## NOTE: 																	##
##		No provided input for genetic code results in Translation with the	##
## 		UNIVERSAL genetic code (as default)									##
##																			##
##		E-mail Xyrus (author) for help if needed: maurerax@gmail.com		##
##																			##
##############################################################################


import argparse, os, sys
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


#----------------------------- Colors For Print Statements ------------------------------#
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
		
	args.Prepped_Folder = '../'+args.input_file.split('/')[1]+'/Prepped'
	args.folder = '../'+args.input_file.split('/')[1]
	args.out_name = args.input_file.split('.Prepped')[0]+'.'+args.genetic_code.title()+'.AA.fasta'
	args.new_ntd_name = args.input_file.split('.Prepped')[0]+'.'+args.genetic_code.title()+'.NTD.fasta'
	
	return args
	
	
###########################################################################################
###------------------------------- Script Usage Message --------------------------------###
###########################################################################################

def usage_msg():
	return (color.BOLD+color.RED+'\n\nExample usage:'+color.CYAN+' python 3g_GCodeTranslate.py'\
	' --input_file ../Stentor_coeruleus.WGS.CDS.Prep/Stentor_coeruleus.WGS.CDS.Prepped.fasta'\
	' --genetic_code Universal'+color.END)


##########################################################################################
###-------- Storage for LARGE (Annoying) Print Statements for Flagged Options ---------###
##########################################################################################

def return_more_info(args):

	valid_arg = 0

	supported_gcodes_names = ['bleph','blepharisma','chilo','chilodonella','condy',\
	'condylostoma','none','eup','euplotes','peritrich','vorticella','ciliate','universal',\
	'taa','tag','tga']

	supported_gcodes_list = ['Blepharisma\t(TGA = W)','Chilodonella\t(TAG/TGA = Q)','Ciliate\t\t(TAR = Q)',\
	'Conylostoma\t(TAR = Q, TGA = W)','Euplotes\t(TGA = C)','Peritrich\t(TAR = E)','None\t\t(TGA/TAG/TAA = X)',\
	'Universal\t(TGA/TAG/TAA = STOP)','TAA\t\t(TAG/TGA = Q)', 'TAG\t\t(TRA = Q)', 'TGA\t\t(TAR = Q)']

	author = (color.BOLD+color.ORANGE+'\n\n\tQuestions/Comments? Email Xyrus (author) at'\
	' maurerax@gmail.com\n\n'+color.END)


	if args.genetic_code != None and args.genetic_code.lower() not in supported_gcodes_names:
		print (color.BOLD+color.RED+'\nProvided genetic code is currently unsupported.\n\n'\
		'If you have a new genetic code, please contact the author (with some evidence).\n\n'\
		'Otherwise, use one of the currently supported genetic codes.\n'+color.END)
		print (color.BOLD+color.ORANGE+'\n'.join(supported_gcodes_list)+'\n\n'+color.END)
		print (author)
		valid_arg += 1
	else:	
		if args.list_codes == True:
			print (color.BOLD+color.RED+'\nThese are the currently supported genetic codes.\n'+color.END)
			print (color.BOLD+color.ORANGE+'\n'.join(supported_gcodes_list)+'\n\n'+color.END)
			valid_arg += 1	

		if args.author == True:
			print (author)
			valid_arg += 1

	if args.input_file != None:
		if os.path.isfile(args.input_file) != False:
			if args.input_file.split('/')[-1] not in os.listdir('/'.join(args.input_file.split('/')[:-1])):
				print (color.BOLD+color.RED+'\nError:'+color.END+color.BOLD+' The provided Fasta file '\
				'('+color.DARKCYAN+args.input_file.split('/')[-1]+color.END+color.BOLD+')\ndoes not'\
				' exist or is incorrectly formatted.\n\nDouble-check then try again!\n\n'+color.END) 
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

	if os.path.isdir(args.Prepped_Folder) != True:
		os.system('mkdir '+args.Prepped_Folder)
		os.system('cp '+args.input_file+' '+args.Prepped_Folder)
		
	if os.path.isdir(args.folder+'/Translated') != True:
		os.system('mkdir '+args.folder+'/Translated')

##########################################################################################
###------------------ Translates CDSs from the Provided Genetic Code ------------------###
##########################################################################################
	
def translate_seqs(args):
	
	inFasta = [i for i in SeqIO.parse(args.input_file,'fasta')]
	
	print (color.BOLD+'\n\n\nTranslating: '+color.CYAN+args.input_file.split('/')[-1]+color.END+\
	color.BOLD+'\nwith the '+color.GREEN+args.genetic_code.upper()+' Genetic Code\n'+color.END)

	
	if args.genetic_code.lower() == 'ciliate' or args.genetic_code.lower() == 'tga':
		translated_seqs = ['>'+seq_rec.description+'\n'+str(seq_rec.seq.translate(table=6)).rstrip('*').replace('*','X')+'\n' for seq_rec in inFasta]

	if args.genetic_code.lower() == 'peritrich' or args.genetic_code.lower() == 'vorticella':
		translated_seqs = ['>'+seq_rec.description+'\n'+str(seq_rec.seq.translate(table=peritrich_table)).rstrip('*').replace('*','X')+'\n' for seq_rec in inFasta]

	if args.genetic_code.lower() == 'tag':
		translated_seqs = ['>'+seq_rec.description+'\n'+str(seq_rec.seq.translate(table=tag_table)).rstrip('*').replace('*','X')+'\n' for seq_rec in inFasta]

	if args.genetic_code.lower() == 'chilo' or args.genetic_code.lower() == 'chilodonella' or args.genetic_code.lower() == 'taa':
		translated_seqs = ['>'+seq_rec.description+'\n'+str(seq_rec.seq.translate(table=c_uncinata_table)).rstrip('*').replace('*','X')+'\n' for seq_rec in inFasta]

	if args.genetic_code.lower() == 'bleph' or args.genetic_code.lower() == 'blepharisma':
		translated_seqs = ['>'+seq_rec.description+'\n'+str(seq_rec.seq.translate(table=blepharisma_table)).rstrip('*').replace('*','X')+'\n' for seq_rec in inFasta]

	if args.genetic_code.lower() == 'eup' or args.genetic_code.lower() == 'euplotes':
		translated_seqs = ['>'+seq_rec.description+'\n'+str(seq_rec.seq.translate(table=euplotes_table)).rstrip('*').replace('*','X')+'\n' for seq_rec in inFasta]

	if args.genetic_code.lower() == 'universal':
		translated_seqs = ['>'+seq_rec.description+'\n'+str(seq_rec.seq.translate(table=1)).rstrip('*').replace('*','X')+'\n' for seq_rec in inFasta]
	
	return translated_seqs


##########################################################################################
###---------------------------- Writes Out Translated CDSs ----------------------------###
##########################################################################################

def write_out(args):

	translated_seqs = translate_seqs(args)
	
	## Keep only ORFs greater than 10 amino acids long
	translated_seqs = [i for i in translated_seqs if len(i.split('\n')[1]) > 10]
	
	print (color.BOLD+'\nTranslated '+color.ORANGE+str(len(translated_seqs))+color.END\
	+color.BOLD+' seqeunces using the '+color.GREEN+args.genetic_code.upper()+' Genetic Code\n\n'+color.END)
	
	with open(args.out_name,'w+') as w:
		w.write(''.join(translated_seqs))


##########################################################################################
###--------------------- Cleans up the Folder and Moves Final Files -------------------###
##########################################################################################

def clean_up(args):
	
	os.system('cp '+args.input_file+' '+args.Prepped_Folder)

	os.system('mv '+args.input_file+' '+args.new_ntd_name)

	os.system('cp '+args.out_name.split('.AA')[0]+'* '+args.folder+'/Translated/')
	
	os.system('mv '+args.input_file.split('.fasta')[0]+'.GeneticCode.txt '+args.folder+'/Translated/')

	
##########################################################################################
###----------------------------- Calls on Above Functions -----------------------------###
##########################################################################################

def main():
	
	args = check_args()
	
	prep_folders(args)
	
	write_out(args)
	
	clean_up(args)
		
	print (color.BOLD+'Next Script is: '+color.PURPLE+' 4g_CountOgsUsearch.py\n\n'+color.END)
	
main()