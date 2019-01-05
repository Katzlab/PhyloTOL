#!/usr/bin/env python

##__Updated__: 18_08_2017
##__Author__: Xyrus Maurer-Alcala; maurerax@gmail.com
##__Usage__: python 4_InFrameStopFreq.py --help


##########################################################################################
## This script is intended to aid in identifying the genetic code of the data given		##
##																						##
## Prior to running this script, ensure the following:									##
##																						##
## 1. You have assembled your transcriptome and COPIED the 'assembly' file 				##
##    (contigs.fasta, or scaffolds.fasta) to the PostAssembly Folder					##
## 2. Removed small sequences (usually sequences < 300bp) with ContigFilterPlusStats.py	##
## 3. Removed SSU/LSU sequences from your Fasta File									##
## 4. Classified your sequences as Strongly Prokaryotic/Eukaryotic or Undetermined		##
## 5. Classified the Non-Strongly Prokaryotic sequences into OGs 						##
##																						##
## 								COMMAND Example Below									##
##							Extra Notes at Bottom of Script								##
##																						##
## 			E-mail Xyrus (author) for help if needed: maurerax@gmail.com				##
##																						##
##								Next Script(s) to Run: 									##
##							 	 5_GCodeTranslate.py									##
##																						##
##########################################################################################


import argparse, os, sys
from argparse import RawTextHelpFormatter,SUPPRESS
from distutils import spawn

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data.CodonTable import CodonTable


#-------------------------- Set-up Codon Tables (Genetic Codes) --------------------------#

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
###---------------------------- UPDATE USEARCH PATH BELOW! -----------------------------###
###########################################################################################
	## IF Usearch is IN YOUR PATH then no updating is needed...

def check_usearch_path():

	usearch_path = ''

	if usearch_path == '':
		usearch_path = spawn.find_executable("usearch")
	else:
		pass

	if usearch_path == None:
		print (color.BOLD + '\n\nPlease open this script and check that you have included'\
		+' the PATH to the'+color.BLUE+' "usearch" '+color.END+color.BOLD+'executable.\n\n'+color.END)
		print (color.BOLD+color.BLUE+'LOOK FOR:\n\n'+color.RED\
		+'#------------------------------ UPDATE USEARCH PATH BELOW! -------------------------------#'\
		+color.BLUE+'\n\nThis is somewhere around lines 50 - 80...\n\n'+color.END)

		sys.exit()
	else:
		pass
	
	return usearch_path


###########################################################################################
###--------------------- Parses and Checks Command-Line Arguments ----------------------###
###########################################################################################

def check_args():

	parser = argparse.ArgumentParser(description=
	color.BOLD+'\n\nThis script is intended to '+color.RED+'AID You '+color.END+color.BOLD\
	+'in determining the '+color.RED+'\nLikely Genetic Code'+color.END+color.BOLD+' of a'\
	' given Fasta File of transcripts\n\nInterpretation of the output (StopFreq.tsv) is difficult \nand so '+color.ORANGE\
	+'TWO EXAMPLES'+color.END+color.BOLD+' can be found in the '+color.CYAN+'NOTES Section'\
	+color.END+color.BOLD+' of\nTHIS Script '+color.GREEN+'(4_InFrameStopFreq.py)\n'+color.END\
	+usage_msg(), usage=SUPPRESS,formatter_class=RawTextHelpFormatter)
	
	required_arg_group = parser.add_argument_group(color.ORANGE+color.BOLD+'Required Options'+color.END)
 
	required_arg_group.add_argument('--input_file','-in', action='store', required=True,
	help=color.BOLD+color.GREEN+'Fasta file of Nucleotide sequences enriched \nwith'\
	' Eukaryotic protein coding transcripts'+color.END)

	optional_arg_group = parser.add_argument_group(color.ORANGE+color.BOLD+'Options'+color.END)
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
		
	return args


###########################################################################################
###------------------------------- Script Usage Message --------------------------------###
###########################################################################################

def usage_msg():
	return (color.BOLD+color.RED+'\n\nExample usage:'+color.CYAN+' python 4_InFrameStopFreq.py'\
	' --input_file ../Op_me_Xxma/Op_me_Xxma_WTA_NBU.Renamed.fasta'+color.END)


##########################################################################################
###-------- Storage for LARGE (Annoying) Print Statements for Flagged Options ---------###
##########################################################################################

def return_more_info(args):

	valid_arg = 0

	author = (color.BOLD+color.ORANGE+'\n\n\tQuestions/Comments? Email Xyrus (author) at'\
	' maurerax@gmail.com\n\n'+color.END)

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
			elif args.input_file.endswith('WTA_NBU.Renamed.fasta') != True:
				print (color.BOLD+'\n\nInvalid Fasta File! Only Fasta Files that were processed'\
				' with '+color.GREEN+'3_CountOGsUsearcy.py '+color.END+color.BOLD+'are valid\n\n'\
				'However, to bypass that issue, Fasta Files MUST end with '+color.CYAN+\
				'"WTA_NBU.Renamed.fasta"\n\n'+color.END)
				valid_arg += 1
		else:
			print (color.BOLD+color.RED+'\nError:'+color.END+color.BOLD+' The provided Fasta file '\
			'('+color.DARKCYAN+args.input_file.split('/')[-1]+color.END+color.BOLD+')\ndoes not'\
			' exist or is incorrectly formatted.\n\nDouble-check then try again!\n\n'+color.END) 
			valid_arg += 1

	if os.path.isdir('../../Databases/db_StopFreq') != True:
		print (color.BOLD+color.RED+'\nError:'+color.END+color.BOLD+' Cannot find the '\
		+color.ORANGE+'db_StopFreq Folder!\n\n'+color.END+color.BOLD+'Ensure that this folder '\
		'can be found in the main '+color.ORANGE+'Databases Folder'+color.END+color.BOLD\
		+'\n\nThen try once again\n\n.'+color.END)
		valid_arg += 1

	elif os.path.isfile('../../Databases/db_StopFreq/RepEukProts.udb') != True:
		print (color.BOLD+color.RED+'\nError:'+color.END+color.BOLD+' Cannot find the '\
		'Usearch formatted '+color.ORANGE+'Representative Eukaryotic Protein Database!\n\n'+color.END+color.BOLD+\
		'Ensure that they can be found in the '+color.ORANGE+'db_StopFreq folder'+color.END+\
		color.BOLD+',\nwhich can be found in the main '+color.ORANGE+'Databases Folder'+\
		color.END+color.BOLD+'\n\nThen try once again.\n\n'+color.END)
		valid_arg += 1

	return valid_arg
	
	
###########################################################################################
###--------------------------- Does the Inital Folder Prep -----------------------------###
###########################################################################################

def prep_folders(args):

	Stop_folder = '../'+args.input_file.split('/')[1]+'/StopCodonFreq/'
	
	if os.path.isdir(Stop_folder) != True:
		os.system('mkdir '+Stop_folder)

	if os.path.isdir(Stop_folder+'StopCodonFastas') != True:
		os.system('mkdir '+Stop_folder+'StopCodonFastas')			

	if os.path.isdir(Stop_folder+'SpreadSheets') != True:
		os.system('mkdir '+Stop_folder+'SpreadSheets')			

	return Stop_folder+'StopCodonFastas/'


###########################################################################################
###--------------------- Translates Sequences with Each Stop Codon ---------------------###
###########################################################################################

def prep_translations(args):
	print (color.BOLD+'\nIdentifying ORFs in the Fasta file based on the output of'\
	' 3_CountOGsUsearch.py\n'+color.END)
	
	intsv = [i for i in open(args.input_file.replace('.fasta','_allOGCleanresults.tsv')).readlines() if i != '\n']

	inFasta = [i for i in SeqIO.parse(args.input_file,'fasta')]

	prot_dict = {}

	for i in intsv:
#		print i
		prot_dict.setdefault(i.split('\t')[0],[])
		if int(i.split('\t')[6]) < int(i.split('\t')[7]):
			prot_dict[i.split('\t')[0]].append('F')
			if (int(i.split('\t')[6])) < 5:
				prot_dict[i.split('\t')[0]].append(int(i.split('\t')[6])-1)
			else:
				prot_dict[i.split('\t')[0]].append(int(i.split('\t')[6])-1)
			prot_dict[i.split('\t')[0]].append(int(i.split('\t')[7])+3)
		if int(i.split('\t')[7]) < int(i.split('\t')[6]):
			prot_dict[i.split('\t')[0]].append('RC')
			prot_dict[i.split('\t')[0]].append(int(i.split('\t')[6]))
			if (int(i.split('\t')[7])-4) < 5:
				prot_dict[i.split('\t')[0]].append(int(i.split('\t')[7]))
			else:
				prot_dict[i.split('\t')[0]].append(int(i.split('\t')[7])-4)


	#------------- Prep translation with 'TAA' as the only Stop -------------#
	
	print (color.BOLD+'\n\nTranslating DNA using'+color.RED+' TAA'+color.END\
	+color.BOLD+' as the sole STOP codon\n'+color.END)

	for key, value in prot_dict.items():
		for seq_rec in inFasta:
			if key in seq_rec.description:
				stop_pos = 0
				if prot_dict[key][0] == 'F':
					temp = seq_rec.seq[prot_dict[key][1]:]
					temp_prot = str(temp.translate(table=c_uncinata_table))
					if '*' in temp_prot:		
						stop_pos = (temp_prot.index('*')+1)*3
						prot_dict[key].append(temp[:stop_pos])
					else:
						prot_dict[key].append(seq_rec.seq[prot_dict[key][1]:prot_dict[key][2]])
				if prot_dict[key][0] == 'RC':
					temp = seq_rec.seq[:prot_dict[key][1]].reverse_complement()
					temp_prot = str(temp.translate(table=c_uncinata_table))
					if '*' in temp_prot:		
						stop_pos = (temp_prot.index('*')+1)*3
						prot_dict[key].append(temp[:stop_pos])
					else:
						prot_dict[key].append(seq_rec.seq[prot_dict[key][2]:prot_dict[key][1]].reverse_complement())


	#------------- Prep translation with 'TGA' as the only Stop -------------#
	print (color.BOLD+'\n\nTranslating DNA using'+color.RED+' TGA'+color.END\
	+color.BOLD+' as the sole STOP codon\n'+color.END)	
	
	for key, value in prot_dict.items():
		for seq_rec in inFasta:
			if key in seq_rec.description:
				stop_pos = 0
				if prot_dict[key][0] == 'F':
					temp = seq_rec.seq[prot_dict[key][1]:]
					temp_prot = str(temp.translate(table=6))
					if '*' in temp_prot:		
						stop_pos = (temp_prot.index('*')+1)*3
						prot_dict[key].append(temp[:stop_pos])
					else:
						prot_dict[key].append(seq_rec.seq[prot_dict[key][1]:prot_dict[key][2]])
				if prot_dict[key][0] == 'RC':
					temp = seq_rec.seq[:prot_dict[key][1]].reverse_complement()
					temp_prot = str(temp.translate(table=6))
					if '*' in temp_prot:		
						stop_pos = (temp_prot.index('*')+1)*3
						prot_dict[key].append(temp[:stop_pos])
					else:
						prot_dict[key].append(seq_rec.seq[prot_dict[key][2]:prot_dict[key][1]].reverse_complement())
	
			
	#------------- Prep translation with 'TAG' as the only Stop -------------#
	print (color.BOLD+'\n\nTranslating DNA using'+color.RED+' TAG'+color.END\
	+color.BOLD+' as the sole STOP codon\n'+color.END)
	
	
	for key, value in prot_dict.items():
		for seq_rec in inFasta:
			if key in seq_rec.description:
				stop_pos = 0
				if prot_dict[key][0] == 'F':
					temp = seq_rec.seq[prot_dict[key][1]:]
					temp_prot = str(temp.translate(table=tag_table))
					if '*' in temp_prot:		
						stop_pos = (temp_prot.index('*')+1)*3
						prot_dict[key].append(temp[:stop_pos])
					else:
						prot_dict[key].append(seq_rec.seq[prot_dict[key][1]:prot_dict[key][2]])
				if prot_dict[key][0] == 'RC':
					temp = seq_rec.seq[:prot_dict[key][1]].reverse_complement()
					temp_prot = str(temp.translate(table=tag_table))
					if '*' in temp_prot:		
						stop_pos = (temp_prot.index('*')+1)*3
						prot_dict[key].append(temp[:stop_pos])
					else:
						prot_dict[key].append(seq_rec.seq[prot_dict[key][2]:prot_dict[key][1]].reverse_complement())
					
	#------------ Parsing through data to maintain OG assignments ------------#
	inOGs = intsv
	inOGs = [i.split('\t')[0]+';'+i.split('\t')[1][-10:] for i in inOGs]
	inOGs2 = []
	for i in inOGs:
		if 'no_group' not in i.split(';')[1]:
			inOGs2.append(i)
		else:
			inOGs2.append(i.split(';')[0]+';no_group')
	inOGs2 = list(set(inOGs2))

	#---------------- Write file with 'TAA' is the only Stop ----------------#
						
	with open(args.input_file.split('.fas')[0]+'_taa_ORF.fasta','w+') as w:
		print (color.BOLD+'\n\nWriting FASTA files with ORF and Protein sequences with'+color.RED\
		+' TAA '+color.END+color.BOLD+'as only STOP codon\n'+color.END)
		
		for key, value in prot_dict.items():
			for j in inOGs2:
				if key == j.split(';')[0]:
					if len(prot_dict[key]) < 4:
						pass
					else:
						w.write('>'+key+'_'+j.split(';')[1]+'\n'+str(value[-3]).upper()+'\n')
						
	with open(args.input_file.split('.fas')[0]+'_taa_ORF.aa.fasta','w+') as w:
		for key, value in prot_dict.items():
			for j in inOGs2:
				if key == j.split(';')[0]:
					if len(prot_dict[key]) < 4:
						pass
					else:
						w.write('>'+key+'_'+j.split(';')[1]+'\n'+str(Seq(str(value[-3])).translate(table=c_uncinata_table)).upper()+'\n')					

	#---------------- Write file with 'TGA' is the only Stop ----------------#
	
	with open(args.input_file.split('.fas')[0]+'_tga_ORF.fasta','w+') as w:
		print (color.BOLD+'\n\nWriting FASTA files with ORF and Protein sequences with'+color.RED\
		+' TGA '+color.END+color.BOLD+'as only STOP codon\n'+color.END)
	
		for key, value in prot_dict.items():
			for j in inOGs2:
				if key == j.split(';')[0]:
					if len(prot_dict[key]) < 4:
						pass
					else:
						w.write('>'+key+'_'+j.split(';')[1]+'\n'+str(value[-2]).upper()+'\n')
					
	with open(args.input_file.split('.fas')[0]+'_tga_ORF.aa.fasta','w+') as w:
		for key, value in prot_dict.items():
			for j in inOGs2:
				if key == j.split(';')[0]:
					if len(prot_dict[key]) < 4:
						pass
					else:
						w.write('>'+key+'_'+j.split(';')[1]+'\n'+str(Seq(str(value[-2])).translate(table=6)).upper()+'\n')
						
	#---------------- Write file with 'TAG' is the only Stop ----------------#
	
	with open(args.input_file.split('.fas')[0]+'_tag_ORF.fasta','w+') as w:
		print (color.BOLD+'\n\nWriting FASTA files with ORF and Protein sequences with'+color.RED\
		+' TAG '+color.END+color.BOLD+'as only STOP codon\n'+color.END)
		
		for key, value in prot_dict.items():
			for j in inOGs2:
				if key == j.split(';')[0]:
					if len(prot_dict[key]) < 4:
						pass
					else:
						w.write('>'+key+'_'+j.split(';')[1]+'\n'+str(value[-1]).upper()+'\n')
	
	with open(args.input_file.split('.fas')[0]+'_tag_ORF.aa.fasta','w+') as w:
		for key, value in prot_dict.items():
			for j in inOGs2:
				if key == j.split(';')[0]:
					if len(prot_dict[key]) < 4:
						pass
					else:
						w.write('>'+key+'_'+j.split(';')[1]+'\n'+str(Seq(str(value[-1])).translate(table=tag_table)).upper()+'\n')				


###########################################################################################		
###---------- USearches the Translations Against a SMALL Euk Protein Database ----------###
###########################################################################################

def ublast_ProtDB(args, usearch_path):
	os.system(usearch_path+' -ublast '+args.input_file.split('.fas')[0]+'_tag_ORF.aa.fasta '\
	'-db ../../Databases/db_StopFreq/RepEukProts.udb -evalue 1e-5 -maxaccepts 1 -blast6out '\
	+args.input_file.split('.fas')[0]+'_tag_ORF.RepEukProts.tsv')

	os.system(usearch_path+' -ublast '+args.input_file.split('.fas')[0]+'_tga_ORF.aa.fasta '\
	'-db ../../Databases/db_StopFreq/RepEukProts.udb -evalue 1e-5 -maxaccepts 1 -blast6out '\
	+args.input_file.split('.fas')[0]+'_tga_ORF.RepEukProts.tsv')

	os.system(usearch_path+' -ublast '+args.input_file.split('.fas')[0]+'_taa_ORF.aa.fasta '\
	'-db ../../Databases/db_StopFreq/RepEukProts.udb -evalue 1e-5 -maxaccepts 1 -blast6out '\
	+args.input_file.split('.fas')[0]+'_taa_ORF.RepEukProts.tsv')


###########################################################################################		
###-------------------- Manages the search for In-Frame Stop Codons --------------------###
###########################################################################################


def hunt_for_stops(args):
		
	#------------------------ Open Fasta Files ------------------------#
	try: 
		TAGinFasta = [i for i in SeqIO.parse(args.input_file.split('.fas')[0]+'_tag_ORF.fasta','fasta') if str(i.seq).endswith('TAG')]
		print (color.BOLD+'\n\nGathering Sequence information from FASTA and TSV files\n'+color.END)
		
	except:
		print (color.BOLD+color.RED+'\n\nMissing Necessary Inputs: Open Script for Usage'\
		' Information\n\n'+color.END)
		sys.exit()
		
	TGAinFasta = [i for i in SeqIO.parse(args.input_file.split('.fas')[0]+'_tga_ORF.fasta','fasta') if str(i.seq).endswith('TGA')]

	TAAinFasta = [i for i in SeqIO.parse(args.input_file.split('.fas')[0]+'_taa_ORF.fasta','fasta') if str(i.seq).endswith('TAA')]
	
	## This section originally ONLY considered sequences WITH OG assignments:
	## TAAinFasta = [i for i in TAAinFasta if 'no_group' not in i.description and str(i.seq).endswith('TAA')]
	## This has been taken out for now
	
	#----------------------- Open BLAST Reports -----------------------#

	TAGinTSV = [i for i in open(args.input_file.split('.fas')[0]+'_tag_ORF.RepEukProts.tsv').read().split('\n') if i != '']

	TGAinTSV = [i for i in open(args.input_file.split('.fas')[0]+'_tga_ORF.RepEukProts.tsv').read().split('\n') if i != '']

	TAAinTSV = [i for i in open(args.input_file.split('.fas')[0]+'_taa_ORF.RepEukProts.tsv').read().split('\n') if i != '']

## This section originally ONLY considered sequences WITH OG assignments:
	## TAAinTSV = i for i in TAAinTSV if i != ''and 'no_group' not in i.split('\t')[0]]
	## This has been taken out for now

	
	#------------ Set-up Genetic Code Specific Dictionaries ------------#
	
	tag_dict = {}
	for i in TAGinTSV:
		tag_dict.setdefault(i.split('\t')[0].replace('_TAG',''),[]).append(int(i.split('\t')[-6]))
		tag_dict.setdefault(i.split('\t')[0].replace('_TAG',''),[]).append(int(i.split('\t')[-5]))

	tga_dict = {}
	for i in TGAinTSV:
		tga_dict.setdefault(i.split('\t')[0].replace('_Ciliate',''),[]).append(int(i.split('\t')[-6]))
		tga_dict.setdefault(i.split('\t')[0].replace('_Ciliate',''),[]).append(int(i.split('\t')[-5]))

	taa_dict = {}
	for i in TAAinTSV:
		taa_dict.setdefault(i.split('\t')[0].replace('_Chilo',''),[]).append(int(i.split('\t')[-6]))
		taa_dict.setdefault(i.split('\t')[0].replace('_Chilo',''),[]).append(int(i.split('\t')[-5]))

	#-------------- Preparing In-Frame Stop Codon Counts --------------#

# All the data when TGA is the sole stop codon
	tga_codons = 0
	tga_data_tag = 0
	tga_data_tga = 0
	tga_data_taa = 0
	tga_seq_count = 0 

# All the data when TAG is the sole stop codon	
	tag_codons = 0
	tag_data_tag = 0
	tag_data_tga = 0
	tag_data_taa = 0
	tag_seq_count = 0
	
# All the data when TAA is the sole stop codon	
	taa_codons = 0
	taa_data_tag = 0
	taa_data_tga = 0
	taa_data_taa = 0
	taa_seq_count = 0
	
# All the data for each stop codon combined
	tga_inframe = 0
	tag_inframe = 0
	taa_inframe = 0
	total_codons = 0
	total_seq_counts = len(open(args.input_file).read().split('>'))-1


	#-------- Gathering In-frame Stop Codon Density Information --------#
	
### Collect in-frame stop information for "TAA" and "TAG" when TGA is the ONLY stop	
	print (color.BOLD+'\nCollecting in-frame stop codon information when'+color.RED\
	+' TGA'+color.END+color.BOLD+' is the only STOP\n'+color.END)
	
	for i in TGAinFasta:
		try:
			if tga_dict[i.description][0] == 1:
				for n in range((tga_dict[i.description][0]-1),((tga_dict[i.description][1])*3)-3,3):
					if str(i.seq).upper()[n:n+3] == 'TAG':
						tga_data_tag += 1
						tag_inframe += 1
					if str(i.seq).upper()[n:n+3].upper() == 'TAA':
						tga_data_taa += 1
						taa_inframe += 1
					if str(i.seq).upper()[n:n+3].upper() == 'TGA':
						tga_data_tga += 1
						tga_inframe += 1
					tga_codons += 1
					total_codons += 1
				tga_seq_count += 1 
					
			else:
				for n in range(((tga_dict[i.description][0]-1)*3),((tga_dict[i.description][1])*3)-3,3):
					if str(i.seq).upper()[n:n+3] == 'TAG':
						tga_data_tag += 1
						tag_inframe += 1
					if str(i.seq).upper()[n:n+3].upper() == 'TAA':
						tga_data_taa += 1
						taa_inframe += 1
					if str(i.seq).upper()[n:n+3].upper() == 'TGA':
						tga_data_tga += 1
						tga_inframe += 1
					tga_codons += 1 
					total_codons += 1
				tga_seq_count += 1
		except:
			pass		
			
### Collect in-frame stop information for "TAA" and "TGA" when TAG is the ONLY stop	
	print (color.BOLD+'\nCollecting in-frame stop codon information when'+color.RED\
	+' TAG'+color.END+color.BOLD+' is the only STOP\n'+color.END)	

	for i in TAGinFasta:
		try:
			if tag_dict[i.description][0] == 1:
				for n in range((tag_dict[i.description][0]-1),((tag_dict[i.description][1])*3)-3,3):
					if str(i.seq).upper()[n:n+3] == 'TAG':
						tag_data_tag += 1
						tag_inframe += 1
					if str(i.seq).upper()[n:n+3].upper() == 'TAA':	
						tag_data_taa += 1
						taa_inframe += 1
					if str(i.seq).upper()[n:n+3].upper() == 'TGA':	
						tag_data_tga += 1
						tga_inframe += 1
					tag_codons += 1
					total_codons += 1 
				tag_seq_count += 1
				
			else:
				for n in range(((tag_dict[i.description][0]-1)*3),(tag_dict[i.description][1]*3)-3,3):
					if str(i.seq).upper()[n:n+3] == 'TAG':
						tag_data_tag += 1
						tag_inframe += 1
					if str(i.seq).upper()[n:n+3].upper() == 'TAA':	
						tag_data_taa += 1
						taa_inframe += 1
					if str(i.seq).upper()[n:n+3].upper() == 'TGA':	
						tag_data_tga += 1
						tga_inframe += 1
					tag_codons += 1 
					total_codons += 1
				tag_seq_count += 1
		except:
			pass


### Collect in-frame stop information for "TGA" and "TAG" when TAA is the ONLY stop	
	print (color.BOLD+'\nCollecting in-frame stop codon information when'+color.RED\
	+' TAA'+color.END+color.BOLD+' is the only STOP\n'+color.END)
		
	for i in TAAinFasta:
		try:
			if taa_dict[i.description][0] == 1:
				for n in range((taa_dict[i.description][0]-1),((taa_dict[i.description][1])*3)-3,3):
					if str(i.seq).upper()[n:n+3] == 'TAG':
						taa_data_tag += 1
						tag_inframe += 1
					if str(i.seq).upper()[n:n+3].upper() == 'TAA':	
						taa_data_taa += 1
						taa_inframe += 1
					if str(i.seq).upper()[n:n+3].upper() == 'TGA':	
						taa_data_tga += 1
						tga_inframe += 1
					taa_codons += 1 
					total_codons += 1
				taa_seq_count += 1
				
			else:
				for n in range(((taa_dict[i.description][0]-1)*3),(taa_dict[i.description][1]*3)-3,3):
					if str(i.seq).upper()[n:n+3] == 'TAG':
						taa_data_tag += 1
						tag_inframe += 1
					if str(i.seq).upper()[n:n+3].upper() == 'TAA':	
						taa_data_taa += 1
						taa_inframe += 1
					if str(i.seq).upper()[n:n+3].upper() == 'TGA':	
						taa_data_tga += 1
						tga_inframe += 1
					tag_codons += 1 
					total_codons += 1
				taa_seq_count += 1
		except:
			pass

	#-------------- Writing Data Out and Print Statement --------------#	

	with open(args.input_file.split('.fas')[0]+'_StopCodonStats.tsv','w+') as w:
		w.write('Stop Codon\tNumber of Seqs Analyzed\tIn-frame TAG\tIn-frame TGA\tIn-frame TAA\tTotal Codons\tIn-frame TAG density\tIn-frame TGA density\tIn-frame TAA density\n')
		if tga_codons != 0:
			w.write('TGA\t'+str(tga_seq_count)+'\t'+str(tga_data_tag)+'\t'+str(tga_data_tga)+'\t'+str(tga_data_taa)+'\t'+str(tga_codons)\
			+'\t'+"%.2f" % ((float(tga_data_tag)*1000)/float(tga_codons))+'\t'+"%.2f" % ((float(tga_data_tga)*1000)/float(tga_codons))+'\t'\
			+"%.2f" % ((float(tga_data_taa)*1000)/float(tga_codons))+'\n')
		else:
			w.write('TGA\t0\t0\t0\t0\t0\t0\t0\t0\n')
	
		if tag_codons != 0:
			w.write('TAG\t'+str(tag_seq_count)+'\t'+str(tag_data_tag)+'\t'+str(tag_data_tga)+'\t'+str(tag_data_taa)+'\t'+str(tag_codons)\
			+'\t'+"%.2f" % ((float(tag_data_tag)*1000)/float(tag_codons))+'\t'+"%.2f" % ((float(tag_data_tga)*1000)/float(tag_codons))+'\t'\
			+"%.2f" % ((float(tag_data_taa)*1000)/float(tag_codons))+'\n')
		else:
			w.write('TAG\t0\t0\t0\t0\t0\t0\t0\t0\n')
		if taa_codons != 0:
			w.write('TAA\t'+str(taa_seq_count)+'\t'+str(taa_data_tag)+'\t'+str(taa_data_tga)+'\t'+str(taa_data_taa)+'\t'+str(taa_codons)\
			+'\t'+"%.2f" % ((float(taa_data_tag)*1000)/float(taa_codons))+'\t'+"%.2f" % ((float(taa_data_tga)*1000)/float(taa_codons))+'\t'\
			+"%.2f" % ((float(taa_data_taa)*1000)/float(taa_codons))+'\n')
		else:
			w.write('TAA\t0\t0\t0\t0\t0\t0\t0\t0\n')
	
		w.write('\n \n')
		w.write('Summary\t'+str(tga_seq_count+tag_seq_count+taa_seq_count)+'\t'+str(tag_inframe)+'\t'+str(tga_inframe)+'\t'+str(taa_inframe)\
		+'\t'+str(total_codons)+'\t'+"%.2f" % ((float(tag_inframe)*1000)/float(total_codons))+'\t'+"%.2f" % ((float(tga_inframe)*1000)/float(total_codons))\
		+'\t'+"%.2f" % ((float(taa_inframe)*1000)/float(total_codons))+'\n')
		w.write('\nTotal Seqs in Fasta\t'+str(total_seq_counts))
	
#	print color.BOLD + color.BLUE + '\nSummary\t'+str(tag_inframe)+'\t'+str(tga_inframe)+'\t'+str(taa_inframe)+'\t'+str(total_codons)+'\t'+"%.2f" % ((float(tag_inframe)*1000)/float(total_codons))+'\t'\
#		+"%.2f" % ((float(tga_inframe)*1000)/float(total_codons))+'\t'+"%.2f" % ((float(taa_inframe)*1000)/float(total_codons))+'\n\n'\
#		+ str(tag_seq_count) + '\t' + str(tga_seq_count) + '\t' + str(taa_seq_count) + color.END


##########################################################################################
###--------------------- Cleans up the Folder and Moves Final Files -------------------###
##########################################################################################

def clean_up(args):
	if os.path.isdir('../'+args.input_file.split('/')[1]+'/StopCodonFreq') != True:
		os.system('mkdir ../'+args.input_file.split('/')[1]+'/StopCodonFreq/')
	else:
		pass
			
	os.system('mkdir ../'+args.input_file.split('/')[1]+'/StopCodonFreq/StopCodonFastas/')
	os.system('mkdir ../'+args.input_file.split('/')[1]+'/StopCodonFreq/SpreadSheets/')
	os.system('mv '+args.input_file.split('.fas')[0]+'_t*_ORF.*fasta ../'+args.input_file.split('/')[1]+'/StopCodonFreq/StopCodonFastas/')
	os.system('mv '+args.input_file.split('.fas')[0]+'_t*Prots.tsv ../'+args.input_file.split('/')[1]+'/StopCodonFreq/SpreadSheets/')


###########################################################################################
###-------------------------------- Next Script Message --------------------------------###
###########################################################################################

def next_script(args):

	home_folder = args.input_file.split('/')[1]+'/'

	print (color.BOLD+'\nLook for '+color.DARKCYAN+args.input_file.split('/')[-1]\
	.replace('.fasta','_StopCodonStats.tsv')+color.END+color.BOLD+' in the '+home_folder\
	+' Folder\n\n' + color.END)
	
	print (color.BOLD+'Next Script is: '+color.GREEN+'5_GCodeTranslate.py\n\n'+ color.END)


##########################################################################################
###--------------- Checks Command Line Arguments and Calls on Functions ---------------###
##########################################################################################

def main():

	usearch_path = check_usearch_path()

	args = check_args()
	
	prep_translations(args)
	
	ublast_ProtDB(args, usearch_path)

	hunt_for_stops(args)
	
	clean_up(args)
	
	next_script(args)
 
main()


#----------------------------------------- NOTES -----------------------------------------#
#
# This script is designed to HELP you make an informed decision about the genetic code being
# used by your particular organism. Be aware that it will be limited by the quality of the
# data given to it! 
#
# You will need:
#
# USearch, BioPython, AND the output from '3_CountOGSUsearch.py'
#
# If you are not using the Author's database, update your database name(s) in lines: 345-360
#
# katzlab$ python StopFrequency.py YourFastaFile.fasta
#
#
#------------------------------- Interpretation of Results -------------------------------#
#
# FORMATTED BELOW WITH TEXTWRANGLER... 
#
# Example output using CILIATE (TGA) genetic Code (NOTE THE In-Frame Densities):
#
# Stop Codon   Number_of_Seqs_Analyzed   In-frame TAG   In-frame TGA   In-frame TAA   Total Codons   In-frame TAG density   In-frame TGA density   In-frame TAA density
# TGA                   341                   14              0			    	  22            113156              1.2                     0                     0.92
# TAG                   424                    0              0			    	  34            140085               0                      0                     0.78
# TAA                   205                   14              0			    	   0             16714              0.84                    0                      0
# Summary               970                   28              0			    	  56            269955              2.04                    0                     1.7
#
# VALUES in summary line (OR SUM of Density) that are > 1.5 likely indicate that the STOP 
# codon has been reassigned... in the case above, TAG and TAA look like they have been
# reassigned.
#
#
# Example output using UNIVERSAL genetic Code (NOTE THE In-Frame Densities):
#
# Stop Codon   Number_of_Seqs_Analyzed   In-frame TAG   In-frame TGA   In-frame TAA   Total Codons   In-frame TAG density   In-frame TGA density   In-frame TAA density
# TGA                   341                    1              0			    	   2            113156              0.2                     0                     0.05
# TAG                   424                    0              2			    	   4            140085               0                      0                     0.08
# TAA                   205                    1              0			    	   0             16714              0.04                    0                      0
# Summary               970                    2              2			    	   6            269955              0.15                    0                     0.06
#
# VALUES in summary line (OR SUM of Density) that are > 0.5 likely indicate that the STOP 
# codon still acts as STOP... in the case above, TAG, TGA and TAA look like they still behave
# as a stop codon.
#
# THIS IS A ROUGH GUIDE FOR INTERPRETING THE RESULTS!!!! BE VERY VERY WARY! NUMBER OF TOTAL
# SEQUENCES AND TOTAL CODONS OBSERVED ARE IMPORTANT (TOO FEW AND ANY INTERPRETATION IS DEVOID
# OF ANY MEANING).