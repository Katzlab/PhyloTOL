#!/usr/bin/env python3.5

##__Updated__: 16_10_2017
##__Author__: Xyrus Maurer-Alcala; maurerax@gmail.com
##__Usage__: python 3_CountOGsUsearch.py --help


##########################################################################################
## This script is intended to classify the STRONGLY Eukaryotic and UNDETERMINED/UNKNOWN ##
## contigs into different OGs (e.g. orthologous gene-families)							##
##																						##
## For more info about the OGs, check out: OrthoMCL.org									##
##																						##
## Prior to running this script, ensure the following:									##
##																						##
## 1. You have assembled your transcriptome and COPIED the 'assembly' file 				##
##    (contigs.fasta, or scaffolds.fasta) to the PostAssembly Folder					##
## 2. Removed small sequences (usually sequences < 300bp) with ContigFilterPlusStats.py	##
## 3. Removed SSU/LSU sequences from your Fasta File									## 
##																						##
## 			E-mail Xyrus (author) for help if needed: maurerax@gmail.com				##
##																						##
##								Next Script(s) to Run: 									##
##							 	 4_StopFrequency.py										##
##																						##
##########################################################################################

import argparse, os, sys
from argparse import RawTextHelpFormatter,SUPPRESS
from distutils import spawn
from Bio import SeqIO


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
	color.BOLD + '\n\nThis script will categorize Contigs into'+color.ORANGE+' "Homologous" '\
	+color.END+color.BOLD+'Gene Families (OGs)\nbased on '+color.RED+'OrthoMCL'+color.END\
	+color.BOLD+"'s Gene Family Grouping\n\n\nNotes on this script and "+color.GREEN+\
	'OrthoMCL Families'+color.END+color.BOLD+' can be found\nat the bottom of '+color.GREEN\
	+'THIS script (3_CountOGsUsearch.py)\n'+color.END+usage_msg(), usage=SUPPRESS,
	formatter_class=RawTextHelpFormatter)
	
	required_arg_group = parser.add_argument_group(color.ORANGE+color.BOLD+'Required Options'+color.END)
 
	required_arg_group.add_argument('--input_file','-in', action='store',
	help=color.BOLD+color.GREEN+'Fasta file of Nucleotide sequences enriched \nwith'\
	' Eukaryotic protein coding transcripts'+color.END)

	optional_arg_group = parser.add_argument_group(color.ORANGE+color.BOLD+'Options'+color.END)
	optional_arg_group.add_argument('--threads','-t', default='2',
	help=color.BOLD+color.GREEN+' Number of threads to use for BLAST\n (default = 2)\n'+color.END)
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
	return (color.BOLD+color.RED+'\n\nExample usage:'+color.CYAN+' python 3_CountOGsUsearch.py'\
	' --input_file ../Op_me_Xxma/Op_me_Xxma_WTA_NBU.fasta'+color.END)


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
			elif args.input_file.endswith('WTA_NBU.fasta') != True:
				print (color.BOLD+'\n\nInvalid Fasta File! Only Fasta Files that were processed'\
				' with '+color.GREEN+'2b_remove_Bact.py '+color.END+color.BOLD+'are valid\n\n'\
				'However, to bypass that issue, Fasta Files MUST end with '+color.CYAN+\
				'"WTA_NBU.fasta"\n\n'+color.END)
				valid_arg += 1
		else:
			print (color.BOLD+color.RED+'\nError:'+color.END+color.BOLD+' The provided Fasta file '\
			'('+color.DARKCYAN+args.input_file.split('/')[-1]+color.END+color.BOLD+')\ndoes not'\
			' exist or is incorrectly formatted.\n\nDouble-check then try again!\n\n'+color.END) 
			valid_arg += 1

	if os.path.isdir('../../Databases/db_OG') != True:
		print (color.BOLD+color.RED+'\nError:'+color.END+color.BOLD+' Cannot find the '\
		+color.ORANGE+'db_OG Folder!\n\n'+color.END+color.BOLD+'Ensure that this folder '\
		'can be found in the main '+color.ORANGE+'Databases Folder'+color.END+color.BOLD\
		+'\n\nThen try once again\n\n.'+color.END)
		valid_arg += 1

	elif os.path.isfile('../../Databases/db_OG/OGSout0.udb') != True:
		print (color.BOLD+color.RED+'\nError:'+color.END+color.BOLD+' Cannot find the '\
		'Usearch formatted '+color.ORANGE+'Gene Family databases!\n\n'+color.END+color.BOLD+\
		'Ensure that they can be found in the '+color.ORANGE+'db_OG folder'+color.END+\
		color.BOLD+',\nwhich can be found in the main '+color.ORANGE+'Databases Folder'+\
		color.END+color.BOLD+'\n\nThen try once again.\n\n'+color.END)
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

		
###########################################################################################
###--------------------- Runs Usearch on Split OrthoMCL Databases ----------------------###
###########################################################################################

def OG_ublast(args, usearch_path):
	
	print (color.BOLD+'\nStarting to "BLAST" against OG databases'+color.END)
		
	OG_folder = '../'+args.input_file.split('/')[1]+'/UsearchOG/OG_CompData/'
	UBlast_output = args.input_file.split('/')[-1].split('.fas')[0]+'_OGSout'
	
	for db in os.listdir(os.curdir + './../Databases/db_OG'):
		if db.startswith('OGS') and db.endswith('.udb'):
			dbnum = db.split('out')[-1].split('.')[0]
	
			print (color.BOLD+'\n\n"BLAST"-ing against OG database: '+color.DARKCYAN+\
			'OGSout'+dbnum+'.udb'+color.END+'\n\n')			
			
			# Update the path to USearch for your computer's architecture (if in the PATH
			# you can just call 'usearch' without the full path).		
			
			OG_usearch_cmd = usearch_path+' -ublast '+args.input_file+' -maxaccepts 1'\
			' -db ../../Databases/db_OG/OGSout'+dbnum+'.udb -evalue 1e-5 -target_cov 0.35 -blast6out '\
			+OG_folder+UBlast_output+dbnum		
			
			os.system(OG_usearch_cmd)
			
			
	os.system('cat '+OG_folder+'*OGSout* > '+OG_folder+args.input_file.split('/')[-1].split('.fas')[0]\
	+'_allOGresults')
	
	os.system('rm '+OG_folder+'*OGSout*')		


###########################################################################################
###--------------- Keeps the Single BEST Hit (HSP-score) Per Transcript ----------------###
###########################################################################################

def keep_best(args):

	print (color.BOLD+color.PURPLE+'\n\nProcessing OG-database results to keep only the'\
	' BEST match for each transcript\n\n'+color.END)
	
	OG_folder = '../'+args.input_file.split('/')[1]+'/UsearchOG/'
		
	inTSV = [i for i in open(OG_folder+'OG_CompData/'+args.input_file.split('/')[-1].split('.fas')[0]\
	+'_allOGresults').readlines()]
	
	inTSV.sort(key = lambda x: -float(x.split('\t')[-1]))
	
	keep = []
	for i in inTSV:
		if any(i.split('\t')[0] in j for j in keep) != True:
			keep.append(i)

	updated_lines = list(set([line.split('\t')[0]+'_'+'_'.join(line.split('\t')[1].split('_')[-2:])+'\t'+'\t'.join(line.split('\t')[1:]) for line in keep]))
		
	with open(args.input_file.replace('.fasta','.Renamed_allOGCleanresults.tsv'), 'w+') as w:
		for i in updated_lines:
			w.write(i)
	

###########################################################################################
###-------- Copies and Updates Names of Transcripts With OG Hits to New Fasta ----------###
###########################################################################################

def update_fasta(args):

	print (color.BOLD+color.PURPLE+'Updating Fasta File Sequence Names with their'\
	' BEST OG hits\n\n'+color.END)

	Renamed_TSV = args.input_file.replace('.fasta','.Renamed_allOGCleanresults.tsv')

	keep = [i for i in open(Renamed_TSV).readlines() if i != '\n']

	keep_dict = {line.split('\t')[0].split('_OG5')[0]:line.split('\t')[0].split('_OG5')[0]+'_OG5_'+line.split('\t')[1].split('_')[-1] for line in keep if 'OG5' in line.split('\t')[1]}
	
	inFasta = [i for i in SeqIO.parse(args.input_file,'fasta')]
	
	updated_seq_name = ['>'+keep_dict[i.description]+'\n'+str(i.seq)+'\n' for i in inFasta if i.description in keep_dict.keys()]

	seqs_without_OG = ['>'+i.description+'\n'+str(i.seq)+'\n' for i in inFasta if i.description not in keep_dict.keys()]
	
	with open(args.input_file.replace('.fasta','.Renamed.fasta'),'w+') as w:
		for i in updated_seq_name:
			w.write(i)		

	with open(args.input_file.replace('.fasta','.LackOG.fasta'),'w+') as x:
		for i in seqs_without_OG:
			x.write(i)		


###########################################################################################
###-------------------- Updates Log With OG Assignment Information ---------------------###
###########################################################################################
	
def update_log(args):

	if os.path.isdir('../PostAssembly_Logs/') != True:
		os.system('mkdir ../PostAssembly_Logs/')
	else:
		pass

	home_folder = '../'+args.input_file.split('/')[1]+'/'
	
	Renamed_TSV = home_folder+args.input_file.split('/')[-1].replace('.fasta','.Renamed_allOGCleanresults.tsv')

	keep = [line for line in open(Renamed_TSV).readlines()]
	all_ogs = [line.split('\t')[1].split('_')[-1] for line in keep if 'OG5' in line.split('\t')[1]]

	total_with_ogs = str(len(all_ogs))
	unique_ogs = str(len(set(all_ogs)))
	
	print (color.BOLD +'There are '+color.BLUE +total_with_ogs+' Contigs'+color.END\
	+color.BOLD+' that hit '+color.DARKCYAN+unique_ogs+' Unique OGs\n'+color.END)

	
	for Logname in os.listdir(os.curdir+'./PostAssembly_Logs/'):
		if Logname.startswith(args.input_file.split('/')[2].split('_WTA')[0]) and Logname.endswith('Log.txt'):
			with open('../PostAssembly_Logs/'+Logname,'a') as LogFile:
				LogFile.write('Contigs With OG\t'+total_with_ogs+'\tn/a\tn/a\n')
				LogFile.write('Unique OGs\t'+unique_ogs+'\tn/a\tn/a\n')


##########################################################################################
###--------------------- Cleans up the Folder and Moves Final Files -------------------###
##########################################################################################

def clean_up(args):

	OG_folder = '../'+args.input_file.split('/')[1]+'/UsearchOG/'
	
	os.system('rm '+args.input_file)
	
	os.system('cp '+args.input_file.replace('.fasta','.Renamed.fasta')+' '+OG_folder)

	os.system('cp '+args.input_file.replace('.fasta','.Renamed_allOGCleanresults.tsv')+' '+OG_folder)


###########################################################################################
###-------------------------------- Next Script Message --------------------------------###
###########################################################################################

def next_script(args):

	home_folder = '../'+args.input_file.split('/')[1]+'/'

	print (color.BOLD+'\nLook for '+color.DARKCYAN+args.input_file.split('/')[-1]\
	.replace('.fasta','WTA_NBU.fasta')+color.END+color.BOLD+' in the '+home_folder\
	+' Folder\n\n' + color.END)
	
	print (color.BOLD+'Next Script is: '+color.GREEN+'4_InFrameStopFreq.py\n\n'+ color.END)


##########################################################################################
###--------------- Checks Command Line Arguments and Calls on Functions ---------------###
##########################################################################################
				
def main():

	usearch_path = check_usearch_path()

	args = check_args()

	prep_folders(args)

	OG_ublast(args, usearch_path)

	keep_best(args)

	update_fasta(args)

	update_log(args)		

	clean_up(args)

	next_script(args)

main()
