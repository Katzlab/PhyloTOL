#!/usr/bin/env python3.5

##__Updated__: 19_09_2017
##__Author__: Xyrus Maurer-Alcala; maurerax@gmail.com
##__Usage__: python 1p_RenameProts.py --help


from Bio import SeqIO
from Bio.SeqUtils import GC
import argparse, os, sys, time
from argparse import RawTextHelpFormatter,SUPPRESS


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
###--------------------- Parses and Checks Command-Line Arguments ----------------------###
###########################################################################################

def check_args():

	parser = argparse.ArgumentParser(description=
	color.BOLD + '\n\nThis script is intended to reformat the naming for '+color.PURPLE+\
	'Protein Sequences\n'+color.END+color.BOLD+'from a provided Genbank formatted file.'\
	+usage_msg(), usage=SUPPRESS, formatter_class=RawTextHelpFormatter)
	
	required_arg_group = parser.add_argument_group(color.ORANGE+color.BOLD+'Required Options'+color.END)
 
	required_arg_group.add_argument('--input_file','-in', action='store',
	help=color.BOLD+color.GREEN+' Fasta file with Protein Sequences\n'+color.END)

	optional_arg_group = parser.add_argument_group(color.ORANGE+color.BOLD+'Options'+color.END)

	optional_arg_group.add_argument('--source','-s', action='store', default='GenBank',
	help=color.BOLD+color.GREEN+' Data Source of CDSs (default = "GenBank")\n'+color.END)

	optional_arg_group.add_argument('--list_source','-lsrc', action='store_true',
	help=color.BOLD+color.GREEN+' Lists supported data sources\n'+color.END)

	optional_arg_group.add_argument('-author', action='store_true',
	help=color.BOLD+color.GREEN+' Prints author contact information\n'+color.END)


	if len(sys.argv[1:]) == 0:
		print (parser.description)
		print ('\n')
		sys.exit()

	args = parser.parse_args()
	
	more_info = return_more_info(args)
	if more_info != None:
		print (parser.description)
		print (more_info)
		sys.exit()
	
	args.folder = args.input_file.split('.fasta')[0]
		
	return args
	
	
###########################################################################################
###------------------------------- Script Usage Message --------------------------------###
###########################################################################################

def usage_msg():
	return (color.BOLD+color.RED+'\n\nExample usage:'+color.CYAN+' python 1p_RenameProts.py'\
	' --input_file ../Plasmopara_viticola_Proteins.fasta --source GenBank'+color.END)


##########################################################################################
###-------- Storage for LARGE (Annoying) Print Statements for Flagged Options ---------###
##########################################################################################

def return_more_info(args):

	acceptable_sources = ['in-house', 'in-lab', 'GenBank', 'gb', 'NCBI']

	author = (color.BOLD+color.ORANGE+'\n\n\tQuestions/Comments? Email Xyrus (author) at'\
	' maurerax@gmail.com\n\n'+color.END)

	if args.author == True:
		return author

	if args.list_source == True:
		print (color.BOLD+color.RED+'\nThese are the currently supported data sources.\n'+color.END)
		print (color.BOLD+color.ORANGE+'\n'.join(acceptable_sources)+'\n\n'+color.END)
		sys.exit()

	if args.source.lower() not in [i.lower() for i in acceptable_sources]:
		print (color.BOLD+color.RED+'\nUnsupported source was provided.\n\nEnsure that '\
		'you are providing a valid data source (see below).\n'+color.END)
		print (color.BOLD+color.ORANGE+'\n'.join(acceptable_sources)+'\n'+color.END)
		sys.exit()
	
	if args.input_file != None:
		if args.input_file.split('/')[-1] not in os.listdir('/'.join(args.input_file.split('/')[:-1])):
			print (color.BOLD+color.RED+'\nError:'+color.END+color.BOLD+' The provided Fasta file '\
			'('+color.DARKCYAN+args.input_file.split('/')[-1]+color.END+color.BOLD+')\ndoes not'\
			' exist or is incorrectly formatted.\n\nDouble-check then try again!\n\n'+color.END) 
			sys.exit()
	

###########################################################################################
###--------------------------- Does the Inital Folder Prep -----------------------------###
###########################################################################################

def prep_folders(args):

	if os.path.isdir(args.folder) != True:
		os.system('mkdir '+args.folder)
	if os.path.isdir(args.folder+'/Original') != True:
		os.system('mkdir '+args.folder+'/Original')		
		os.system('cp '+args.input_file+' '+args.folder+'/Original/'+args.input_file.split('/')[-1]\
		.split('.fas')[0]+'.Original.fasta')

###########################################################################################
###------------- Renames Protein-Coding CDS Sequences to Standard Format ---------------###
###########################################################################################

def renamed_GenomeCDS(args):
		
	print (color.BOLD+'\n\nPrepping to rename '+color.GREEN+args.input_file.split('/')[-1]+\
	color.END+color.BOLD+"'s Protein sequences"+color.END)
	inFasta = sorted((i for i in SeqIO.parse(args.input_file,'fasta')),key=lambda seq_rec: -len(seq_rec.seq))

	renamed_seqs = []
	seq_code_dict = {}

	count = 1
	if '-' not in args.source:
		for seq_rec in inFasta:
			seq_code_dict.setdefault(seq_rec.description,[]).append(seq_rec.description.split('_')[-1].split('.')[0]+\
			'_Contig_'+str(count)+'_Len'+str(len(seq_rec.seq)*3))
			seq_code_dict[seq_rec.description].append(str(seq_rec.seq).upper())
			renamed_seqs.append('>'+seq_rec.description.split('_')[-1].split('.')[0]+'_Contig_'\
			+str(count)+'_Len'+str(len(seq_rec.seq))+'\n'+str(seq_rec.seq).upper())
			count += 1
	else:
		for seq_rec in inFasta:
			seq_code_dict.setdefault(seq_rec.description,[]).append('LKH'+fasta_file.split('LKH')[-1].split('_')[0]+'_'+str(count)+'_Len_'+str(len(seq_rec.seq)))
			seq_code_dict[seq_rec.description].append(str(seq_rec.seq).upper())
			renamed_seqs.append('>LKH'+fasta_file.split('LKH')[-1].split('_')[0]+'_'+str(count)+'_Len_'+str(len(seq_rec.seq)))
			count += 1
	
	## keeps only CDSs that are greater than 30 bp (10 AA --> This is a cut-off in the 
	## phylogenomic pipeline too!)
	renamed_seqs = [i for i in renamed_seqs if len(i.split('\n')[-1]) > 30]
	
	print (color.BOLD+'\n\nFor '+color.DARKCYAN+args.input_file.split('/')[-1]+color.END+\
	color.BOLD+', '+color.RED+str(len(renamed_seqs))+' Protein sequences\n'+color.END+color.BOLD+\
	'were renamed while preserving the '+color.ORANGE+args.source+color.END+color.BOLD+' formatting'\
	+color.END+'\n')
	
	with open(args.input_file.replace('.fasta','.Prepped.fasta'),'w+') as w:
		w.write('\n'.join(renamed_seqs))

	with open('../'+args.input_file.split('/')[-1].replace('.fasta','.SeqCodes.tsv'),'w+') as w:
		w.write('Original Name\tNew Name\tSeq Length\n')
		for k, v in seq_code_dict.items():
			w.write(k+'\t'+v[0]+'\t'+str(len(v[1]))+'\n')


###########################################################################################
###--------------------- Cleans up the Folder and Moves Final Files --------------------###
###########################################################################################

def clean_up(args):

	os.system('rm '+args.input_file)
	os.system('mv ../'+args.input_file.split('/')[-1].replace('.fasta','.SeqCodes.tsv')+' '+\
	args.folder+'/Original/')
	os.system('mv '+args.input_file.replace('.fasta','.Prepped.fasta')+' '+args.folder)


###########################################################################################
###-------------------------------- Next Script Message --------------------------------###
###########################################################################################

def next_script(args):

	print (color.BOLD+'\nLook for '+color.DARKCYAN+args.input_file.split('/')[-1].replace('.fasta','.Renamed.fasta')\
	+'.fasta'+color.END+color.BOLD+'\nin the '+color.ORANGE+args.folder.split('/')[-1]+\
	' Folder\n\n'+color.END+color.BOLD)

	print ('Next Script is:\n\n'+color.PURPLE+'2p_CountOGsUsearch.py\n\n'+color.END)
	
	
##########################################################################################
###----------------------------- Calls on Above Functions -----------------------------###
##########################################################################################

def main():

	args = check_args()
	
	prep_folders(args)
		
	renamed_GenomeCDS(args)
	
	clean_up(args)
	
	next_script(args)
	
main()