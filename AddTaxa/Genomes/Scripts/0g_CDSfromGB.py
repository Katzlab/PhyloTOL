#!/usr/bin/env python3.5

##__Updated__: 31_08_2017
##__Author__: Xyrus Maurer-Alcala; maurerax@gmail.com
##__Usage__: python 0g_CDSfromGB.py --help


##############################################################################################
##																							##
## Script is designed to Extract CDS sequences (Nuc) from genbank entries from WGS projects	##
##																							##
## Special thanks to Jean-David Grattepanche for fixes and updates to this script!!			##
##																							##
##			E-mail Xyrus (author) for help if needed: maurerax@gmail.com					##
##																							##
##############################################################################################


### from __future__ is just in case someone is running python 2...

from __future__ import print_function
from Bio import SeqIO
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
	color.BOLD + '\n\nThis script is intended to extract '+color.RED+'Annotated ORFs '+\
	color.PURPLE+'ORFS\n'+color.END+color.BOLD+'from a provided Genbank formatted file.'\
	+usage_msg(), usage=SUPPRESS, formatter_class=RawTextHelpFormatter)
	
	required_arg_group = parser.add_argument_group(color.ORANGE+color.BOLD+'Required Options'+color.END)
 
	required_arg_group.add_argument('--input_file','-in', action='store',
	help=color.BOLD+color.GREEN+' Genbank formatted file (gbff or gb)\n'+color.END)

	optional_arg_group = parser.add_argument_group(color.ORANGE+color.BOLD+'Options'+color.END)

	
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
				
	return args
	
	
###########################################################################################
###------------------------------- Script Usage Message --------------------------------###
###########################################################################################

def usage_msg():
	return (color.BOLD+color.RED+'\n\nExample usage:'+color.CYAN+' python 0g_CDSfromGB.py'\
	' --input_file ../Stentor_coeruleus.gbff.gz'+color.END)


##########################################################################################
###-------- Storage for LARGE (Annoying) Print Statements for Flagged Options ---------###
##########################################################################################

def return_more_info(args):

	file_extensions = ['gbff','gbk','gb','genbank']

	author = (color.BOLD+color.ORANGE+'\n\n\tQuestions/Comments? Email Xyrus (author) at'\
	' maurerax@gmail.com\n\n'+color.END)

	if args.author == True:
		return author

	if args.input_file != None:
		if args.input_file.split('/')[-1] not in os.listdir(os.curdir+'./'):
			print (color.BOLD+'\n\nNo valid files were provided\n\nFor usage information (and example)'\
			' run this script with just: "--help"\n\n'+color.END)
			sys.exit()			
		elif args.input_file.split('.gz')[0].split('.')[-1] in file_extensions:
			pass
		else:
			print (color.BOLD+'\n\nDouble check that your file ends with an appropriate\n'\
			'file extension:\n\n'+color.GREEN+'\n'.join(file_extensions)+color.END+color.BOLD+\
			+'\n\nFor usage information (and example) run this script with just: "--help"\n\n'+\
			color.END)
	else:
		print (color.BOLD+'\n\nNo valid files were provided\n\nFor usage information (and example)'\
		' run this script with just: "--help"\n\n'+color.END)
		sys.exit()

##########################################################################################
###------------------------- Creates Folders For Storing Data -------------------------###
##########################################################################################

def prep_folders(args):
	
	home_folder = args.input_file.split('.gb')[0]+'.WGS.CDS.Prep'
	
	if os.path.isdir(home_folder) != True:
		os.system('mkdir '+home_folder)
		
	if os.path.isdir(home_folder+'/Original') != True:
		os.system('mkdir '+home_folder+'/Original')
		
	os.system('cp '+args.input_file+' '+home_folder+'/Original')	

	if args.input_file.endswith('.gz'):
		print (color.BOLD+'\n\nUnzipping the compressed GenBank file')
		os.system('gunzip '+args.input_file)
		args.input_file = args.input_file.rstrip('.gz')
	else:
		print (color.BOLD+'\n\nCompressing the GenBank file to save for later'+color.END)
		os.system('gzip '+home_folder+'/Original/'+args.input_file.split('/')[-1])


###########################################################################################
###---------------------- Extracts the CDSs from the GenBank File ----------------------###
###########################################################################################

def grab_CDS_entries(args):

	home_folder = args.input_file.split('.gb')[0]+'.WGS.CDS.Prep'

# List is where the CDSs will be stored
	ExtractedCDS = []

# open up the GB file then save the data to list... for ease
	print (color.BOLD+'\n\nParsing through data in the '+color.PURPLE+args.input_file\
	.split('/')[-1]+color.END+color.BOLD+'\nGenBank Formatted File\n\n'+color.END)

	inGB_file = [i for i in SeqIO.parse(args.input_file,'genbank')]

# runs through the genbank records and extracts the CDS, the locus name and the protein name
	print (color.BOLD+'Attemping to grab CDS-tagged sequences from the\n'+color.PURPLE+\
	args.input_file.split('/')[-1]+color.END+color.BOLD+' GenBank File\n\n'+color.END)
	
	for gb_rec in inGB_file:
		for feature in gb_rec.features:
			if feature.type == 'CDS':
				if 'translation' in feature.qualifiers:
					try:
						ExtractedCDS.append('>'+str(feature.qualifiers['locus_tag'][0])+\
						'_'+str(feature.qualifiers['protein_id'][0])+'\n'+str(feature.location.extract(gb_rec).seq))
						print (color.BOLD+'Parsed through '+color.ORANGE+str(len(ExtractedCDS))+\
						color.END+color.BOLD+' CDSs so far!', end='\r')
					except:
						try:
							ExtractedCDS.append('>'+str(feature.qualifiers['gene'][0])+\
							'_'+str(feature.qualifiers['protein_id'][0])+'\n'+str(feature.location.extract(gb_rec).seq))
							print (color.BOLD+'Parsed through '+color.ORANGE+str(len(ExtractedCDS))\
							+color.END+color.BOLD+' CDSs so far!', end='\r')
						except:
							print ("ERROR for:", feature.qualifiers['gene'][0], end='\r')

				
# internal check to make sure that there actually are data to be written out
	if len(ExtractedCDS) == 0:
		print (color.BOLD+color.RED+'\n\nError: No CDSs wer found in your GB file\n\n'\
		+'Check by eye that CDS features are actually present\n\n'+color.END)
		sys.exit()
	else:
		print (color.BOLD+'There are '+color.ORANGE+str(len(ExtractedCDS))\
			+' CDS '+color.END+color.BOLD+'records in the '+color.PURPLE\
			+args.input_file.split('/')[-1]+' GenBank'+color.END+color.BOLD\
			+' File provided\n\n'+color.END)

	with open(home_folder+'/'+home_folder.split('/')[-1].split('.Prep')[0]+'.fasta','w+') as w:
		for CDS_seq in ExtractedCDS:
			w.write(CDS_seq+'\n')

##########################################################################################
###--------------- Checks Command Line Arguments and Calls on Functions ---------------###
##########################################################################################
		
def main():
	
	args = check_args()
	
	prep_folders(args)
	
	grab_CDS_entries(args)

	os.system('rm '+args.input_file)
	
 	print (color.BOLD +'Next Script is: '+ color.PURPLE + '1g_RenameCDS.py\n\n'+ color.END)
	
main()