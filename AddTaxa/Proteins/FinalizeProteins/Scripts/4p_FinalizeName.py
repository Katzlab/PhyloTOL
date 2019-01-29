#!/usr/bin/env python3.5

##__Updated__: 21_09_2017
##__Author__: Xyrus Maurer-Alcala; maurerax@gmail.com
##__Usage__: python 4p_FinalizeName.py --help

##################################################################################################
## This script is intended to rename the outputs of the FilterPartials script					##
## to a given 10-character that is used in the Katz lab Phylogenomic Tree building methods		##
##																								##
## Prior to running this script, ensure the following:											##
##																								##
## 1. You have assembled your transcriptome and COPIED the 'assembly' file 						##
##    (contigs.fasta, or scaffolds.fasta) to the PostAssembly Folder							##
## 2. Removed small sequences (usually sequences < 300bp) with ContigFilterPlusStats.py			##
## 3. Removed SSU/LSU sequences from your Fasta File											##
## 4. Classified your sequences as Strongly Prokaryotic/Eukaryotic or Undetermined				##
## 5. Classified the Non-Strongly Prokaryotic sequences into OGs 								##
## 6. You either know (or have inferred) the genetic code of the organism						##
## 7. You have translated the sequences and checked for the data in the RemovePartials folder	##
## 8. Partial sequences have been removed from the transcriptomic data sets						##
##																								##
## 										COMMAND Example Below									##
##									Extra Notes at Bottom of Script								##
##																								##
## 					E-mail Xyrus (author) for help if needed: maurerax@gmail.com				##
##																								##
##										Next Script(s) to Run: 									##
##                                     NONE! You're FINISHED! :D                                ##
##																								##
##################################################################################################

import argparse, os, sys
from argparse import RawTextHelpFormatter,SUPPRESS

#----------------------- Solely to Make Print Statements Colorful -----------------------#

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
	color.BOLD + '\n\nThis script is intended to '+color.RED+'Rename '+color.END\
	+color.BOLD+'the core set of '+color.PURPLE+'ORFS\n'+color.END+color.BOLD+'with a valid '\
	+color.RED+'10-character code'+color.END+color.BOLD+' for use in the KatzLab\nPhylogenomic Pipeline'\
	+usage_msg(), usage=SUPPRESS, formatter_class=RawTextHelpFormatter)
	
	required_arg_group = parser.add_argument_group(color.ORANGE+color.BOLD+'Required Options'+color.END)
 
	required_arg_group.add_argument('--input_file','-in', action='store',
	help=color.BOLD+color.GREEN+' The Fasta file to rename\n'+color.END)
	required_arg_group.add_argument('--name','-n', action='store',
	help=color.BOLD+color.GREEN+' A valid 10-Character code for updating the data\n'+color.END)


	optional_arg_group = parser.add_argument_group(color.ORANGE+color.BOLD+'Options'+color.END)

	optional_arg_group.add_argument('-author', action='store_true',
	help=color.BOLD+color.GREEN+' Prints author contact information\n'+color.END)

	if len(sys.argv[1:]) == 0:
		print (parser.description)
		print ('\n')
		sys.exit()

	args = parser.parse_args()
	
	if args.name == None:
		print (parser.description)
		print ('\n')
		sys.exit()	
	
	quit_eval = return_more_info(args)
	if quit_eval > 0:
		print ('\n')
		sys.exit()

	args.out_XML = args.name+'_XX_'+args.input_file.split('/')[-1]+'_1e-10keepall_BlastOutall.oneHit'
	
	args.file_prefix = args.input_file.split('/')[-1].split('_Filtered.Final')[0]

	if 'fasta' in args.file_prefix:
		args.file_prefix = args.name

	args.r2g_aa = '../../../ReadyToGo/ReadyToGo_AA/'
	args.r2g_ntd = '../../../ReadyToGo/ReadyToGo_NTD/'
	args.r2g_tsv = '../../../ReadyToGo/ReadyToGo_TSV/'
	args.r2g_xml = '../../../ReadyToGo/ReadyToGo_XML/'
	

	return args


###########################################################################################
###------------------------------- Script Usage Message --------------------------------###
###########################################################################################

def usage_msg():
	return (color.BOLD+color.RED+'\n\nExample usage:'+color.CYAN+' python 4p_FinalizeName.py'\
	' --input_file ../ToRename/Plasmopara_viticola_Filtered.Final.AA.fasta --name Sr_st_Pvit'+color.END)


##########################################################################################
###-------- Storage for LARGE (Annoying) Print Statements for Flagged Options ---------###
##########################################################################################

def return_more_info(args):

	valid_args = 0

	author = (color.BOLD+color.ORANGE+'\n\n\tQuestions/Comments? Email Xyrus (author) at'\
	' maurerax@gmail.com\n\n'+color.END)

	if args.author == True:
		print (author)
		valid_args += 1

	args.input_AA = args.input_file
	args.input_TSV = args.input_file.replace('.fasta','.allOGCleanresults.tsv')

	if os.path.isfile(args.input_AA) != True:
		print (color.BOLD+color.RED+'\nError:'+color.END+color.BOLD+' The provided Protein '\
		'Fasta file ('+color.DARKCYAN+args.input_AA.split('/')[-1]+color.END+color.BOLD+')\ndoes not'\
		' exist or is incorrectly formatted.\n\nDouble-check then try again!\n\n'+color.END) 
		valid_args += 1

	if os.path.isfile(args.input_TSV) != True:
		print (color.BOLD+color.RED+'\nError:'+color.END+color.BOLD+' The provided Nucleotide '\
		'Fasta file ('+color.DARKCYAN+args.input_TSV.split('/')[-1]+color.END+color.BOLD+')\ndoes not'\
		' exist or is incorrectly formatted.\n\nDouble-check then try again!\n\n'+color.END) 
		valid_args += 1

	return valid_args

###########################################################################################
###-------------------- Double Checks Format for 10-Character Code ---------------------###
###########################################################################################

def check_code(args):
	
	check_name = args.name.split('_')
	
	if len(args.name) != 10:
		print (color.BOLD+'\n\nNew Species Prefix is not 10 characters long\n\n')
		print ('Three examples below:\n'+color.CYAN+'\n\tSr_ci_Cunc\n\n\tOp_me_Hsap\n\n\t'\
		'Am_ar_Ehis\n\n'+color.END) 
		sys.exit()

	elif args.name.count('_') != 2:
		print (color.BOLD+'\n\nCheck the format of your Species Prefix!\n\n')
		print ('Three examples below:\n'+color.CYAN+'\n\tSr_ci_Cunc\n\n\tOp_me_Hsap\n\n\t'\
		'Am_ar_Ehis\n\n'+color.END) 

		sys.exit()

	if len(check_name[0]) == 2 and len(check_name[1]) == 2 and len(check_name[2]) == 4:
		print (color.BOLD+"\n\nRenaming "+color.ORANGE+args.input_file.split('/')[-1]\
		.split('_Filtered')[0]+color.END+color.BOLD+"'s files with the following 10-character\n"\
		"code: "+color.CYAN+args.name+color.END+'\n')
	else:
		print (color.BOLD+'\n\nCheck the format of your Species Prefix!\n\n')
		print ('Three examples below:\n'+color.CYAN+'\n\tSr_ci_Cunc\n\n\tOp_me_Hsap\n\n\t'\
		'Am_ar_Ehis\n\n'+color.END) 
		sys.exit()

			
##########################################################################################
###------------------------- Creates Folders For Storing Data -------------------------###
##########################################################################################

def prep_folders(args):
	

	if os.path.isdir(os.curdir+'./../../ReadyToGo/') != True:
		os.system('mkdir ../../../ReadyToGo')

	if os.path.isdir(os.curdir+'./../../ReadyToGo/ReadyToGo_NTD/') != True:
		os.system('mkdir '+args.r2g_ntd)
	if os.path.isdir(os.curdir+'./../../ReadyToGo/ReadyToGo_AA/') != True:
		os.system('mkdir '+args.r2g_aa)
	if os.path.isdir(os.curdir+'./../../ReadyToGo/ReadyToGo_TSV/') != True:
		os.system('mkdir '+args.r2g_tsv)
	if os.path.isdir(os.curdir+'./../../ReadyToGo/ReadyToGo_XML/') != True:
		os.system('mkdir '+args.r2g_xml)
		
	if os.path.isdir('../'+args.file_prefix) != True:
		os.system('mkdir ../'+args.file_prefix)	

	if os.path.isdir('../'+args.file_prefix+'/Renamed') != True:
		os.system('mkdir ../'+args.file_prefix+'/Renamed')
		
	if os.path.isdir('../../Finished/') != True:
		os.system('mkdir ../../Finished/')

###########################################################################################
###----------- Renames the NTD and AA CDSs with the Given 10-Character Code ------------###
###########################################################################################

def rename_paralogs(args):

	home_folder = '../'+args.file_prefix+'/Renamed/'

	print (color.BOLD+'\nRenaming Translated '+color.PURPLE+'Protein Sequences\n'+color.END)
	renamed_Final_Prots = open(args.input_AA).read().replace('>','>'+args.name+'_XX_')
		
	print (color.BOLD+'\nUpdating CDS Names in the Spreadsheet'+color.END)
	if '\n\n' in open(args.input_TSV).read():
		renamed_Final_tsv = args.name+'_XX_'+open(args.input_TSV).read().rstrip('\n')\
		.replace('\n\n','\n'+args.name+'_XX_')
	else:
		renamed_Final_tsv = args.name+'_XX_'+open(args.input_TSV).read().rstrip('\n')\
		.replace('\n','\n'+args.name+'_XX_')
		
	with open(home_folder+args.name+'_XX_'+args.input_AA.split('/')[-1],'w+') as w:
		w.write(renamed_Final_Prots)
	
	with open(home_folder+args.name+'_XX_'+args.input_TSV.split('/')[-1],'w+') as y:
		y.write(renamed_Final_tsv)


###########################################################################################
###--------------------------------- Header/Tail Lines ---------------------------------###
###########################################################################################

def header_tail():
	header = '<?xml version="1.0"?>\n<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">\n'\
		'<BlastOutput>\n  <BlastOutput_program>blastp</BlastOutput_program>\n  <BlastOutput_version>BLASTP 2.2.29+</BlastOutput_version>\n'\
		'  <BlastOutput_reference>Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&amp;auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), &quot;Gapped BLAST and PSI-BLAST: a new generation of protein database search programs&quot;, Nucleic Acids Res. 25:3389-3402.</BlastOutput_reference>\n'\
		'  <BlastOutput_db>../OGBlastDB/renamed_aa_seqs_OrthoMCL-5_12653.fasta</BlastOutput_db>\n  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>\n'

	tail = '</BlastOutput_iterations>\n</BlastOutput>'
	return header, tail


###########################################################################################
###------------------------------- TSV to XML Conversion -------------------------------###
###########################################################################################

def convert_TSV_data(args):

	home_folder = '../'+args.file_prefix+'/Renamed/'

	TSVforConvert = home_folder+args.name+'_XX_'+args.input_TSV.split('/')[-1]

	inTSV = [line.rstrip('\n') for line in open(TSVforConvert).readlines() if line != '\n']

	iterations = []

	for n in range(len(inTSV)):
		if n == 0:
			iterations.append('  <BlastOutput_query-def>'+inTSV[n].split('\t')[0]+'</BlastOutput_query-def>\n  <BlastOutput_query-len>'+str(abs(int(inTSV[n].split('\t')[-3])-int(inTSV[n].split('\t')[-4])+1))+'</BlastOutput_query-len>\n'\
				'  <BlastOutput_param>\n    <Parameters>\n      <Parameters_matrix>BLOSUM62</Parameters_matrix>\n      <Parameters_expect>1e-10</Parameters_expect>\n'\
				'      <Parameters_gap-open>11</Parameters_gap-open>\n      <Parameters_gap-extend>1</Parameters_gap-extend>\n      <Parameters_filter>F</Parameters_filter>\n'\
				'    </Parameters>\n  </BlastOutput_param>\n<BlastOutput_iterations>\n<Iteration>\n  <Iteration_iter-num>1</Iteration_iter-num>\n  <Iteration_query-ID>Query_1</Iteration_query-ID>\n'\
				'  <Iteration_query-def>'+inTSV[n].split('\t')[0]+'</Iteration_query-def>\n  <Iteration_query-len>'+str(abs(int(inTSV[n].split('\t')[-3])-int(inTSV[n].split('\t')[-4])+1))+'</Iteration_query-len>\n'\
				'<Iteration_hits>\n<Hit>\n  <Hit_num>1</Hit_num>\n  <Hit_id>Fake_Entry</Hit_id>\n  <Hit_def>'+inTSV[n].split('\t')[1]+'</Hit_def>\n  <Hit_accession>Fake_Accession</Hit_accession>\n'\
				'  <Hit_len>'+str(abs(int(inTSV[n].split('\t')[-3])-int(inTSV[n].split('\t')[-4])+1))+'</Hit_len>\n  <Hit_hsps>\n    <Hsp>\n      <Hsp_num>1</Hsp_num>\n      <Hsp_bit-score>1234</Hsp_bit-score>\n'\
				'      <Hsp_score>'+inTSV[n].split('\t')[-1]+'</Hsp_score>\n      <Hsp_evalue>'+inTSV[n].split('\t')[-2]+'</Hsp_evalue>\n      <Hsp_query-from>'+inTSV[n].split('\t')[-4]+'</Hsp_query-from>\n'\
				'      <Hsp_query-to>'+inTSV[n].split('\t')[-3]+'</Hsp_query-to>\n      <Hsp_hit-from>'+inTSV[n].split('\t')[-4]+'</Hsp_hit-from>\n      <Hsp_hit-to>'+inTSV[n].split('\t')[-3]+'</Hsp_hit-to>\n'\
				'      <Hsp_query-frame>0</Hsp_query-frame>\n      <Hsp_hit-frame>0</Hsp_hit-frame>\n      <Hsp_identity>'+str(abs(int(inTSV[n].split('\t')[-3])-int(inTSV[n].split('\t')[-4])))+'</Hsp_identity>\n'\
				'      <Hsp_positive>'+str(abs(int(inTSV[n].split('\t')[-3])-int(inTSV[n].split('\t')[-4])))+'</Hsp_positive>\n      <Hsp_gaps>0</Hsp_gaps>\n      <Hsp_align-len>'+str(abs(int(inTSV[n].split('\t')[-3])-int(inTSV[n].split('\t')[-4])))+'</Hsp_align-len>\n'\
				'      <Hsp_qseq></Hsp_qseq>\n      <Hsp_hseq></Hsp_hseq>\n      <Hsp_midline></Hsp_midline>\n    </Hsp>\n  </Hit_hsps>\n</Hit>\n'\
				'\n</Iteration_hits>\n  <Iteration_stat>\n    <Statistics>\n      <Statistics_db-num>379660</Statistics_db-num>\n      <Statistics_db-len>197499634</Statistics_db-len>\n'\
				'      <Statistics_hsp-len>123</Statistics_hsp-len>\n      <Statistics_eff-space>184705217500</Statistics_eff-space>\n      <Statistics_kappa>0.041</Statistics_kappa>\n'\
				'      <Statistics_lambda>0.267</Statistics_lambda>\n      <Statistics_entropy>0.14</Statistics_entropy>\n    </Statistics>\n  </Iteration_stat>\n</Iteration>\n')
		else:
			iterations.append('<Iteration>\n  <Iteration_iter-num>'+str(n+1)+'</Iteration_iter-num>\n  <Iteration_query-ID>Query_'+str(n+1)+'</Iteration_query-ID>\n'\
				'  <Iteration_query-def>'+inTSV[n].split('\t')[0]+'</Iteration_query-def>\n  <Iteration_query-len>'+str(abs(int(inTSV[n].split('\t')[-3])-int(inTSV[n].split('\t')[-4])+1))+'</Iteration_query-len>\n'\
				'<Iteration_hits>\n<Hit>\n  <Hit_num>1</Hit_num>\n  <Hit_id>Fake_Entry</Hit_id>\n  <Hit_def>'+inTSV[n].split('\t')[1]+'</Hit_def>\n  <Hit_accession>Fake_Accession</Hit_accession>\n'\
				'  <Hit_len>'+str(abs(int(inTSV[n].split('\t')[-3])-int(inTSV[n].split('\t')[-4])+1))+'</Hit_len>\n  <Hit_hsps>\n    <Hsp>\n      <Hsp_num>1</Hsp_num>\n      <Hsp_bit-score>1234</Hsp_bit-score>\n'\
				'      <Hsp_score>'+inTSV[n].split('\t')[-1]+'</Hsp_score>\n      <Hsp_evalue>'+inTSV[n].split('\t')[-2]+'</Hsp_evalue>\n      <Hsp_query-from>'+inTSV[n].split('\t')[-4]+'</Hsp_query-from>\n'\
				'      <Hsp_query-to>'+inTSV[n].split('\t')[-3]+'</Hsp_query-to>\n      <Hsp_hit-from>'+inTSV[n].split('\t')[-4]+'</Hsp_hit-from>\n      <Hsp_hit-to>'+inTSV[n].split('\t')[-3]+'</Hsp_hit-to>\n'\
				'      <Hsp_query-frame>0</Hsp_query-frame>\n      <Hsp_hit-frame>0</Hsp_hit-frame>\n      <Hsp_identity>'+str(abs(int(inTSV[n].split('\t')[-3])-int(inTSV[n].split('\t')[-4])))+'</Hsp_identity>\n'\
				'      <Hsp_positive>'+str(abs(int(inTSV[n].split('\t')[-3])-int(inTSV[n].split('\t')[-4])))+'</Hsp_positive>\n      <Hsp_gaps>0</Hsp_gaps>\n      <Hsp_align-len>'+str(abs(int(inTSV[n].split('\t')[-3])-int(inTSV[n].split('\t')[-4])))+'</Hsp_align-len>\n'\
				'      <Hsp_qseq></Hsp_qseq>\n      <Hsp_hseq></Hsp_hseq>\n      <Hsp_midline></Hsp_midline>\n    </Hsp>\n  </Hit_hsps>\n</Hit>\n'\
				'\n</Iteration_hits>\n  <Iteration_stat>\n    <Statistics>\n      <Statistics_db-num>379660</Statistics_db-num>\n      <Statistics_db-len>197499634</Statistics_db-len>\n'\
				'      <Statistics_hsp-len>123</Statistics_hsp-len>\n      <Statistics_eff-space>184705217500</Statistics_eff-space>\n      <Statistics_kappa>0.041</Statistics_kappa>\n'\
				'      <Statistics_lambda>0.267</Statistics_lambda>\n      <Statistics_entropy>0.14</Statistics_entropy>\n    </Statistics>\n  </Iteration_stat>\n</Iteration>\n')

	return iterations


###########################################################################################
###--------------------------- Writes Out the Fake XML File ----------------------------###
###########################################################################################

def write_Fake_XML(args):

	home_folder = '../'+args.file_prefix+'/'

	print (color.BOLD+'\n\nConverting '+color.ORANGE+args.input_file.split('/')[-1]+color.END\
	+color.BOLD+' to XML format\n'+color.END)

	header, tail = header_tail()
	
	iterations = convert_TSV_data(args)	
	
	with open(home_folder+args.out_XML,'w+') as w:
		w.write(header)
		w.write(''.join(iterations))
		w.write(tail)
		
##########################################################################################
###-------------------- Cleans up the Folder and Moves Final Files --------------------###
##########################################################################################
def clean_up(args):

	folders = []
	for f in os.listdir(os.curdir+'./'+args.file_prefix+'/Original/'):
		if os.path.isdir('../'+args.file_prefix+'/Original/'+f) and args.file_prefix in f:  
			folders.append(f)

	home_folder = '../'+args.file_prefix+'/Renamed/'
	
	os.system('cp '+home_folder+'*fasta ../'+args.file_prefix)
	os.system('cp '+home_folder+'*tsv ../'+args.file_prefix)

	os.system('cp ../'+args.file_prefix+'/*_XX_*fasta '+args.r2g_aa)
	os.system('cp ../'+args.file_prefix+'/'+args.out_XML+' '+args.r2g_xml)
	os.system('cp ../'+args.file_prefix+'/*tsv '+args.r2g_tsv)
	
	os.system('cp -r ../'+args.file_prefix+' ../../Finished/'+args.name+'_XX_'+args.file_prefix)
	
	os.system('rm -rf ../'+args.file_prefix)
	
 	os.system('rm ../ToRename/*'+args.file_prefix+'*')

	if len(os.listdir('../ToRename/')) == 0:
		os.system('rm -r ../ToRename')
	elif '.DS_Store' in os.listdir('../ToRename/') and len(os.listdir('../ToRename/')) == 1:
		os.system('rm -r ../ToRename')

	for f in os.listdir(os.curdir+'./../ToFinalize/'):
		if f in folders:
			os.system('rm -r ../../ToFinalize/'+f)
			
	if len(os.listdir(os.curdir+'./../ToFinalize/')) == 0:
		os.system('rm -r ../../ToFinalize')
	elif '.DS_Store' in os.listdir('../../ToFinalize/') and len(os.listdir('../../ToFinalize/')) == 1:
		os.system('rm -r ../../ToFinalize')


	
###########################################################################################
###-------------------------------- Next Script Message --------------------------------###
###########################################################################################

def next_script(args):

	print (color.BOLD+'\nThere is no next script! The final '+color.ORANGE+args.out_XML\
	.split('_XX')[0]+color.END+color.BOLD+' files can be\nfound in the '+color.RED+\
	args.out_XML.split('_XX_')[-1].split('_Filtered')[0]+color.END+color.BOLD+' and '\
	+color.RED+'ReadyToGo folders'+color.END+color.BOLD+' and are ready\n'\
	'for the KatzLab Phylogenomic Tree-Building Steps!\n\n'+color.END)

##########################################################################################
###--------------- Checks Command Line Arguments and Calls on Functions ---------------###
##########################################################################################
			
def main():

	args = check_args()
	
	check_code(args)
	
	prep_folders(args)
	
	rename_paralogs(args)
	
	write_Fake_XML(args)
	
	clean_up(args)
	
	next_script(args)
	
main()