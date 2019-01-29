#!/usr/bin/env python3.5

##__Updated__: 19_09_2017
##__Author__: Xyrus Maurer-Alcala; maurerax@gmail.com
##__Usage__: python 2g_GCodeEval.py --help


#############################################################################################
#                                                                                           #
# Suggests which Genetic Code to use based upon Presence/Absence of Specific Stop Codons    #
# at the end of the CDS sequences. This is to provide a ROUGH gauge for the user.           #
#                                                                                           #
#############################################################################################


import argparse, os, sys
from argparse import RawTextHelpFormatter,SUPPRESS
from Bio import SeqIO
from Bio.Seq import Seq

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
	color.BOLD + '\n\nThis script is intended to aid you with '+color.RED+'evaluating\n(or checking) '+\
	color.END+color.BOLD+'the putative '+color.PURPLE+'Genetic Code'+color.END+color.BOLD+\
	' for a given\nFasta file of annotated (and untranslated) CDSs.\n\nTo do so, this script'\
	' checks for stop codon usages,\n'+color.RED+'suggesting '+color.END+color.BOLD+'the use of'\
	+color.PURPLE+' published and well-known\nalternate genetic codes'+color.END+color.BOLD+\
	' that are supported by the\nnext script: '+color.END+color.BOLD+color.PURPLE+'3g_GCodeTranslate.py'\
	+usage_msg(), usage=SUPPRESS, formatter_class=RawTextHelpFormatter)

	
	required_arg_group = parser.add_argument_group(color.ORANGE+color.BOLD+'Required Options'+color.END)
 
	required_arg_group.add_argument('--input_file','-in', action='store',
	help=color.BOLD+color.GREEN+' Fasta file with CDSs\n'+color.END)

	optional_arg_group = parser.add_argument_group(color.ORANGE+color.BOLD+'Options'+color.END)

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
		
	args.folder = '../'+args.input_file.split('/')[1]
		
	return args
	
	
###########################################################################################
###------------------------------- Script Usage Message --------------------------------###
###########################################################################################

def usage_msg():
	return (color.BOLD+color.RED+'\n\nExample usage:'+color.CYAN+' python 2g_GCodeEval.py'\
	' --input_file ../Stentor_coeruleus.WGS.CDS.Prep/Stentor_coeruleus.WGS.CDS.Renamed.fasta'+color.END)


##########################################################################################
###-------- Storage for LARGE (Annoying) Print Statements for Flagged Options ---------###
##########################################################################################

def return_more_info(args):

	valid_arg = 0

	supported_gcodes = ['Blepharisma\t(TGA = W)','Chilodonella\t(TAG/TGA = Q)','Ciliate\t\t(TAR = Q)',\
	'Conylostoma\t(TAR = Q, TGA = W)','Euplotes\t(TGA = C)','Peritrich\t(TAR = E)','None\t\t(TGA/TAG/TAA = X)',\
	'Universal\t(TGA/TAG/TAA = STOP)','TAA\t\t(TAG/TGA = Q)', 'TAG\t\t(TRA = Q)', 'TGA\t\t(TAR = Q)']

	author = (color.BOLD+color.ORANGE+'\n\n\tQuestions/Comments? Email Xyrus (author) at'\
	' maurerax@gmail.com\n\n'+color.END)

	if args.list_codes == True:
		print (color.BOLD+color.RED+'\nThese are the currently supported genetic codes.\n'+color.END)
		print (color.BOLD+color.ORANGE+'\n'.join(supported_gcodes)+'\n\n'+color.END)
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
###-------------------- Counts Several Metrics of Stop Codon Usage ---------------------###
###########################################################################################

def count_stops(args):

	print (color.BOLD+'\n\nScanning CDSs for In-Frame Stop Codons and Tracking\nFINAL '\
	'(Terminal) stop codon usage\n\n'+color.END)
	
	inFasta = [i for i in SeqIO.parse(args.input_file,'fasta')]
	seq_ends = [str(i.seq)[-3:].lower() for i in inFasta]
	inFrame_stops_raw = [str(i.seq[:-3].translate()).count('*') for i in inFasta]
	inFrame_stops_summary = [i for i in inFrame_stops_raw if i != 0]
				
	tga_end = seq_ends.count('tga')
	tag_end = seq_ends.count('tag')
	taa_end = seq_ends.count('taa')
	
	end_stop_freq = [tga_end, tag_end, taa_end]
	
	if max(end_stop_freq) > 0.95*sum(end_stop_freq):
		pos_to_keep = [i for i, j in enumerate(end_stop_freq) if j == max(end_stop_freq)][0]
	try:
		if pos_to_keep == 0:
			end_stop_freq = [end_stop_freq[0],0,0]
		elif pos_to_keep == 1:
			end_stop_freq = [0,end_stop_freq[1],0]
		elif pos_to_keep == 2:
			end_stop_freq = [0,0,end_stop_freq[2]]
	except:
		pass
		
	inFrame_stop_info = [len(inFrame_stops_summary), int(round(len(inFrame_stops_raw)*0.05)), sum(inFrame_stops_summary)]
	return end_stop_freq, inFrame_stop_info


###########################################################################################
###-------------------- Suggests Genetic Code Given Stop Codon Usage -------------------###
###########################################################################################

def suggest_code(args):

	stop_freq, inFrames = count_stops(args)

	genetic_code = ''

	if stop_freq.count(0) == 3:
		print (color.BOLD + color.RED + '\n\nNO Stop Codons Present in Data-set\n\n')
		genetic_code = 'None (UNDETERMINED -- NO STOP CODONS)'
	else:
	## DUMB way of checking if there are a significant (> 5%) number of CDSs with IN-FRAME stop codons 
		if inFrames[0] < inFrames[1]:
			print (color.BOLD + '\n\nSuggested Genetic Code is: '+color.CYAN+' Universal (table = 1)'+color.END)
			genetic_code = 'Universal (table = 1)'	
		else:
		
			if stop_freq[0] != 0 and stop_freq[1] != 0 and stop_freq[2] != 0:
				print (color.BOLD + '\n\nSuggested Genetic Code is: '+color.CYAN+' Condylostoma-Code'\
				' (No Dedicated Stops) OR None (all stops = "X")'+color.END)
				genetic_code = 'Condylostoma or None'
			if stop_freq[0] == 0 and stop_freq[1] == 0:
				print (color.BOLD + '\n\nSuggested Genetic Code is: '+color.CYAN+' Chilodonella-Code'\
				+' (Only Stop = TAA)'+color.END)
				genetic_code = 'Chilodonella or TAA'
			if stop_freq[0] == 0 and stop_freq[2] == 0:
				print (color.BOLD + '\n\nSuggested Genetic Code is: '+color.CYAN+' TAG-Code'\
				+' (Only Stop = TAG)'+color.END)
				genetic_code = 'TAG'
			if stop_freq[1] == 0 and stop_freq[2] == 0:
				print (color.BOLD + '\n\nSuggested Genetic Code is: '+color.CYAN+' Ciliate-Code'\
				+' (table = 6)'+color.END)
				genetic_code = 'Ciliate (table = 6)'
			if stop_freq[0] != 0 and stop_freq[1] != 0 and stop_freq[2] == 0:
				print (color.BOLD + '\n\nSuggested Genetic Code is: '+color.CYAN+' TGA/TAG are STOP'+color.END)
				genetic_code = 'TGA/TAG'
			if stop_freq[0] != 0 and stop_freq[1] == 0 and stop_freq[2] != 0:
				print (color.BOLD + '\n\nSuggested Genetic Code is: '+color.CYAN+' TGA/TAA are STOP'+color.END)
				genetic_code = 'TGA/TAA'
			if stop_freq[0] == 0 and stop_freq[1] != 0 and stop_freq[2] != 0:
				print (color.BOLD + '\n\nSuggested Genetic Code is: '+color.CYAN+' Blepharisma/Euplotes-Codes'\
				+color.END + color.BOLD+'\n--- NOTE: '+color.RED+' Stop-Codon Reassignments'\
				+' differ! (TGA = W or TGA = C)' + color.END)
				genetic_code = 'Blepharisma (TGA = W) or Euplotes (TGA = C)'

	return genetic_code, stop_freq
	

###########################################################################################
###---------------- Writes Out Currently Crummy Summary of Genetic Codes ---------------###
###########################################################################################

def summarize(args):

	suggestion, stop_freq = suggest_code(args)

	with open(args.input_file.split('.fa')[0]+'.GeneticCode.txt','w+') as w:
		w.write('Stop Codon\tFrequency\n')
		w.write('TGA\t'+str(stop_freq[0])+'\n')
		w.write('TAG\t'+str(stop_freq[1])+'\n')
		w.write('TAA\t'+str(stop_freq[2])+'\n\n')
		w.write('Suggestion For Genetic Code:\t'+suggestion+'\n\n')


##########################################################################################
###--------------- Checks Command Line Arguments and Calls on Functions ---------------###
##########################################################################################
		
def main():
	
	args = check_args()
		
	summarize(args)
		
	print (color.BOLD+'\nNext Script is: '+color.PURPLE+' 3g_GCodeTranslate.py\n\n'+color.END)

main()