#!/usr/bin/python
from Pipeline import Pipeline
import os, sys
import Utilities

def get_parameters():
	PathtoFiles = ''					##Path to folders ncbiFiles, BlastFiles, allOGFiles and the OGs of interest file
	testPipelineList = ''			##Name of the file with the OGs of interest -  this file needs to be in the same location as your datafile folders
	listTaxaInterest = ''
	blastCutOff = 1e-15				##Blast e-value cutoff (also depends on how blasts were run)
	seqLenCompCutOff = 20			##We remove very short sequences during the Taxon step - this is the amino acid length cutoff
	tooSimCutOff = 0.98				##Sequences are removed if they are more similar than this cutoff (default = 0.98)
	guidanceIter = 'n'				##number of Guidance iterations for sequences removal. the default is running until no bad scores remain (default = n)
	seqcutoff = 0.4 				##During guidance, taxa are removed if their score is below this cutoff (default as of 11/2014 = .35)
	colcutoff = 0.4					##During guidance, taxa are removed if their score is below this cutoff (default as of 11/2014 = .40)
	rescutoff = 0.4					##During guidance, residues are removed if their score is below this cutoff (default as of 03/2017 = .40)  # MAC---I added this v3
	concatAlignment = ''					##Do you want to remove paralogs and generate an alignment for concatenation? (y/n)  # MAC---I added this v3

	#try:
	infile = open('pipeline_parameter_file.txt','r')
	
	for line in infile:
		if line[0] == '#':
			attribute = line.split()[0].strip('#')
			value = line.split()[2].strip()
			if attribute == 'PathtoFiles':
				PathtoFiles = value
			elif attribute == 'testPipelineList':
				testPipelineList = value
			elif attribute == 'listTaxaInterest':
				listTaxaInterest = value
			elif attribute == 'blastCutOff':
				blastCutOff = float(value)
			elif attribute == 'seqLenCompCutOff':
				seqLenCompCutOff = int(value)			
			elif attribute == 'tooSimCutOff':
				tooSimCutOff = float(value)						
			elif attribute == 'guidanceIter':
				guidanceIter = value
			elif attribute == 'seqcutoff':
				seqcutoff = float(value)			
			elif attribute == 'colcutoff':
				colcutoff = float(value)
			elif attribute == 'rescutoff':
				rescutoff = float(value)
			elif attribute == 'concatAlignment':
				concatAlignment = value	

	#except:
	#	print('no parameter file found - using default parameters. ')
		
	return PathtoFiles,testPipelineList,listTaxaInterest,blastCutOff,seqLenCompCutOff,tooSimCutOff,guidanceIter,seqcutoff,colcutoff,rescutoff,concatAlignment


def writelog(PathtoOutput,string):
	logfile = open(PathtoOutput + '/logfile', 'a')
	logfile.write(string + '\n')
	logfile.close()
	
def main():
	arg = sys.argv
	ct = 'n'
	if len(arg) == 2:
		mode = sys.argv[1]
		allowed_modes = ['ng', 'nr']
		if mode == 'ct' : ct = 'y'
		if mode not in allowed_modes : mode = 'df'
	else:
		mode = 'df'
	print("\n** mode -> %s **" % mode)
		
	PathtoFiles,testPipelineList,listTaxaInterest,blastCutOff,seqLenCompCutOff,tooSimCutOff,guidanceIter,seqcutoff,colcutoff,rescutoff,concatAlignment = get_parameters()
	paramList = [blastCutOff,seqLenCompCutOff,tooSimCutOff,guidanceIter,seqcutoff,colcutoff,rescutoff,concatAlignment]

	print('################################################################################')
	print('')
	print('KATZLAB PHYLOGENOMICS PIPELINE')
	print('')
	print('This script assumes your data files and scripts are in the folders they came in. ')
	print('It also assumes there is a list of OGs you are interested in is in the Files folder.')
	print('')
	print('')
	print('PARAMETERS:')
	print('name of OG list = %s' % testPipelineList )
	print('name of list of taxa of interest = %s' % listTaxaInterest )
	print('Blast cutoff = %s' %  blastCutOff)
	print('Sequence length cutoff = %s' % seqLenCompCutOff)
	print('Cluster cutoff (too similar) = %s' % tooSimCutOff)
	print('Number of Guidance iterations = %s' % guidanceIter)
	print('Guidance sequence cutoff = %s' % seqcutoff)
	print('Guidance colum cutoff = %s' % colcutoff)
	print('Guidance residue cutoff = %s' % rescutoff)
	print('Alignment for concatenation = %s' % concatAlignment)
	print('################################################################################')
		
	if concatAlignment is not 'y' and concatAlignment is not 'n':
		print("\n*** your answer concatAlignment = " + concatAlignment + " is not correct. The pipeline takes 'n' as default ***")
	
	if not ct == 'y':
		if os.path.exists('../' + testPipelineList + '_results2keep'):
			print('terminating PhyloTOL: the folder ' + '../' + testPipelineList + '_results2keep exists. Choose another name for your OG list\n\n')
			quit()

	infile = open(PathtoFiles + '/' + testPipelineList,'r').readlines()  #list of ogs of interest
	if infile == []: 
		print('terminating PhyloTOL: Your list of OGs is empty\n\n')
		quit()
	
	'''
	MACR - For pipeleline 3
 
	The next part of the code produces two lists of taxa. One list is the list of taxa that the user wants 
	to process through the pipeline. The other list is the list of whole genomes that is used in some other 
	parts of the pipeline in the gene step. The data that is required is a filed called taxaDBpipeline3 that 
	should be in DataFiles and a list of taxa of interest provided for the user. The file taxaDBpipeline3
	says if every taxon of the pipeline is whole genome and if it is prokaryote or eukaryote.
	'''
		
	taxaDBfile = open(PathtoFiles + 'taxaDBpipeline3', 'r')
	taxaDBfile = taxaDBfile.readlines()
	taxaDB = [taxon.strip('\n') for taxon in taxaDBfile]
	wholegenomeDB = [taxon[3:] for taxon in taxaDB if taxon.split(',')[0] == 'wg']

	
	if not os.path.exists(PathtoFiles + listTaxaInterest):
		print("you need to have a list of taxa of interest")
		quit()
	else:
		listTOI = open('%s%s' % (PathtoFiles, listTaxaInterest), 'r').readlines()
		taxa2SF = []
		if listTOI[0] == "all\n":
			taxa2analyze = 'all'
			print("you chose to run your analysis with all taxa\n\n")
		else:
			taxaInterest = []
			sf = ''
			
			## MACR - from the taxa list specified for user take all taxa that match the database
			## until it finds '#' as 'taxa to be analysed'. Then, take all taxa that follow # as 
			## 'taxa to apply similarity filter, SF'
			
			for taxon in listTOI:
				taxon = taxon.strip('\n')
				if '#' in taxon:
					sf = 'y'
				else:
					if sf is not 'y':
						if not taxon.startswith('-'):
							for taxonINdb in taxaDB:
								if taxon in taxonINdb:
									taxaInterest.append(taxonINdb)
					else:
						for taxonINdb in taxaDB:
							if taxon in taxonINdb:
								taxa2SF.append(taxonINdb.split(',')[2])
						
			taxaInterest2 = list(taxaInterest)
			sf = ''
			for taxon in listTOI:
				if '#' in taxon:
					sf = 'y'
				else:
					if taxon.startswith('-'):
						taxon = taxon[1:].strip('\n')
						for taxonInterest in taxaInterest2:
							if taxon in taxonInterest:
								taxaInterest.remove(taxonInterest)
					
			taxaInterest = list(set(taxaInterest))

			if taxaInterest:
				taxa2analyze = []
				for taxon in taxaInterest:
					taxon2analyze = taxon.split(',')[2]
					taxa2analyze.append(taxon2analyze)
			else:
				taxa2analyze = []
		
			if len(taxa2analyze) == 0:
				print("none of your taxa of interest are in the pipeline database\n\n")
				quit()
			else:
				print("%s taxa will be analized\n\n%s\n\n" % (len(taxa2analyze), taxa2analyze))
				if taxa2SF:
					print("Similarity filter will be applied to these taxa:\n\n%s\n\n" % taxa2SF)
					

# MACR - Creating files, folders and writing logfiles.
				
	PathtoOutput =  '../my-data/' + testPipelineList + '_results/Output/'
	os.system('mkdir ../my-data/')
	os.system('mkdir ../my-data/' + testPipelineList + '_results')
	os.system('mkdir ../my-data/' + PathtoOutput)
	writelog(PathtoOutput,'testPipelineList = ' + testPipelineList)
	writelog(PathtoOutput,'blastCutOff = ' + str(blastCutOff))
	writelog(PathtoOutput,'seqLenCompCutOff = ' + str(seqLenCompCutOff))
	writelog(PathtoOutput,'tooSimCutOff = ' + str(tooSimCutOff))
	writelog(PathtoOutput,'guidanceIter = ' + str(guidanceIter))
	writelog(PathtoOutput,'seqcutoff = ' + str(seqcutoff))
	writelog(PathtoOutput,'colcutoff = ' + str(colcutoff))
	writelog(PathtoOutput,'rescutoff = ' + str(rescutoff))
	writelog(PathtoOutput,'concatAlignment = ' + concatAlignment + ' (y = remove paralogs and make alignment, n = keep paralogs and do not make alignment)')
	
	# MACR 03/04/19 -- added this for calculating og average length for OF and SF
	oglengths = open(PathtoOutput + "oglengths", "a")
	ogs = open(PathtoFiles + "/" + testPipelineList, "r").readlines()
	for og in ogs:
		og = og.strip()
		seq_len = {}
		ogFile = open(PathtoFiles + "/allOG5Files/" + og, "r").readlines()
			
		for line in ogFile:
			line = line.strip()

			if line.startswith(">"):
				tag = line
				seq_len[tag] = 0
			else:
				seq_len[tag] += len(line)
			
		og_totalLength = 0
		for seqLength in seq_len.values() : og_totalLength += seqLength
		averageLength = og_totalLength / len(seq_len.values())
		oglengths.write("%s\t%s\n" % (og, averageLength))
	oglengths.close()
	
	'''
	MACR - Taxon step with changes for Pipeline 3

	The next part of the code calls the Taxon class. The aim is to generate a folder the folder fasta2keep. That folder
	contains sequences from non-orthomcl taxa categorized as OGs. Here there are modifications for processing only the 
	non-orthomcl taxa that are included in the list of taxa interest of the user. 
	'''

	for f in os.listdir(PathtoFiles + '/ncbiFiles'):

	# MACR - Pipeline 3: Given that the use provides a list of taxa of interest. Only takes the Blast reports of the taxa that match the list

		if taxa2analyze is not 'all':		
			taxonBlast = f[:10]
			if taxonBlast in taxa2analyze:
				print '\n' + f + '\n'
				if f[0] != '.':
					try:
						newPipe = Pipeline(PathtoFiles + testPipelineList, PathtoFiles, ('queueTaxa',f),paramList,taxa2analyze,taxa2SF,wholegenomeDB,mode)
					except Exception, e:
						elog = open('errorlog','a')
						elog.write(f + " failed on %s with: %s" % (f, e.message))
						elog.close()
						print ("failed on %s with: %s" % (f, e.message))
		else:
	
 	# MACR - Pipeline 3: Given that the user does not provide a list of taxa of interest. Run pipeline for all taxa
	
			print '\n' + f + '\n'
			if f[0] != '.':
				try:
					newPipe = Pipeline(PathtoFiles + testPipelineList, PathtoFiles, ('queueTaxa',f),paramList,taxa2analyze,taxa2SF,wholegenomeDB,mode)
				except Exception, e:
					elog = open('errorlog','a')
					elog.write(f + " failed on %s with: %s" % (f, e.message))
					elog.close()
					print ("failed on %s with: %s" % (f, e.message))
		
	'''
	JG - Pipeline 3
	
	run Gene step - just up to guidance.
	To do this, I made a new pipeline method (two, actually) called 'test_keep' and 'test_remove'
	keep and remove are for ingroup paralogs, as I wasn't sure which you wanted.
	
	MACR - Pipeline 3
	
	Besides the changes described by Jessica above, there 3 major changes for pipeline 3:
	- A new Guidance method: This is a bash script made by Miguel Fonseca
	- A "helper" for the new Guidance method. This helper is a perl script that edits the intermediary
	files of each guidance run in order to allow the looping.
	- A new module called Utilities. I made this module to take the preguidance files and run Miguel's 
	scripts. Once Guidance and raxml are done. Some functions of the module allow to continue to produce the alignments
	for concatenation. Finally it organizes all output files.
	
	MACR - Pipeline 3.1

	- Neddle step and ingroup paralogs removal were replaced. Now we have an overlap filter (OF) and a similarity filter (SF). Both are performed using 
	Usearch-Ublast. We removed this step from the taxon class and added as a separated script called 'iterUblast.py'. This 
	script is called from the module 'Utilities'.
	
	The logic now is: All sequences per taxa should be 1.5 times smaller tha the average OG legth and pass the overlap filter. Then the user specify 
	if she/he wants to run a similarity filter for the sequences.
	'''	

	os.system('mkdir ' + PathtoFiles+ '/FileLists_' + testPipelineList)	
	count = 0
	li = []
	for line in infile:
		if 'OG5_' in line:
			outfile = open(PathtoFiles+ '/FileLists_' + testPipelineList + '/list' + str(count),'w')
			li.append('list' + str(count))
			outfile.write(line)
			outfile.close()
			count = count + 1

	# MAC - for pipeline 3.1 all methods for gene step were replaced by this one
	try:
		for f in os.listdir(PathtoFiles+ '/FileLists_' + testPipelineList):			
			newPipe = Pipeline(PathtoFiles +'/'+ testPipelineList, PathtoFiles, ('geneStep',f),paramList,taxa2analyze,taxa2SF,wholegenomeDB,mode)	
		
		answer_Cleaner = ''
		valid_answers = ['y', 'n']

		if ct == 'y':
			answer_Cleaner = 'y'
		else:
			while (answer_Cleaner not in valid_answers):  
				answer_Cleaner = raw_input("\n\nDo you want to execute the cleaner? (y/n): ")
	
				if (answer_Cleaner not in valid_answers): 
					print("\n\nplease answer y or n")
	
		if answer_Cleaner is 'y':
			Utilities.cleaner(testPipelineList, PathtoFiles, PathtoOutput)
		
	except Exception, e:
		elog = open('errorlog','a')
		line = open(PathtoFiles+ '/FileLists_' + testPipelineList + '/' + f,'r').read()
		elog.write(line + " failed on %s with: %s" % (f, e.message))
		elog.close()
		print ("failed on %s with: %s" % (f, e.message))
	return True				
main()
