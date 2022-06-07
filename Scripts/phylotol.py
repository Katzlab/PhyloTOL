#!/usr/bin/python
from Pipeline import Pipeline
import os, sys
import Utilities
from ConcatUtils import remove_paralogs


infile = open('pipeline_parameter_file_testing.txt','r')

for line in infile:
	if line[0] == '#':
		attribute = line.split()[0].strip('#')
		try:
			value = line.split()[2].strip()
		except IndexError:
			continue
		if attribute == 'PathtoFiles':
			PathtoFiles = value
		elif attribute == 'resume':
			resumer_param = value
		elif attribute == 'run_until':
			run_until_stage = value
		elif attribute == 'resume_status':
			resumer_status_param = value
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
		elif attribute == 'majorClades':
			majorClades = value
		
		elif(attribute == 'query_clades'):
			cont_query_clades = value
		elif(attribute == 'sister_clades'):
			cont_sister_clades = value
		elif(attribute == 'break_up'):
			cont_break_up = value
		elif(attribute == 'branch_length_filter'):
			cont_branch_length_filter = value
		elif(attribute == 'single_sister_only'):
			cont_single_sister_only = value
			
		elif(attribute == 'use_contloop'):
			use_contloop = value
		elif(attribute == 'seq_or_clade'):
			cont_seq_or_clade = value
		elif(attribute == 'n_loops'):
			max_num_contloops = value
		
		elif(attribute == 'rules_file'):
			cont_rules_file = value
		elif(attribute == 'remove_nm'):
			remove_nm = value
		
		elif(attribute == 'target_clade'):
			cont_target_clade = value
		elif(attribute == 'num_contams'):
			cont_num_contams = value
		elif(attribute == 'min_target_presence'):
			cont_min_target_presence = value
		elif(attribute == 'target_taxa_file'):
			cont_target_taxa_file = value
		elif(attribute == 'at_least_file'):
			cont_at_least_sisters_file = value
		elif(attribute == 'at_least_sisters_num'):
			cont_at_least_sisters_num = value
		elif(attribute == 'include_outgroups'):
			include_outgroups_contloop = value
		elif(attribute == 'return_subtrees'):
			return_cladegrabbing_subtrees = value
		elif(attribute == 'filter_differential_coverage'):
			cont_filter_differential_coverage = value
		elif(attribute == 'coverage_diff_OM'):
			cont_coverage_diff_OM = value
		elif(attribute == 'coverage_diff_abs'):
			cont_coverage_diff_abs = value

		elif(attribute == 'tree_font_size'):
			tree_font_size = value		

def writelog(PathtoOutput,string):
	logfile = open(PathtoOutput + '/logfile', 'a')
	logfile.write(string + '\n')
	logfile.close()	

	#run_resumer(resumer_param, 'unaligned', testPipelineList, guidanceIter, seqcutoff, colcutoff, rescutoff, contloop_mode = True, ogList = trees2reprocess)
def run_resumer(resumer_param, resumer_status_param, oglistName, guidanceIter, seqcutoff, colcutoff, rescutoff, contloop_mode = False, ogList = None):
	
	if(resumer_status_param == 'unaligned'):
		raxmlTemp = resumer_param + "/RAxML/"	
		PathtoOutput = resumer_param + "/out_resume/"
		os.system('mkdir ' + PathtoOutput)
		os.system('mkdir ' + PathtoOutput + "/Guidance")

		for file in os.listdir(resumer_param):
			if 'postguidance' not in file and (file.endswith('.fasta') or file.endswith('.fas') or file.endswith('.fna') or file.endswith('.fa')):					
				if(contloop_mode):
					keyw = 'OG5_' + file.split('OG5_')[-1][:6]
				else:
					keyw = file

				if((contloop_mode and ogList != None and keyw in ogList) or not contloop_mode):
					os.system("cp " + resumer_param + '/' + file + ' '  + PathtoOutput + '/Guidance/'+ keyw + 'forGuidance.fas')
					path2og = PathtoOutput + '/Guidance/'+ keyw + 'forGuidance.fas'
					Utilities.iterGuidance(oglistName, keyw, PathtoOutput, guidanceIter, seqcutoff, colcutoff, rescutoff, 'df')
			else:
				print('\nInput un-aligned file ' + file + ' is improperly formatted. It must be a .fasta, .fas, .fna, or .fa file. Skipping this alignment\n')

	elif(resumer_status_param == 'aligned'):
		raxmlTemp = resumer_param + "/RAxML/"
		if(not os.path.isdir(raxmlTemp)):
			os.mkdir(raxmlTemp)	
		for aln in os.listdir(resumer_param):
			if aln.endswith('.fasta') or aln.endswith('.fas') or aln.endswith('.fna') or aln.endswith('.fa'):
				os.system("cp " + resumer_param + '/' + aln + " " + raxmlTemp)
			else:
				print('\nInput alignment ' + aln + ' is improperly formatted. It must be a .fasta, .fas, .fna, or .fa file. Skipping this alignment\n')
	
		for aln in os.listdir(raxmlTemp):
			if aln.endswith('.fasta') or aln.endswith('.fas') or aln.endswith('.fna') or aln.endswith('.fa'):
				input = raxmlTemp + aln
				output =  aln[:10] + '_postguidance.fas' 
				#os.system("raxmlHPC-AVX2 -s " + input + " -m PROTGAMMALG -f d -p 12345 -# 1 -w " + raxmlTemp + " -n " + output)
				os.system("raxmlHPC-PTHREADS-AVX2 -s " + input + " -m PROTGAMMALG -f d -p 12345 -# 1 -w " + raxmlTemp + " -n " + output + " -T 3")
				os.system("mv " + raxmlTemp + "RAxML_bestTree." + output + " " + resumer_param + "/" + "RAxML_bestTree." + output + ".tre")
	elif(resumer_status_param != 'trees' and resumer_status_param != 'unaligned' and resumer_status_param != 'aligned'):
		print('\nInvalid resumer status. Please either turn off the resumer by setting the "resume" parameter to "NA" or input "aligned" or "unaligned" as the resumer status.\n')
		quit()

	
def main():

	allowed_modes = ['ng', 'nr', 'all', 'clean']
	if(run_until_stage not in allowed_modes):
		print('\nYour input run_until_stage is invalid. Please enter "ng" to run through preguidance only, "nr" to run through postguidance only, "all" to run through tree building, or "clean" to clean up the file structure\n')
		exit()
		
	paramList = [blastCutOff,seqLenCompCutOff,tooSimCutOff,guidanceIter,seqcutoff,colcutoff,rescutoff,concatAlignment]

	print("\n** mode -> %s **" % run_until_stage)

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
	print('Major clades = %s' % majorClades)
	print('################################################################################')

	if concatAlignment is not 'y' and concatAlignment is not 'n':
		print("\n*** your answer concatAlignment = " + concatAlignment + " is not correct. The pipeline takes 'n' as default ***")
	
	#ACL -- added the resumer feature to this main script instead of having two separate scripts
	if resumer_param != 'NA':
		if(os.path.isdir(resumer_param)):
			run_resumer(resumer_param, resumer_status_param, testPipelineList, guidanceIter, seqcutoff, colcutoff, rescutoff)
		else:
			print('\nTrying to resume PhyloToL, but the input folder of intermediate could not be found. Try fixing this path in the parameter file or setting the "resume" parameter to "NA"\n')
			quit()
	else:
		if os.path.exists('../' + testPipelineList + '_results2keep'):
			print('terminating PhyloTOL: the folder ' + '../' + testPipelineList + '_results2keep exists. Choose another name for your OG list\n\n')
			quit()

		infile = open(PathtoFiles + '/' + testPipelineList,'r').readlines()  #list of ogs of interest
		if infile == []: 
			print('terminating PhyloTOL: Your list of OGs is empty\n\n')
			quit()		

		# MACR - Creating files, folders and writing logfiles.
				
		PathtoOutput =  '../my-data/' + testPipelineList + '_results/Output/'
		os.system('mkdir ../my-data/')
		os.system('mkdir ../my-data/' + testPipelineList + '_results')
		os.system('mkdir ../my-data/' + PathtoOutput)
		writelog(PathtoOutput,'mode = ' + run_until_stage)
		writelog(PathtoOutput,'testPipelineList = ' + testPipelineList)
		writelog(PathtoOutput,'blastCutOff = ' + str(blastCutOff))
		writelog(PathtoOutput,'seqLenCompCutOff = ' + str(seqLenCompCutOff))
		writelog(PathtoOutput,'tooSimCutOff = ' + str(tooSimCutOff))
		writelog(PathtoOutput,'guidanceIter = ' + str(guidanceIter))
		writelog(PathtoOutput,'seqcutoff = ' + str(seqcutoff))
		writelog(PathtoOutput,'colcutoff = ' + str(colcutoff))
		writelog(PathtoOutput,'rescutoff = ' + str(rescutoff))
		writelog(PathtoOutput,'concatAlignment = ' + concatAlignment + ' (y = remove paralogs and make alignment, n = keep paralogs and do not make alignment)')
		writelog(PathtoOutput,'majorClades = ' + str(majorClades))
	
		# MACR -- V4 -- Added this method here for cleaning intermediary files and logs (for instance, after incomplete run or forced stoppage) using phylotol
		if run_until_stage == "clean": 
			Utilities.cleaner(testPipelineList, PathtoFiles, PathtoOutput)
			print("Cleaning folders -- done!")
			quit()
		
		try:
			taxaDBfile = open(PathtoFiles + 'taxaDBpipeline4', 'r')
		except IOError:
			print('\nYour taxaDBpipeline4 file could not be found! This should contain a list of all of the available taxa and their attributes (nw/wg/om, euk/prok)\n.')
			exit()
		taxaDBfile = taxaDBfile.readlines()
		taxaDB = [taxon.strip('\n') for taxon in taxaDBfile]

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
					print("%s taxa will be analyzed\n\n%s\n\n" % (len(taxa2analyze), taxa2analyze))
					if taxa2SF:
						print("Similarity filter will be applied to these taxa:\n\n%s\n\n" % taxa2SF)
	
		# MACR 03/04/19 -- added this for calculating og average length in the HookDB for OF and SF
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
					print('\n' + f + '\n')
					if f[0] != '.':
						try:
							newPipe = Pipeline(PathtoFiles + testPipelineList, PathtoFiles, ('queueTaxa',f),paramList,taxa2analyze,taxa2SF,majorClades,run_until_stage)
						except Exception as e:
							elog = open('errorlog','a')
							elog.write(f + " failed on %s with: %s" % (f, e))
							elog.close()
							print ("failed on %s with: %s" % (f, e))
			else:
	
		# MACR - Pipeline 3: Given that the user does not provide a list of taxa of interest. Run pipeline for all taxa
	
				print('\n' + f + '\n')
				if f[0] != '.':
					try:
						newPipe = Pipeline(PathtoFiles + testPipelineList, PathtoFiles, ('queueTaxa',f),paramList,taxa2analyze,taxa2SF,majorClades,run_until_stage)
					except Exception as e:
						elog = open('errorlog','a')
						elog.write(f + " failed on %s with: %s" % (f, e))
						elog.close()
						print ("failed on %s with: %s" % (f, e))
		
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
			if 'OG' in line:
				outfile = open(PathtoFiles+ '/FileLists_' + testPipelineList + '/list' + str(count),'w')
				li.append('list' + str(count))
				outfile.write(line)
				outfile.close()
				count = count + 1

		# MAC - for pipeline 3.1 all methods for gene step were replaced by this one
		try:
			for f in os.listdir(PathtoFiles+ '/FileLists_' + testPipelineList):			
				newPipe = Pipeline(PathtoFiles +'/'+ testPipelineList, PathtoFiles, ('geneStep',f),paramList,taxa2analyze,taxa2SF,majorClades,run_until_stage)
		
			### By inactivating the next line you can have access to all intermediary files
			Utilities.cleaner(testPipelineList, PathtoFiles, PathtoOutput)
		except Exception as e:
			elog = open('errorlog','a')
			line = open(PathtoFiles+ '/FileLists_' + testPipelineList + '/' + f,'r').read()
			elog.write(line + " failed on %s with: %s" % (f, e))
			elog.close()
			print ("failed on %s with: %s" % (f, e))
			
	#ACL -- Post-tree-building analyses and contamination loop
			
	if(resumer_param.lower() != 'na'):
		treeFolder = resumer_param
	else:
		treeFolder = '../' + testPipelineList + '_results2keep/'
	
	if(os.path.isdir(treeFolder)):
		if(use_contloop == 'y'):
			if(not os.path.isdir('../' + 'ContLoopOut')):
				os.system('mkdir ../' + 'ContLoopOut')
			ContLoopOut = '../ContLoopOut'
			restart = 'y'
			run = 0

			global cont_rules_file
			if "/" not in cont_rules_file:
				cont_rules_file = PathtoFiles + "/" + cont_rules_file
			else:
				cont_rules_file = cont_rules_file.replace(" ", "")
				
			seqs2remove_all = open(PathtoFiles + 'seqs2remove_all', 'w')
			
			taxaDBfile = open(PathtoFiles + 'taxaDBpipeline4', 'r').readlines()
			taxaDB = [taxon.strip('\n') for taxon in taxaDBfile]
			homologDB = [taxon[7:] for taxon in taxaDB if taxon.split(',')[0] == 'om']
			
			while restart == 'y': 
				# Run contamination removal. Check contaminationRemoval() in Utilities for details.
				print("Running the contamination loop: removing contamination and non-homologs ...\n")
				
				Utilities.contaminationRemoval(PathtoFiles, treeFolder, cont_rules_file, homologDB, remove_nm, cont_seq_or_clade, cont_num_contams, cont_target_clade, cont_min_target_presence, cont_at_least_sisters_file, cont_at_least_sisters_num, cont_target_taxa_file, include_outgroups_contloop, cont_filter_differential_coverage, cont_coverage_diff_OM, cont_coverage_diff_abs, return_cladegrabbing_subtrees)
		
				print("run " + str(run) + " PhyloTOL-conCleaner: contamination and non-homologs removal done ...")
				
				# Inside the ContLoopOut folder create an exclusive folder for temporary files in current run and move temporary files
				runOuts_path = "../ContLoopOut/" + testPipelineList + "_ContRun" + str(run)
				if(not os.path.isdir(runOuts_path)):
					os.system("mkdir " + runOuts_path)		
				os.system('cp -r ' + treeFolder + " " + runOuts_path)
				if(cont_seq_or_clade == 'seq'):
					os.system('mv sisterReport ' +  runOuts_path)
				os.system('mv nonHomologs ' +  runOuts_path)
				os.system('mv seqs2remove_out ' +  runOuts_path)
				os.system('mv nonHomol_treeWnhom ' +  runOuts_path)
				os.system('mv seqs2remove_out_treesWcont ' +  runOuts_path)
				if(os.path.isdir('CladeGrabbingSubtrees')):
					os.system('mv CladeGrabbingSubtrees ' + runOuts_path)
		
				# From outputs of contaminationRemoval() take contaminant sequences of the current run and attach them for the the list of ALL contaminant sequences.
				seqs2remove = open(runOuts_path + "/seqs2remove_out").readlines()
		
				for seq in seqs2remove : 
					seqs2remove_all.write("%s" % seq)
					print("run " + str(run) + " contamination loop: copying sequnce to seqs2remove_all --> " + seq.replace("\n", ""))
				print("run " + str(run) + " contamination loop: contamination sequeces appended to final list (seqs2remove_all)")
		
				# From outputs of contaminationRemoval() take list of trees with contamination, and use it for running next iteration of contamination removal.
				treesWcont = open(runOuts_path + "/seqs2remove_out_treesWcont").readlines()
				treesWnhom = open(runOuts_path + "/nonHomol_treeWnhom").readlines()
				os.system('mv ' + PathtoFiles + testPipelineList + " " + runOuts_path)
				newListOGs = open(PathtoFiles + testPipelineList, 'w')
		
				trees2reprocess = [val.strip() for val in list(set(treesWcont + treesWnhom))]

				for t, tree in enumerate(trees2reprocess): 
					print("run " + str(run) + " contamination loop: This tree has contamination or non-homologs and needs to be re-processed --> " + tree.replace("\n", ""))
					if(t != len(trees2reprocess) - 1):
						newListOGs.write(tree + '\n')
					else:
						newListOGs.write(tree)
				newListOGs.close()
			
				print("run " + str(run) + " contamination loop: new list of OGs to run generated")
				#exit()
		
				# If there were no trees with contamination stop loop of contamination removal (change restart to "n")
				if trees2reprocess == []:
					print("run " + str(run) + " contamination loop: no more contamination or non-homologs, this is the last run")
					restart = 'n'
				else:
					# Run PhyloTOL-resumer
					run += 1
					print("run " + str(run) + " contamination loop: running PhyloTOL with contamination removal\n\n")
					run_resumer(resumer_param, 'unaligned', testPipelineList, guidanceIter, seqcutoff, colcutoff, rescutoff, contloop_mode = True, ogList = trees2reprocess)
					print("run " + str(run) + " contamination loop finished")

			final_paths = { }

			for i, iter_dir in enumerate(os.listdir(ContLoopOut)):
				if('ContRun' in iter_dir):
					for file in os.listdir(ContLoopOut + '/' + iter_dir + '/' + treeFolder.split('/')[-1]):
						final_paths.update({ file : ContLoopOut + '/' + iter_dir + '/' + treeFolder.split('/')[-1] })
						
			if(not os.path.isdir('../Trees_ContamRemoved')):
				os.mkdir('../Trees_ContamRemoved')
				
			for file in final_paths:
				os.system('cp ' + final_paths[file] + '/' + file + ' ../Trees_ContamRemoved/' + file)

			folder_final_trees = '../Trees_ContamRemoved/'
		else:
			folder_final_trees = treeFolder

		if(resumer_param.lower() == 'na' or use_contloop == 'y'):
			os.system('python walk_tree_contamination_single.py ' + folder_final_trees + ' sisterReport ' + cont_query_clades + ' ' + cont_sister_clades + ' ' + cont_branch_length_filter + ' ' + cont_break_up + ' ' + cont_single_sister_only)
		
		if(treeFolder.endswith('/')):
			treeFolder = treeFolder[:-1]
		os.system('python color_by_clade.py ' + treeFolder + ' ' + tree_font_size)

		if(resumer_param.lower() != 'na' and concatAlignment == 'y'):
			if(not os.path.isdir(resumer_param + '/forConcatenation')):
				os.mkdir(resumer_param + '/forConcatenation')
			remove_paralogs(resumer_param + '/', resumer_param + '/forConcatenation/')

	else:
		print('\nSomething went wrong locating the output trees... look over the folder topology to make sure that you have a folder in the same directory as the "Scripts" folder that ends with "_results2keep" and contains your trees. Once things are fixed, you might want to set the resumer to start at "trees".\n')
		quit()
	
		
main()
























