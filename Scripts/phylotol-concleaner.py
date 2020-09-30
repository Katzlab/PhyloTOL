import os, sys
import Utilities

script, rules, mode = sys.argv

"""
running as ...
python phylotol-concleaner.py path_to_rules mode

modes:
c1 = cleaning contamination in gene families DB
c2 = cleaning contamination in added taxa DB
c3 = cleaning contamination in gene families DB and added taxa DB
"""

print("\n\n ** runing PhyloTOL-conCleaner **\n\n")

# read parameters file and take the path to fasta files and the name of the OG list. 
infile = open('pipeline_parameter_file.txt','r').readlines()
for line in infile:
	if line[0] == '#':
		attribute = line.split()[0].strip('#')
		value = line.split()[2].strip()
		if attribute == 'PathtoFiles':
			PathtoFiles = value + "/"
		elif attribute == 'testPipelineList':
			testPipelineList = value


# Start variables for looping contamination removal, including creating the folder for temporal files and creating the file that will contain ALL sequences of contamination
os.system('mkdir ../' + 'Temp')
temp = '../Temp'
restart = 'y'
run = 0

if "/" not in rules:
	rules = PathtoFiles + "/" + rules
else:
	rules = rules.replace(" ", "")

seqs2remove_all = open(PathtoFiles + 'seqs2remove_all', 'w')
treeFolder = '../' + testPipelineList + '_results2keep/'
if os.path.exists(treeFolder) : print('terminating PhyloTOL-conCleaner: the folder ' + '../' + testPipelineList + '_results2keep exists. You should start with clean folders\n\n')


# grabbing taxa included in the homologs database
taxaDBfile = open(PathtoFiles + 'taxaDBpipeline3', 'r').readlines()
taxaDB = [taxon.strip('\n') for taxon in taxaDBfile]
homologDB = [taxon[7:] for taxon in taxaDB if taxon.split(',')[0] == 'om']

# making sure that a valid mode was chosen... Valid modes are c1, c2 or c3

if mode not in ["c1", "c2", "c3"]:
	print("prease choose a correct mode of cleaning contamination:\nc1: in gene families DB\nc2: in added taxa DB\nc3: in both")
	quit()
else:
	print("mode of contamination cleaning: %s\n" % mode)

# Start looping ... The loop will only stop when restart changes to 'n'
while restart == 'y': 
	
	if os.path.exists(treeFolder):
	
		# See if there are results from PhyloTOL. If so, do contamination removal. If there are not results, check 'else'.

		
		# Run contamination removal. Check contaminationRemoval() in Utilities for details.  
		print("run " + str(run) + " PhyloTOL-conCleaner: removing contamination and non-homologs ...\n")
		Utilities.contaminationRemoval(treeFolder, PathtoFiles, rules, homologDB, mode)
		print("run " + str(run) + " PhyloTOL-conCleaner: contamination and non-homologs removal done ...")

		# Inside the Temp folder create an exclusive folder for temporal files in current run and move temporary files
		runOuts_path = "../Temp/" + testPipelineList + "_ContRun" + str(run) 
		os.system("mkdir " + runOuts_path)		
		os.system('mv ' + treeFolder + " " + runOuts_path)
		os.system('mv sisterReport ' +  runOuts_path)
		os.system('mv nonHomologs ' +  runOuts_path)
		os.system('mv seqs2remove_out ' +  runOuts_path)
		os.system('mv nonHomol_treeWnhom ' +  runOuts_path)
		os.system('mv seqs2remove_out_treesWcont ' +  runOuts_path)
		
		# From outputs of contaminationRemoval() take contaminant sequences of the current run and attach them for the the list of ALL contaminant sequences.
		seqs2remove = open(runOuts_path + "/seqs2remove_out").readlines()
		
		for seq in seqs2remove : 
			seqs2remove_all.write("%s" % seq)
			print("run " + str(run) + " PhyloTOL-conCleaner: copying sequnce to seqs2remove_all --> " + seq.replace("\n", ""))
		print("run " + str(run) + " PhyloTOL-conCleaner: contamination sequeces appended to final list (seqs2remove_all)")
		
		# From outputs of contaminationRemoval() take list of trees with contamination, and use it for running next iteration of contamination removal.
		treesWcont = open(runOuts_path + "/seqs2remove_out_treesWcont").readlines()
		treesWnhom = open(runOuts_path + "/nonHomol_treeWnhom").readlines()
		os.system('mv ' + PathtoFiles + testPipelineList + " " + runOuts_path)
		newListOGs = open(PathtoFiles + testPipelineList, 'w')
		
		trees2reprocess = list(set(treesWcont + treesWnhom))
		
#		for tree in treesWcont: 
#			print "run " + str(run) + " PhyloTOL-conCleaner: This tree has contamination and needs to be re-processed --> " + tree.replace("\n", "")
#			newListOGs.write("%s" % tree)
#		newListOGs.close()

		for tree in trees2reprocess: 
			print("run " + str(run) + " PhyloTOL-conCleaner: This tree has contamination or non-homologs and needs to be re-processed --> " + tree.replace("\n", ""))
			newListOGs.write("%s" % tree)
		newListOGs.close()
			
		print("run " + str(run) + " PhyloTOL-conCleaner: new list of OGs to run generated")
		
		# If there  were no trees with contamination stop loop of contamination removal (change restart to "n")
#		if treesWcont == []:
		if trees2reprocess == []:
			print("run " + str(run) + " PhyloTOL-conCleaner: no more contamination or non-homologs, this is the last run")
			restart = 'n'
		else:
			# Run PhyloTOL
			run += 1
			print("run " + str(run) + " PhyloTOL-conCleaner: runing PhyloTOL with contamination removal\n\n")
			os.system("python phylotol.py " + mode)  # run the pipeline with mode 'ct'
			print("run " + str(run) + " PhyloTOL-conCleaner: PhyloTOL done")
	
	else:
		
		# There are not results from PhyloTOL, then runs it for first time.
		
		run += 1
		print("run " + str(run) + " PhyloTOL-conCleaner: runing PhyloTOL - pre-contamination removal \n\n")
		os.system("python phylotol.py " + mode)
		print("run " + str(run) + " PhyloTOL-conCleaner: PhyloTOL done")

seqs2remove_all.close()