import os
import sys
import re
import time

""""
iterUblast() is a function to run JD's script for filtering sequences based on similarity and using 'Cluster'. It runs from the Gene class
and uses the files OG_taxa as input and the cutoff set by the user in pipeline_parameters.txt
"""

def iterUblast(file4ublast, ubCutoff, taxa2SF, taxon):
	
	if taxon in taxa2SF:
		os.system('python iterUblast.py ' + file4ublast + ' ' + str(ubCutoff))
	else:
		os.system('python iterUblast.py ' + file4ublast + ' n')
	
"""
iterGuidance() is a function that calls Miguel's script too iterate guidance. Guiven that
Miguel's script also runs raxml, this function would end up performing both steps
guidance and raxml. It needs the name of the lists of OGs, the OG that is currently
being processed the Gene step (class Gene), the path to 'my-data' and the sequence and column
cutoffs for guidance. It creates a series of folders to preform guidance and raxml and puts 
the most important files from results to a folder called 'results2keep'. This function is being 
called from the class Gene (function __init__.py) and it retrieves the path to 'results2keep' 
"""

def iterGuidance(oglist, og, PathtoOutput, guidanceIter, seqcutoff, colcutoff, rescutoff, mode):

	PathtoOutput = os.path.abspath(PathtoOutput)
	tokeepdir = PathtoOutput + '/' + oglist + '_results2keep/'
	results = PathtoOutput + '/' + oglist + '_postguidance' +  '/'
	outdir = results + 'out/'
	tempdir = results + 'temp/'
	

	if not os.path.exists(results):
		os.system("mkdir " + results)
	if not os.path.exists(outdir):
		os.system("mkdir " + outdir)
	if not os.path.exists(tempdir):			
		os.system("mkdir " + tempdir)
	if not os.path.exists(tokeepdir):
		os.system("mkdir " + tokeepdir)
	
	if 'OG' in og:
		path2og = PathtoOutput + '/Guidance/'+ og + 'forGuidance.fas'
				
		if str(mode) != "ng":
			print("\n" + og + ": Running Guidance")
			print('bash exeGuidance_2.2.sh -i ' + path2og + ' -o ' + tempdir + ' -t 20 -c ' + outdir + ' -g ' + str(guidanceIter) + ' -s ' + str(seqcutoff) + ' -l ' + str(colcutoff) + ' -r ' + str(rescutoff) + ' -m ' + str(mode))
			os.system('bash exeGuidance_2.2.sh -i ' + path2og + ' -o ' + tempdir + ' -t 20 -c ' + outdir + ' -g ' + str(guidanceIter) + ' -s ' + str(seqcutoff) + ' -l ' + str(colcutoff) + ' -r ' + str(rescutoff) + ' -m ' + str(mode))
			os.system('cp ' + path2og + ' ' + tokeepdir + og + '_preguidance.fas') # non-aligned pre-guidance file		
		
#			if not os.path.exists(outdir + og + 'forGuidance.fas.output/NInitialSeqBelow4') or not os.path.exists(outdir + og + 'forGuidance.fas.output/NAboveCutoffBelow4'):	
			if not os.path.exists(outdir + og + 'forGuidance.fas.output/NSeqBelow4'):	
				guidanceFail = 'n'
				os.system('mv ' + outdir + og + 'forGuidance.fas.output/' + og + 'forGuidance.fas.95gapTrimmed.fas ' + tokeepdir + og + '_postguidance.fas.95gapTrimmed.fas') # trimmed Post-guidance alignment in fasta format
				os.system('mv ' + outdir + og + 'forGuidance.fas.output/' + og + 'forGuidance.fas.Post ' + tokeepdir) # Post-guidance alignment in phylip format
#				os.system('mv ' + outdir + og + 'forGuidance.fas.output/' + og + 'forGuidance.fas ' + tokeepdir) # non-aligned file before last guidance run
				if str(mode) != "nr":
					os.system('mv ' + outdir + og + 'forGuidance.fas.output/' + 'RAxML_bestTree.' + og + 'forGuidance.fas ' + tokeepdir + 'RAxML_bestTree.' + og + '_postguidance.fas')
			
				print("convert postguidance alignment from .phy to .fas")
				postGuid2Fasta = open(tokeepdir + og + 'forGuidance.fas.Post', 'r').readlines()
				postGuidFasta = open(tokeepdir + og + '_postguidance.fas', 'a') # new postguidance alignment in fasta
				os.system('rm ' + tokeepdir + og + 'forGuidance.fas.Post')
				
				for line in postGuid2Fasta:
					if "_" in line: # if line contains a taxon
						line = re.sub(" +", " ", line) # replacing variable number of spaces for a space between OTU and aligned seq in phylip
						taxon = line.split(" ")[0]
						seq = line.split(" ")[1]
						postGuidFasta.write(">%s\n" % taxon)
						postGuidFasta.write("%s\n" % seq)
			else:
				guidanceFail = 'y'
		else:
			os.system('cp ' + path2og + ' ' + tokeepdir + og + '_preguidance.fas') # non-aligned pre-guidance file
			guidanceFail = 'y'
			
	return tokeepdir,guidanceFail

	
"""
The function build_codesDic() takes the file OG_seqcodes.txt that contains the key for the names
(e.g., uniqueNumber_mc_Spcs : MC_mc_Spcs_SeqID) and uses it for biulding a dictionary. This function 
is being called from class Gene (function __init__.py) and returns the dictionary, which is used latter 
by the function renamefiles().  
"""

def build_codesDic(og, PathtoOutput):
	if os.path.exists(PathtoOutput + '/seqcodeFiles/' + og + '_seqcodes.txt'):
		codes = open (PathtoOutput + '/seqcodeFiles/' + og + '_seqcodes.txt')
		codes = codes.readlines()
		codesDic = {}

		for code in codes:
			code = code.strip("\n")
			taxon = code.split(":")[0]
			number = code.split(":")[1]
			codesDic[number] = taxon
	else:	
		codesDic = ""	

	return codesDic
	
"""
The function renamefiles() takes the dictionary that contains key of names
(e.g., uniqueNumber_mc_Spcs : MC_mc_Spcs_SeqID) and uses it for renaming all
files from folder results2keep. This function is being called after calling 
build_codesDic() in class Gene (function __init__.py).
"""

def renamefiles(og, tokeepdir, codesDic):

	print("\n" + og + ": Renaming files:")
	for file2rename in os.listdir(tokeepdir):
		if (og in file2rename):
			if ('renamed' not in file2rename):
				if not 'Tree' in file2rename:
					to_rename = open(tokeepdir + file2rename, 'r')
					to_rename = to_rename.readlines()
					file_renamed = open(tokeepdir + file2rename + '_renamed.fas', 'a')
					for renameline in to_rename:
						renameline = renameline.strip('\n')
						if '>' in renameline:
							code = renameline.replace('>', '')
							renameline = renameline.replace(code, codesDic[code])
							file_renamed.write("%s\n" % renameline)
						else:
							file_renamed.write("%s\n" % renameline)
				else:
					tree2r = []
					tree2rename = open(tokeepdir + file2rename, 'r')
					tree2rename = tree2rename.readline()
					tree2rename = tree2rename.split(',')

					for leave in tree2rename:
						taxon = leave.split(":")[0]
						taxon = taxon.replace("(", "")

						for code in codesDic:
							if code == taxon:
								leave = leave.replace(code, codesDic[code])
								tree2r.append(leave)
					
					treeRenamed_line = ",".join(tree2r)
					treeRenamed = open(tokeepdir + file2rename + '_renamed.tre', 'w')
					treeRenamed.write("%s" % treeRenamed_line)
				
				print(file2rename +" : renamed")
				
				os.system('rm ' + tokeepdir + file2rename) # everything that IS NOT a renamed tree or a renamed fasta file is removed
			
"""
The function cleaner() removes all files that you don't need. It is called 
from the script PhyloTOL.py
"""

def cleaner(testPipelineList, PathtoFiles, PathtoOutput):
	
	os.system('rm -f *log*')
	os.system('rm -f *RAxML_*')
	os.system('rm -f Utilities.pyc')
	os.system('rm -f Gene/__init__.pyc')
	os.system('rm -f Taxon/__init__.pyc')
	os.system('rm -f Pipeline/__init__.pyc')
	os.system('rm -rf __pycache__')
	os.system('rm -rf Gene/__pycache__')
	os.system('rm -rf Taxon/__pycache__')
	os.system('rm -rf Pipeline/__pycache__')
	
	if os.path.exists(PathtoOutput + '/' + testPipelineList + '_results2keep/'):
		os.system('mv ' + PathtoOutput + '/' + testPipelineList + '_results2keep/ ' + '../')

	os.system('rm -r ' + "../my-data/")
	os.system('rm -r ' + PathtoFiles + 'FileLists_' + testPipelineList + '/')

"""
The function contaminationRemoval removes contamination from sequence data (ncbiFiles). It is called 
from the script PhyloTOL_Ccleaner.py.

IMPORTANT :
- Databases are only useful for current experiment. They have to be replaced because sequences removed by
GUIDANCE (non-homologs by context) are removed in every iteration
- The ouput "seqs2remove_out" contains contamination sequences (bassed on rules provided by user). This 
file can be used for cleanning permanent databases. This file is re-written in each run, but PhyloTOL_conCleaner.py
concatenates the file of every run in the file seqs2remove_all
"""

def contaminationRemoval(treeFolder, PathtoFiles, rules, homologDB, mode):
			
	nonHomol_out = open("nonHomologs", "w")
	nonHomol_treeWnhom = open("nonHomol_treeWnhom", "w")
	
	for file in os.listdir(treeFolder):
		nonHomINtree = 'n'
		if file.endswith('postguidance.fas_renamed.fas'):
			og = file.replace("_postguidance.fas_renamed.fas", "")
			postseqs = open(treeFolder + file, "r").readlines()
			preseqs = open(treeFolder + og + '_preguidance.fas_renamed.fas', "r").readlines()
			for preseq in preseqs:
				if '>' in preseq:
					if preseq not in postseqs:
						preseq = preseq.replace(">", "")
						if preseq[0:10] not in homologDB:
							nonHomINtree = 'y'
							nonHomol_out.write("%s" % preseq)
							print("Utilities.py - contaminationRemoval(): this sequences is not homolog --> " + preseq.replace("\n", ""))
						else:
							print("Utilities.py - contaminationRemoval(): this sequences is not homolog --> " + preseq.replace("\n", ""))
			if nonHomINtree == 'y' : nonHomol_treeWnhom.write("%s\n" % og)
			
	nonHomol_out.close()
	nonHomol_treeWnhom.close() 
			
	os.system('python3 walk_tree_contamination_single.py ' + treeFolder + ' sisterReport')
	print("\nUtilities.py - contaminationRemoval(): Writing sister report ...\n")
	time.sleep(60)
	os.system('ruby seqs2remove.rb ' + PathtoFiles + ' sisterReport ' + rules + ' seqs2remove_out ' + 'nonHomologs ' + mode)
	print("\nUtilities.py - contaminationRemoval(): Replacing fasta files ...\n")
	time.sleep(60)