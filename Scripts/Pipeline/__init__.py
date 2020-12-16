import os
import re
import Utilities
from Taxon import Taxon
from Gene import Gene
from datetime import date
from Bio import SeqIO

##############################################################

class Pipeline: #given a file with a list of OGs, makes the interestList and
					#a fasta file for the OGs.  Also a series of paths for output
					#and a list of Taxon objects
	
	def __init__(self,OGofInterestFile, PathtoFiles, restart, paramList,taxa2analyze,taxa2SF,majorClades,mode):
		self.blastCutOff,self.seqLenCompCutOff,self.tooSimCutOff,self.guidanceIter,self.seqcutoff,self.colcutoff, self.rescutoff, self.concatAlignment = paramList
		testlist = OGofInterestFile.strip('/').split('/')[-1]
		self.startdate = str(date.today())
		self.PathtoFiles = PathtoFiles
		self.OGofInterestFile = OGofInterestFile
		self.majorClades = majorClades
		self.taxa2SF = taxa2SF
		self.taxa2analyze = taxa2analyze
		self.OGList = [] 
		self.OGsofInterestList = [] #list of Gene instances
		self.TaxonList = [] # list of Taxon instances
		self.OGswithTaxofInterest = []
		self.mode = mode
		#make directories and set paths:
		PATH =  '../my-data/' + OGofInterestFile.split('/')[-1] + '_results/'
		if not os.path.exists(PATH):
			os.system('mkdir ' +  PATH)
		if not os.path.exists(PATH  + '/Temp'):
			os.system('mkdir ' + PATH  + '/Temp')
		if not os.path.exists(PATH  + '/Output'):
			os.system('mkdir ' + PATH  + '/Output')	
		self.PathtoTemp = PATH + 'Temp'
		self.PathtoOutput = PATH + 'Output'
		self.PathtoSeqFiles = PathtoFiles + '/ncbiFiles/'
		self.PathtoBlastFiles = PathtoFiles + '/BlastFiles/'
		self.PathtoOGFiles = PathtoFiles + '/allOG5Files/'
		self.PathtoPartialFiles = PathtoFiles + '/FileLists_' + testlist + '/'
		if not os.path.exists(self.PathtoOutput + '/' + 'fasta2keep'):
			os.system('mkdir ' + self.PathtoOutput + '/' + 'fasta2keep')
		if not os.path.exists(self.PathtoOutput + '/intermediate_and_logfiles/'):
			os.system('mkdir ' + self.PathtoOutput + '/intermediate_and_logfiles/')
		if not os.path.isdir(self.PathtoOutput + '/seqcodeFiles/'):
			os.system('mkdir ' + self.PathtoOutput + '/seqcodeFiles/')
		
# Calling methods ----- MAC

		if restart == 'no':
			self.setInterestList()
			self.writelog('Pipeline:' + OGofInterestFile + ',' + PathtoFiles)
			self.makedb()
		else:
			if restart[0] == 'queueTaxa':
				self.setInterestList()
				self.queuetax(restart[1])
				return	
			elif restart[0] == 'geneStep': ##NEW
				self.restart_geneStep(self.PathtoOutput,OGofInterestFile,PathtoFiles, restart[1], self.guidanceIter, self.seqcutoff,self.colcutoff,self.rescutoff,self.concatAlignment) ##NEW
			else:
				self.setInterestList()
				self.restart(self.PathtoOutput,OGofInterestFile,PathtoFiles)

	def setInterestList(self):	
		
	 	# JG - to allow running of different Pipelines, output files names have to be different
		#get the list of OGs for this pipeline and build gene instances
		
		for line in open(self.OGofInterestFile,'r'):
			self.OGList.append(line.strip())
			NewGene = Gene(line,self)
			self.OGsofInterestList.append(NewGene)

	def writelog(self,string):
		self.logfile = open(self.PathtoOutput + '/logfile', 'a') #this will hold info as the pipeline progresses, hopefully to help restart
		self.logfile.write(string + '\n')
	
##############################################################
##NEW stuff just up to guidance for Mario, 3/14/15
##############################################################

	def restart_geneStep(self,PathtoOutput,OGofInterestFile,PathtoFiles,FileListOG, guidanceIter, seqcutoff, colcutoff, rescutoff, concatAlignment,):
		OGswithTaxofInterestfile = open(self.PathtoPartialFiles + '/' + FileListOG,'r')
		for line in OGswithTaxofInterestfile:
			self.OGList.append(line.strip())
			NewGene = Gene(line,self)
			self.OGsofInterestList.append(NewGene)
			for gene in self.OGsofInterestList:

				gene.getAllSeqs()
				if not os.path.exists(self.PathtoOutput + '/Guidance/') : os.system('mkdir ' + self.PathtoOutput + '/Guidance/') 
				gene.getSeqCodes()
				gene.fixGuidFile()
				
				# MACR ---- Here are some changes for Pipeline3 ---------------	
				OGofInterestFile = self.OGofInterestFile.split('/')[-1]	
				results_guidance = Utilities.iterGuidance(OGofInterestFile, line.strip(), self.PathtoOutput, self.guidanceIter, self.seqcutoff, self.colcutoff, self.rescutoff, self.mode)
				
				if results_guidance[1] == 'n':
					if self.concatAlignment == 'y':
						if self.mode != "nr":
							forConcatenation_path = results_guidance[0] + 'forConcatenation/'
							if not os.path.exists(forConcatenation_path):
								os.system('mkdir ' + forConcatenation_path)
							os.system('cp ' + results_guidance[0] + 'RAxML_bestTree.' + line.strip() + '_postguidance.fas ' + forConcatenation_path)
							codesDic = Utilities.build_codesDic(line.strip(), self.PathtoOutput)
							Utilities.renamefiles(line.strip(), results_guidance[0], codesDic)
							gene.removeParalogs(results_guidance[0], forConcatenation_path)
						else:
							codesDic = Utilities.build_codesDic(line.strip(), self.PathtoOutput)
							Utilities.renamefiles(line.strip(), results_guidance[0], codesDic)
					else:
						codesDic = Utilities.build_codesDic(line.strip(), self.PathtoOutput)
						Utilities.renamefiles(line.strip(), results_guidance[0], codesDic)
				else:
					codesDic = Utilities.build_codesDic(line.strip(), self.PathtoOutput)
					if codesDic:
						Utilities.renamefiles(line.strip(), results_guidance[0], codesDic)

##############################################################
#queue taxa
##############################################################		
	def queuetax(self, taxfile):
		#read the taxon files and build Taxon Instances
		flag = 0
		for file in os.listdir(self.PathtoSeqFiles):
			if file == taxfile:
				flag = 1
				NewTaxon = Taxon(file, self, 'no')
				self.TaxonList.append(NewTaxon)
				self.writelog('Taxon:' + NewTaxon.code)
				NewTaxon.clearmem()
		if flag == 0: #never finds it...
			self.writelog('Taxon: Cant find ' + NewTaxon.code)
