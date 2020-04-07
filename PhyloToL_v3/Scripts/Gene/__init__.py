import sys
import os,re
import time
import operator
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio import AlignIO
from Bio import Align
from Bio.Align import MultipleSeqAlignment
from Bio import Phylo
from Bio import Entrez
Entrez.email = 'jgrant@smith.edu'
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import dendropy

##############################################################

##############################################################
			
class Gene:
	def __init__(self,OG,Pipeline):
		self.PathtoTemp = Pipeline.PathtoTemp
		self.PathtoOutput = Pipeline.PathtoOutput
		self.PathtoOGFiles = Pipeline.PathtoOGFiles
		self.OG = OG.strip() #OGnumber that identifies the gene
		self.taxa2analyze = Pipeline.taxa2analyze
		self.seqCodes = {} #dictionary holding numeric seq codes (keys) and original sequence names (values)
		self.paralogDict = {}
		self.MClist = {}
		self.sequenceDict = {}
		self.seqtoDelete = []
		self.backDict = {}
		self.seqLenCompCutOff = Pipeline.seqLenCompCutOff
	
	def getAllSeqs(self):   # Important for looping - MAC - Here we may need to take a look to the elifs
		
		if self.taxa2analyze is not 'all':
			fileOG = (self.PathtoOGFiles + '/' + self.OG)
			fileOG = SeqIO.parse(fileOG, 'fasta')
			filteredOG = open(self.PathtoOutput + '/' + self.OG + '_filtered.fas', 'w')
			
			for seq in fileOG:
				if seq.id[:10] in self.taxa2analyze: 
						filteredOG.write('>' + seq.id + '\n' + str(seq.seq) + '\n')
			
			filteredOG.close()
			
		else:
			os.system('cp ' + self.PathtoOGFiles + '/' + self.OG + ' ' +  self.PathtoOutput + '/' + self.OG + '_filtered.fas')

		os.system ('ls ' + self.PathtoOutput + '/fasta2keep/' + self.OG + '* > ' + self.PathtoOutput + '/' + self.OG + '_list')
		seqsList = open(self.PathtoOutput + '/' + self.OG + '_list', 'r')
		seqsList = seqsList.readlines()
		if seqsList: 
			os.system('cat ' + self.PathtoOutput + '/' + self.OG + '_filtered.fas' + ' ' +  self.PathtoOutput + '/fasta2keep/' + self.OG + '* > ' + self.PathtoOutput + '/' + self.OG + '_all.fas')
		else:
			os.system('cp '+ self.PathtoOutput + '/' + self.OG + '_filtered.fas ' + self.PathtoOutput + '/' + self.OG + '_all.fas')

##############################################################################
#renames seqs so information doesn't get lost in truncation latter in Guidance
##############################################################################

	def getSeqCodes(self):
		infile = open(self.PathtoOutput + '/' + self.OG + '_all.fas','r')		
		outfile = open(self.PathtoOutput + '/Guidance/'  + self.OG + 'forGuidance.fas','w')		
		codeout = open(self.PathtoOutput + '/seqcodeFiles/' + self.OG + '_seqcodes.txt','w')
		
		i = 0
		for line in infile:
			if line[0] == '>':
				seqid = line[1:].strip()
				
				newseqid = '>' + str(i + 1) + '_' +seqid.split('_')[1] + '_' +  seqid.split('_')[2]
				i = i + 1
				outfile.write(newseqid + '\n')
				codeout.write(seqid + ':' + newseqid[1:].strip() + '\n')
				self.seqCodes[seqid] = newseqid[1:].strip()
				
			else:
				line = re.sub('J','X',line)
				line = re.sub('O','X',line)
				line = re.sub('Z','X',line)
				outfile.write(line)
				
		outfile.write('\n')
		infile.close()
		outfile.close()
		codeout.close()
			
#####################################################################
# Prepare pre-guidance file
#####################################################################
	def fixGuidFile(self): 
		'''fixes some things that causes problems in guidance'''
		inortho = open(self.PathtoOutput + '/Guidance/'+ self.OG + 'forGuidance.fas','r')
		inseq = SeqIO.parse(inortho,'fasta')
		outortho = open(self.PathtoOutput + '/Guidance/'+ self.OG + 'forGuidance.fas2','w')
		for seq in inseq:
			 newline = re.sub('U','X',str(seq.seq))
			 newline = newline.replace('*', 'X')   # MACR --- for pipeline V3, we saw some of these cases in new data
			 if (newline) > self.seqLenCompCutOff:
			 	outortho.write('>' + seq.id + '\n' + newline + '\n')
		outortho.close()	 
		os.system('mv  ' + self.PathtoOutput + '/Guidance/'+ self.OG + 'forGuidance.fas2 ' + self.PathtoOutput + '/Guidance/'+ self.OG + 'forGuidance.fas')
		

##############################################################
# paralog removal for alignment concatenation
##############################################################		

	def getseqsfromCodeFile(self):
		for mc in ['Op','Pl','EE','Ex','Sr','Ba','Za','Am','Me']:      ## IMPORTANT ----> THIS IS TEMPORARY FOR RUNNING OGs WITH JD's DATA. REMOVE 'Me' AFTER
			self.MClist[mc]  = []
		infile = open(self.PathtoOutput + '/seqcodeFiles/' + self.OG + '_seqcodes.txt','r')
		for line in infile:
			mc = line.split('_')[0]
			ui = line.split('_')[0] + '_' + line.split('_')[1] + '_' +  line.split('_')[2]
			seqCode = line.split(':')[-1].split('_')[0] #just use number
			fullseq = line.split(':')[0].strip()

			self.MClist[mc].append(seqCode)

			self.sequenceDict[seqCode] = (mc,ui,fullseq)
			self.backDict[fullseq] = seqCode			
			self.paralogDict[ui] = [] #instantiate the dictionary as a dict of lists

		infile.close()		
	
	def removeParalogs(self, tokeepdir, forConcatenation_path):
		self.getseqsfromCodeFile()		
		self.uilist = []
		
		if os.path.exists(forConcatenation_path + 'RAxML_bestTree.' + self.OG + '_postguidance.fas'):
			self.tree_in = Phylo.read(forConcatenation_path + 'RAxML_bestTree.' + self.OG + '_postguidance.fas','newick')
			os.system('rm ' + forConcatenation_path + 'RAxML_bestTree.' + self.OG + '_postguidance.fas')   
		else: 
			print "cannot open" + forConcatenation_path + 'RAxML_bestTree.' + self.OG + '_postguidance.fas'
		  
		self.alignment = open(tokeepdir + self.OG + '_postguidance.fas.95gapTrimmed.fas_renamed.fas','r') 
		
		for seq in self.tree_in.get_terminals():
			try:
				ui = self.sequenceDict[str(seq).split('_')[0]][1] #ui is MC_mc_code
				self.paralogDict[ui].append(str(seq)) # so len is # of paralogs per taxon
				if ui not in self.uilist:
					self.uilist.append(ui)
			except:
				pass

		for ui in self.uilist:	
			if len(self.paralogDict[ui]) > 1:
				self.pickParalog(ui)
		self.deleteSeqsFromAlignment(forConcatenation_path)
		self.alignment.close()
	def deleteSeqsFromAlignment(self, forConcatenation_path):  #self.shortCode[seqCode] = (MC,ui,fullseq)
		nametoDelete = []
	
		infile2 = SeqIO.parse(self.alignment ,'fasta')
		outfile = open(forConcatenation_path + self.OG + '_paralog_purged.fas','w')
	
		for code in self.seqtoDelete:
			nametoDelete.append(self.sequenceDict[code][2])
		outfile2 = open(self.PathtoOutput + '/' + self.OG + '_sequencesKept.txt','w')
		for seq in infile2:
			if seq.id  not in nametoDelete: #if not in too delete, keep.
				newname = seq.id.split('_')[0] + '_' + seq.id.split('_')[1] + '_' + seq.id.split('_')[2]
				outfile.write('>' + newname + '\n' + str(seq.seq) + '\n')
				
				outfile2.write(newname + ':' + seq.id + '\n')
		outfile2.close()
	
		outfile.close
	def pickParalog(self, ui): #ui = unique identifier of taxa that have paralogs
		for sequence in self.paralogDict[ui]:
			clade = self.get_clade(sequence)
			try:
				seq_clade= self.get_parent(self.tree_in,clade)
				if seq_clade == None:
					seq_clade = clade
			except:
				seq_clade = clade
			
			if seq_clade != None:
				self.is_monophyletic(seq_clade,0,"None",0,sequence,[],[])
			else:
				self.seqtoDelete.append(sequence[0:9])
				self.seqtoDelete.append(sequence[0:10 ])
			self.paralogDict[ui] = []
		self.sort_clade_sizes()
	
	def sort_clade_sizes(self):
	#############################################################	
	#for each taxon, makes a list of seq/num pairs sorted py num
	#and passes it to get_seq_to_keep which will take the best seq
	#############################################################	
	
		inFile = open(self.PathtoTemp + '/clade_size.txt','r')
		infile = inFile.readlines()
		inFile.close()
		seqDict = {}
		seqnumDict = {}
		
		for line in infile:
			
			seq = line.split(':')[0]
			num = line.split(':')[1]
			seqDict[seq] = int(num)


		sorted_seqDict = sorted(seqDict.iteritems(), key=operator.itemgetter(1), reverse = True)	
		xx = open(self.PathtoTemp + '/clade_size.txt','w') #clear file
		xx.close()	
		self.get_seq_to_keep(sorted_seqDict) #call next method, not written yet

	def get_seq_to_keep(self,sorted_seqDict):
	#############################################################
	#checks to see if there is one sequence in a larger monophyletic clade
	#if not, makes a list of seqs and takes the one with the smallest 
	#branch length
	#############################################################
		check_dist_list=[]

		if len(sorted_seqDict) > 1:			# sorted_seqDict, list of (seqs,sizeofclade) sorted by size of clade
			if float(sorted_seqDict[0][1]) > float(sorted_seqDict[1][1]): #if the first is bigger than the next then it is the biggest
				sorted_seqDict.pop(0) #keep the first, get rid of the rest
				for double in sorted_seqDict:  
					self.seqtoDelete.append(double[0])  # put all the rest in seqtoDelete
			else: #More than seq in clade with equal numbers
				for double in sorted_seqDict:
					if float(double[1]) == float(sorted_seqDict[0][1]):
						check_dist_list.append(double) #get list of seqs to check for distance
					else:
						self.seqtoDelete.append(double[0])
				self.check_distance(check_dist_list)


	def check_distance(self,check_dist_list):
		distance_dict = {}
		for double in check_dist_list:
			for element in (self.tree_in.find_elements(name = double[0])):
				distance_dict[element] = self.tree_in.distance(element)
	
		sorted_distance_dict = sorted(distance_dict.iteritems(), key=operator.itemgetter(1))
		
		sorted_distance_dict.pop(0) #keep the shortest, remove others
		for element in sorted_distance_dict:
				self.seqtoDelete.append(str(element[0]).split('_')[0]) #want just the number, or SeqCode


	def is_monophyletic(self,clade,recursion_depth,child,length,sequence,mcCladeList,seqCladeList):
		#############################################################
		#checks to see if the clade is monophyletic for major clade and, if it is, 
		#recursively checks the parent clade until it finds a clade that is not monophyletic
		#############################################################
		out = open(self.PathtoTemp + '/clade_size.txt','a')
		try:
			for seq in clade.get_terminals(): #for each terminal sequence in the clade
		
				mc = self.sequenceDict[str(seq)][0] #ui is the major clade identifier i.e Op, Pl, etc
			
				if str(seq) not in seqCladeList:
					seqCladeList.append(str(seq)) #add the sequence to the seq list
					mcCladeList.append(mc) #add the unique identifier to the ui list
		
			parent = clade

			if len(mcCladeList) >= 1:
				if self.all_same(mcCladeList): #if the taxa are all the same in the clade, check to see if parent clade is also monophyletic
					length = len(mcCladeList)
					child = clade			
					parent = self.get_parent(self.tree_in,clade)			
					if parent != 'None':			
						recursion_depth = recursion_depth + 1				
						self.is_monophyletic(parent,recursion_depth,child,length,sequence,mcCladeList,seqCladeList)
				else:
					
					if recursion_depth < 1:
						out.write(sequence.strip() + ':' + str(1) + '\n')

					else:	
						out.write(sequence.strip() + ':' + str(length) + '\n')			

		except:
			out.write(sequence.strip() + ':' + str(0) + '\n')
		out.close()
						
	def all_same(self,uiList):
		return all(x == uiList[0] for x in uiList)	# returns True if all items are the same in the list
	
	def get_clade(self,sequence):
		cladelist = [] #leaf.names.
		for clade in self.tree_in.get_terminals():		
			if str(clade.name).strip() == str(sequence).strip():
				cladelist.append(clade)
		try:
			return cladelist[0]
		except:
			pass

	def get_parent(self,tree, child_clade):
		node_path = tree.get_path(child_clade)
		try:
			return node_path[-2]				
		except:
			return 'None' 