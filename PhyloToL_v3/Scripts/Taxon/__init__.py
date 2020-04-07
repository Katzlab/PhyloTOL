import os,re,gc
import subprocess
import Utilities

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from Bio import AlignIO 


##############################################################

class Taxon: 
	def __init__(self, fileName, Pipeline, restart):
		self.PathtoSeqFiles = Pipeline.PathtoSeqFiles
		self.PathtoBlastFiles = Pipeline.PathtoBlastFiles
		self.PathtoOutput = Pipeline.PathtoOutput
		self.OGList  = Pipeline.OGList 
		self.taxa2SF = Pipeline.taxa2SF
		self.fileName = fileName
		self.seqs = [] #a list of seqeunces from the taxon
		self.blast = [] #a list of blast records		
		self.code = "" #taxon code (get from file?
		self.OGlist = [] # list of Ogs of interest in the taxon--'of interest' depends on the pipeline
		self.OGxSeqHashProt = {} #dictionary with translated seqeunces for each og of interest are stored
		self.blastCutOff = Pipeline.blastCutOff
		self.tooSimCutOff = Pipeline.tooSimCutOff

		##############################################################
		# Step 0 - This is where the different steps, coded further down, are called
		##############################################################
		
		if restart == 'no':
			self.set_seqs(self.PathtoSeqFiles,fileName)
			self.set_code()
			
			self.set_blast(fileName, self.PathtoBlastFiles)
			self.parseblast(fileName)
			outfile = open('log_'+ self.code,'a') #****
			if any(a != [] for a in self.OGxSeqHashProt.values()) == True:    #  or self.OGxSeqHashProt != []:
				self.of_sf(fileName) 
			outfile.write(str(self.OGxSeqHashProt) + '\n')
			for i in range(len(self.blast)): #tried just setting self.blast = [] but some lists toolong? this removes a big list freeing memory
				self.blast.pop()	
			outfile.close()

			if os.path.isdir(self.PathtoOutput + '/intermediate_and_logfiles/'):
				os.system('mv log_' + self.code  +  ' ' + self.PathtoOutput + '/intermediate_and_logfiles/')
			gc.collect()
			
		else:
			infas = SeqIO.parse(open(self.PathtoOutput + '/fasta2keep/' + filename + 'fastatokeep.fas','r'),'fasta')
			if infas != []:
				self.code = infas[0].id.split('_')[0] + '_' + infas[0].id.split('_')[1]
        
##############################################################
#  methods
##############################################################
	
	def clearmem(self):
		#try to free up memory...			
		self.seqs = [] #a list of seqeunces from the taxon
		self.blast = [] #a list of blast records		
		self.code = "" #taxon code (get from file?
		self.OGlist = [] # #list of Ogs of interest in the taxon--'of interest' depends on the pipeline
		self.OGxSeqHashProt = {} #dictionary with translated seqeunces for each og of interest are stored
		
	def set_seqs(self,PathtoSeqFiles, filename):
		print filename[:10] + ": Parsing sequences" 
		infile = SeqIO.parse(open(PathtoSeqFiles + filename,'r'),'fasta')
		for seq in infile:
			self.seqs.append(seq)	
	
	def set_code(self):
		self.code = self.seqs[0].id.split('_')[0] + '_' + self.seqs[0].id.split('_')[1]+ '_' + self.seqs[0].id.split('_')[2]		
			
	def set_blast(self,fileName, PathtoBlastFiles):	
		print fileName[:10] + ': Parsing blast report'
		test = ''
		for blastFileName in os.listdir(PathtoBlastFiles):
			if re.search(fileName[:10], blastFileName):
				blastrecords = NCBIXML.parse(open(PathtoBlastFiles + blastFileName,'r'))		
				for blast in blastrecords:
					if blast.descriptions:
						self.blast.append(blast)
		
		if test == '':
			outfile = open('log_' + self.code,'a')
			outfile.write('no BLAST file found for ' + self.code)
			outfile.close()
			
		
###################################################################################################	
#	parseblast, getgrpp and makeFasta								  							  #
#	Looks through the blasts and makes files for the taxon for each OG							  #
###################################################################################################	

	def parseblast(self, fileName): 
		print fileName[:10] + ': Collecting blast results that match OGs of interest' 
		interestList = self.OGList 
		
		for OGGroup in interestList:
			self.OGxSeqHashProt[(OGGroup,self.code)] = []

		for record in self.blast:
			if record.descriptions:

				if record.alignments[0].hsps[0].expect < self.blastCutOff:

					outfile = open('log_'+ self.code,'a') #****
					outfile.write('e is ' + str(record.alignments[0].hsps[0].expect) + ' for ' + self.code + ' cot off is ' + str(self.blastCutOff) + '\n')
					outfile.close()
					self.getgrpp(record,interestList) 					

				else:
					self.blast.remove(record)
					outfile = open('log_'+ self.code,'a') #****
					outfile.write('e is ' + str(record.alignments[0].hsps[0].expect) + ' for ' + self.code + ' below cutoff' + str(self.blastCutOff) + '\n')
					outfile.close()
			else:
				self.blast.remove(record)	
	
	def getgrpp(self, record,interestList):
	
		if record.descriptions:
			try:
				OGGroup = 'OG5_' +  record.alignments[0].hit_def.split('_')[-1].strip()
		
			except:
				OGGroup = record.alignments[0].hit_def

			if OGGroup in interestList:
				self.makeFasta(record.query,OGGroup)				

	# MACR -- Here we produce the fasta files OG_taxon containing all sequences that are going to be filtered using iterUblast.py. 

	def makeFasta(self, query, OGGroup):	#also populates OGxSeqHash
		outname = OGGroup + '_' + self.code  + '_out.txt'
		self.OGlist.append(OGGroup)
		outpathname = self.PathtoOutput + '/' + outname
		
		for seq in self.seqs:
			if seq.name == query:
				self.OGxSeqHashProt[(OGGroup,self.code)].append(seq)	
				out1 = open(outpathname,'a')	
				out1.write('>' +  seq.name + '\n' + str(seq.seq)  + '\n')
				out1.close()		
###################################################################################################	

			
###################################################################################################	
#	overlap filter (OF) and similarity filter (SF)
###################################################################################################
	
	def of_sf(self, fileName):
		if not os.path.exists(self.PathtoOutput + '/fasta2keep'):
			os.system('mkdir ' + self.PathtoOutput + '/fasta2keep')
			
		for OG in self.OGList:
			filename = OG + '_' + self.code			
			
			file4ublast = self.PathtoOutput + '/' + filename + '_out.txt'
			
			if os.path.exists(file4ublast):
				Utilities.iterUblast(file4ublast, self.tooSimCutOff, self.taxa2SF, self.code)
				
			outfile = open('log_'+ self.code,'a') #****
			outfile.close()	
