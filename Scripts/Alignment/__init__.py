from Bio import SeqIO
import time

class Alignment: #reads in an alignment file and makes an Alignment instance

	def __init__(self,f):
		self.name = f
		self.numchar = 0
		self.numseq = 0
 		self.seqs = []
 		self.model = ""
 		self.idlist = []
		
		#self.populate()
		
		
	def populate(self):
		#print self.name
		infile = open(self.name,'r')
		inSeq = SeqIO.parse(infile,'fasta')
		print infile
		
		for seq in inSeq:
				self.seqs.append(seq)
				self.idlist.append(seq.id)
#			else:
#				print 'duplicate seqs in ' + self.name
				
		self.numchar = len(self.seqs[0].seq)
		self.numseq = len(self.seqs)
		time.sleep(10)
		#print self.numseq



	def setmodel(self,model):
		self.model = model
	
	
	def write(self, name):
		duplicates = []
		outfile = open(name + '_alignment.fas','a')
		for seq in self.seqs:
			if str(seq) not in duplicates:
				duplicates.append(str(seq))
				outfile.write('>' + seq.id + '\n' + str(seq.seq) + '\n')