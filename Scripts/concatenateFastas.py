import os,re,time
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Alignment import Alignment

numchar = 0
numseq = 0
seqidlist = []

def writelog(string):
	logfile = open('logfile.txt','a')
	logfile.write(string + '\n')
	logfile.close()


def write_error(string):
	logfile = open('errorfile.txt','a')
	logfile.write(string + '\n')
	logfile.close()		


def concat(newAlignment,f,model):
	print(f)
	#inSeq2 = SeqIO.parse(open('concatenated.fas','r'),'fasta')
	
	a = Alignment(f)	#make an alignment instance from the file f
	a.setmodel(model)	#set the model for the alignment
	a.populate()
	
	
	fixFastas(newAlignment, a)
	tempAlignment = Alignment('')
	
	writelog(f + ':' + str(a.numchar) + ' ' + model)
	for seq in newAlignment.seqs:
		flag = 0
		for aseq in a.seqs:
			
			if seq.id == aseq.id:
				flag = 1
				newseq = SeqRecord(seq = seq.seq + aseq.seq, id = seq.id)

				tempAlignment.seqs.append(newseq)
				tempAlignment.idlist.append(newseq.id)

	tempAlignment.numchar = len(tempAlignment.seqs[0].seq)
	tempAlignment.numseq = len(tempAlignment.seqs)	
	#print tempAlignment.numchar
	#print tempAlignment.numseq
	
	return tempAlignment

	
def fixFastas(newAlignment, a):

	for seqid in newAlignment.idlist: #for all sequences in ssu, take their seqID
		if seqid not in a.idlist:
			
			newseq = SeqRecord(seq = '-' * a.numchar, id = seqid)
			a.seqs.append(newseq)
			a.idlist.append(newseq.id)

	for seqid in a.idlist: #for all sequences in ssu, take their seqID
		if seqid not in newAlignment.idlist:
			
			newseq = SeqRecord(seq = '-' * newAlignment.numchar, id = seqid)
			newAlignment.seqs.append(newseq)
			newAlignment.idlist.append(newseq.id)
			
			
			
			
def makephy(f):
	inFile = []
	infile = open(f,'r')
	for line in infile:
		try:
			inFile.append(line.strip())
		except:
			print line
	name = f[:-4] + '.phy'
	outfile = open(name,'a')
	seqCount = 0	
	for line in inFile:
		if line[0] == '>':
			seqCount = seqCount + 1
		else:
			charCount = len(line.strip())
		outfile.write(str(seqCount) + ' ' + str(charCount) + '\n')
	for line in inFile:
		if line[0] == '>':	
			seqID = line.strip()[1:] + '         '
			seqID2 = seqID[0:10] + ' '
			outfile.write('\n' + seqID2)
		else:
			outfile.write(line.strip())
	

def main():
	print('##############################################################')
	print('This script must be in a directory with the folder called "Alignment" and the fastas you want to concatenate.')
	print('The fastas must have the extension ".fasta" or it wont work.')
	print('The output files are the concatenated seqs as a fasta and as a phylip for raxml-ing and a logfile.')
	print('The log file tells the length of the added sequence and the order in which they are added, so you can figure out your breakpoints. ')
	print('##############################################################')
	 
	path2alignments = (raw_input("path to the purged alignments: ")).strip() + "/"
	outfile = open('concat_alignment.fasta','w')
	newAlignment = Alignment('')
#	newAlignment = ''

	for f in os.listdir(path2alignments):
#		if re.search('.fas',f):
		if f.endswith('.fas'):
			tempAlignment = concat(newAlignment, path2alignments + f,'LG')	
			newAlignment  = tempAlignment
			

	newAlignment.write('concat')
	time.sleep(20)
	makephy('concat_alignment.fasta')
main()