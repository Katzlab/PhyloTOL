#!/usr/bin/python
from __future__ import print_function

__author__ = "Jean-David Grattepanche"
__version__ = "2, August 28, 2017"
__email__ = "jeandavid.grattepanche@gmail.com"



import sys
import os
import re
import time
import string
import os.path
from Bio import SeqIO
from sys import argv


def Addcoverage(code):
	for Folder in os.listdir(os.curdir+'./'):
# 		print(Folder)
		if code in Folder:
			seqfolder = os.curdir+'./'+Folder+'/Processed'	
			covupd = {}
			try:
				for seqcoll in open(os.curdir+'./'+Folder+'/'+Folder+'_SeqPairsAbove98.txt','r'):
					CL = 0
					for transc in seqcoll.split('\t'):
						if CL == 0:
							reftrans = ('_').join(transc.split('_')[1:])
						coverage = int(transc.split('Cov')[1].split('_')[0])
						Length = int(transc.split('Len')[1].split('_')[0])
						CL += coverage * Length
					covupd[reftrans] = CL 
			#	print(('/').join(seqfolder.split('/')[:-1]) +'/Updated_Coverage/')
			except:
				print("NO PARTIAL FOR", Folder)					

			if os.path.isdir(('/').join(seqfolder.split('/')[:-1]) +'/Updated_Coverage/') != True:
				os.system('mkdir '+('/').join(seqfolder.split('/')[:-1]) +'/Updated_Coverage/')
			if os.path.isdir(('/').join(seqfolder.split('/')[:-1]) +'/Updated_Coverage/SpreadSheets/') != True:
				os.system('mkdir '+('/').join(seqfolder.split('/')[:-1]) +'/Updated_Coverage/SpreadSheets/')
			for spreadsh in os.listdir(seqfolder+'/SpreadSheets/'):
				if spreadsh.endswith('.tsv'):
					outtsvtokeep = open(('/').join(seqfolder.split('/')[:-1])+'/Updated_Coverage/SpreadSheets/'+spreadsh.split('Final')[0] +'UC.Final'+spreadsh.split('Final')[1],'w+')
					for row in open(seqfolder+'/SpreadSheets/'+ spreadsh, 'r'):
						if row.split('_Trans')[0] in covupd:
							newcov2= covupd[row.split('_Trans')[0]] / int(row.split('_Len')[1].split('_')[0])
	# 							print(row.split('_Trans')[0], newcov2)
							outtsvtokeep.write(row.split('Cov')[0]+'Cov'+str(newcov2)+'_OG5' +row.split('OG5')[1].split('_Trans')[0] +'\t' +('\t').join(row.split('\t')[1:]))
						else: 
							if 'Trans' in row:
								outtsvtokeep.write(row.split('_Trans')[0]+ '\t' +('\t').join(row.split('\t')[1:]))
		 					else:
								outtsvtokeep.write(row)
#		 						print(row)
					outtsvtokeep.close()
			for seqfile in os.listdir(seqfolder):
				if seqfile.endswith('.fasta'):
	# 					print(seqfile)
					outseqtokeep = open(('/').join(seqfolder.split('/')[:-1])+'/Updated_Coverage/'+seqfile.split('Final')[0] +'UC.Final'+seqfile.split('Final')[1],'w+')
					for Seq in SeqIO.parse(seqfolder+'/'+seqfile ,'fasta'):
						if Seq.description.split('_Trans')[0] not in covupd:
							outseqtokeep.write('>'+Seq.description.split('_Trans')[0]+ '\n'+str(Seq.seq) +'\n')
						else:
							newcov= covupd[Seq.description.split('_Trans')[0]] / int(Seq.description.split('_Len')[1].split('_')[0])
	# 							print(Seq.description,'\t',newcov)
							outseqtokeep.write('>'+Seq.description.split('Cov')[0]+'Cov'+str(newcov)+'_OG5' +Seq.description.split('OG5')[1].split('_Trans')[0]+ '\n'+str(Seq.seq) +'\n')
					outseqtokeep.close()	
			if os.path.isdir('../ToRename') != True:
				os.system('mkdir ../ToRename')
			os.system('cp ../'+Folder+'/Updated_Coverage/*fasta ../ToRename/')
			os.system('cp ../'+Folder+'/Updated_Coverage/SpreadSheets/*tsv ../ToRename/')
			
			
def main():
	print("\n*******************************************************************")
	print("This script:  \n\t  python update_cov_post_removepartials.py 10digitcode\n")
	print("\n\nExample:\n python update_cov_post_removepartials.py Am_tu_Avu \t ")
	print("*******************************************************************\n")
	script, code = argv
	Addcoverage(code)
main()