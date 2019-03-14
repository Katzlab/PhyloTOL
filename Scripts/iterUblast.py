import os, sys
from Bio import SeqIO

def sort_cluster(forpairComp, toosim):
	
	# set variables and path
	Path = ('/').join(forpairComp.split('/')[:-1]) + '/'
	OGlenDB_file = open(Path + "/oglengths", "r").readlines()
	
	OGlenDB = {}
	for og_len in OGlenDB_file:
		if og_len.startswith("OG5") : OGlenDB[og_len.split("\t")[0]] = (og_len.split("\t")[1]).replace("\n", "")
	
	OG_taxon = forpairComp.split('/')[-1].replace("_out.txt", "")
	OG = OG_taxon[:10]
	taxon = OG_taxon[11:]
	if toosim == 'n':
		print taxon + ": overlap filter (OF) for " + OG
	else:
		print taxon + ": overlap and similarity filters (OF and SF) for " + OG	
	
	# Create files and folders 
	if not os.path.exists(Path+'/UBlastFiles/'): #  working folder for uBlast
		os.makedirs(Path+'/UBlastFiles/')
	os.system('cp ' + forpairComp + ' ' + Path + '/UBlastFiles/' + OG_taxon + '_sorted.fasta')
	fasta2keep_file = open(Path + '/fasta2keep/' + OG_taxon + 'fastatokeep.fas','a') 
	tsvFiltered_file = open(Path + '/UBlastFiles/'+ OG_taxon + '_resultsFiltered.tsv','a')
		
	run = 0  		
	n = 'exe'		
	
	while n == 'exe': # while the key word is 'exe', it will run iterations of Ublast. 
		run += 1
		master_seq = []		
 
 		# ----------------------- MASTER AND QUERY SEQUENCES  -----------------------------
 		# For the first run take both the master sequence (db for Ublast) 
 		# and the remaining sequences (query for Ublast) from '_out.txt' sorted by size
		if run == 1:		
			FastaDict = {}	
			
			for fastaRecord in sorted(SeqIO.parse(Path + '/UBlastFiles/' + OG_taxon + '_sorted.fasta','fasta'), key=lambda fastaRecord: -len(fastaRecord.seq)):
				if len(fastaRecord.seq) <= 1.5 * float(OGlenDB[OG]):
					FastaDict[fastaRecord.description] = fastaRecord.seq
					if not master_seq:				
						master_seq = [fastaRecord.description, fastaRecord.seq]
			
			starting_seqsnum = len(FastaDict)

		# for each subsequent run of SF, the filtered sequences from previous run will be saved 
		# on FastaDict (the variable is re-set every time). So, the master and query sequences 
		# are picked from FastaDict sorted by size.		
		if run > 1:
		
			if len(FastaDict) >= 1:
				largest = sorted(FastaDict, key=lambda seq_Name: -len(FastaDict[seq_Name]))[0]
				master_seq = [largest, FastaDict[largest]]

			if str(toosim) == 'n':
				for seq in FastaDict:
					fasta2keep_file.write('>%s\n%s\n' % (seq, FastaDict[seq]))
					master_seq = ''
					
		# after every run, the master sequence is saved in fasta2keep. So, there should be at
		# least 1 master sequence or fasta2keep is not produced. 		
		if master_seq:		
			del FastaDict[master_seq[0]]
			fasta2keep_file.write('>%s\n%s\n' % (master_seq[0], master_seq[1]))
			
			# If not more sequences remain after an iteration, then stop the loop
			if FastaDict:
				masterseq = open(Path + '/UBlastFiles/'+ OG_taxon + '_master.fasta','w')
				masterseq.write('>%s\n%s' % (master_seq[0], master_seq[1]))
				masterseq.close()		
				sortedSeqs = open(Path + '/UBlastFiles/' + OG_taxon + '_sorted.fasta', 'w')
				for record in FastaDict:		
					sortedSeqs.write(">%s\n%s\n" % (record, FastaDict[record]))
				sortedSeqs.close()		

		# -------------------------------- RUN Ublast ------------------------------------------
				os.system('usearch -ublast ' + Path + '/UBlastFiles/' + OG_taxon + '_sorted.fasta -db ' + Path + '/UBlastFiles/' + OG_taxon +'_master.fasta -evalue 1 -blast6out ' + Path + '/UBlastFiles/results_' + OG_taxon + str(run) + '.tsv') # running usearch
				tsv = open(Path + '/UBlastFiles/results_' + OG_taxon + str(run) + '.tsv', 'r').readlines()
						
				tsvFiltered = [] # Here we collect Ublast results after every iteration
				filteredRecords = {} # there are going to be the filtered sequences. After every iteration this variable will replace FastaDict
				retained_OF = 1	# We start the counter of retain sequences (after OF) in 1, because first master is already picked
				if tsv:  # If tsv file is empty is because there are not matches after Usearch

					# For first run do OF and SF 
					if run == 1:

						for tsvLine in tsv:
							tsvLine = tsvLine.replace("\n", "")
							alignmentLength = int(tsvLine.split('\t')[3])
							gaps = int(tsvLine.split('\t')[5])
							seqName = str(tsvLine.split('\t')[0]) 
							identity = float(tsvLine.split('\t')[2])
							
							if (alignmentLength - gaps) >= 0.35 * len(master_seq[1]):
								retained_OF += 1
								if str(toosim) == 'n':
									tsvFiltered.append(tsvLine + '\tOF+\t' + str(run) + 'SF_NA')
									filteredRecords[seqName] = FastaDict[seqName]
								else:
									if float(identity) < float(toosim):
										tsvFiltered.append(tsvLine + '\tOF+\t' + 'SF_' + str(run) + '+')
										filteredRecords[seqName] = FastaDict[seqName]
									else:
										tsvFiltered.append(tsvLine + '\tOF+\t' + 'SF_' + str(run) + '-')
							else:								
								tsvFiltered.append(tsvLine + '\tOF-\t' + str(run) + 'SF_NA')
						print "%s: %s/%s sequences retained for %s after OF" % (taxon, retained_OF, starting_seqsnum, OG)
					
					if run > 1:

						for tsvLine in tsv:
							tsvLine = tsvLine.replace("\n", "")
							seqName = str(tsvLine.split('\t')[0]) 
							identity = float(tsvLine.split('\t')[2])
							if float(identity) < float(toosim):
								tsvFiltered.append(tsvLine + '\tOF_NA\t' + 'SF_' + str(run) + '+')
								filteredRecords[seqName] = FastaDict[seqName]
							else:
								tsvFiltered.append(tsvLine + '\tOF_NA\t' + 'SF_' + str(run) + '-')

					os.system('rm ' + Path + '/UBlastFiles/results_' + OG_taxon + str(run) + '.tsv')
					FastaDict = filteredRecords  # Here FastaDict is re-set to contain filtered sequences

		# --------------------------------- Stop iterations --------------------------------------------
				else:  # when no sequences matched the master sequence (in any run). Also delete the empty .tsv
					n = 'stop'
					fasta2keep_file.close()
					tsvFiltered.append("No sequences matched " + master_seq[0])
					os.system('rm ' + Path + '/UBlastFiles/results_' + OG_taxon + str(run) + '.tsv')
				for result in tsvFiltered:
					tsvFiltered_file.write("%s\n" % result)
			
			else:	# when there are not more sequences to iterate
				n = 'stop' 
				fasta2keep_file.close()
				tsvFiltered_file.write("No sequences to compare against " + master_seq[0])
				tsvFiltered_file.close()

		else:	# when there are not sequences that can be picked as master (i.e., don't meet criterion of > 150% of average OG length or not sequences retained after last run)
			n = 'stop'
			fasta2keep_file.close()
			tsvFiltered_file.close()

	retained_seqs = open(Path + '/fasta2keep/' + OG_taxon + 'fastatokeep.fas', 'r').readlines()			
	if not retained_seqs:
		os.system('rm ' + Path + '/fasta2keep/' + OG_taxon + 'fastatokeep.fas')
		print "%s: no sequences met criterion of being smaller than %s, 1.5 of the average length of %s" % (taxon, OGlenDB[OG], OG)		
	else:
		if not str(toosim) == 'n':
			print "%s: %s/%s sequences retained for %s OF and SF" % (taxon, len(retained_seqs)/2, starting_seqsnum, OG)

def main():
	script, forpairComp, toosim = sys.argv
	sort_cluster(forpairComp, toosim)

main()
