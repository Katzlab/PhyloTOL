import os
import sys
from Bio import Phylo, SeqIO
import operator

def get_clade(tree, sequence):
	for clade in tree.get_terminals():		
		if str(clade.name).strip() == str(sequence.description).strip():
			return clade
		
def get_parent(tree, child_clade):
	node_path = tree.get_path(child_clade)
	try:
		return node_path[-2]				
	except:
		return None


def is_monophyletic(tree,clade,recursion_depth,child,length,sequence,mcCladeList,seqCladeList):
	#############################################################
	#checks to see if the clade is monophyletic for major clade and, if it is, 
	#recursively checks the parent clade until it finds a clade that is not monophyletic
	#############################################################
			
	out = open('clade_size.txt','a')
	try:
		for seq in clade.get_terminals(): #for each terminal sequence in the clade
			mc = str(seq)[:4] #ui is the major clade identifier i.e Op, Pl, etc
			
			if str(seq) not in seqCladeList:
				seqCladeList.append(str(seq)) #add the sequence to the seq list
				mcCladeList.append(mc) #add the unique identifier to the ui list
	
		parent = clade
						
		if len(mcCladeList) >= 1:
			if all(x == mcCladeList[0] for x in mcCladeList):
				length = len(mcCladeList)
				child = clade			
				parent = get_parent(tree,clade)	
				recursion_depth = recursion_depth + 1			
				is_monophyletic(tree,parent,recursion_depth,child,length,sequence,mcCladeList,seqCladeList)
			else:
				if recursion_depth < 1:
					out.write(sequence.strip() + ':' + str(1) + '\n')
				else:	
					out.write(sequence.strip() + ':' + str(length) + '\n')
	except:
		out.write(sequence.strip() + ':' + str(0) + '\n')
		
	out.close()
		
		
def get_seq_to_keep(tree):
	#############################################################	
	#for each taxon, makes a list of seq/num pairs sorted py num
	#and passes it to get_seq_to_keep which will take the best seq
	#############################################################
	
	inFile = open('clade_size.txt','r')
	infile = inFile.readlines()
	inFile.close()
	seqDict = {}
	seqnumDict = {}
	
	for line in infile:
		seq = line.split(':')[0]
		num = line.split(':')[1]
		seqDict[seq] = int(num)

	sorted_seqDict = sorted(seqDict.items(), key=operator.itemgetter(1), reverse = True)
	xx = open('clade_size.txt','w') #clear file
	xx.close()	
	#############################################################
	#checks to see if there is one sequence in a larger monophyletic clade
	#if not, makes a list of seqs and takes the one with the smallest 
	#branch length
	#############################################################
	
	check_dist_list=[]
	seqs2remove = []
	if len(sorted_seqDict) > 1:			# sorted_seqDict, list of (seqs,sizeofclade) sorted by size of clade
		if float(sorted_seqDict[0][1]) > float(sorted_seqDict[1][1]): #if the first is bigger than the next then it is the biggest
			sorted_seqDict.pop(0) #keep the first, get rid of the rest
			for double in sorted_seqDict:  
				seqs2remove.append(double[0])  # put all the rest in seqs2remove
		else: #More than seq in clade with equal numbers
			for double in sorted_seqDict:
				if float(double[1]) == float(sorted_seqDict[0][1]):
					check_dist_list.append(double) #get list of seqs to check for distance
				else:
					seqs2remove.append(double[0])
			distance_dict = {}
			for double in check_dist_list:
				for element in (tree.find_elements(name = double[0])):
					distance_dict[element] = tree.distance(element)
		
			sorted_distance_dict = sorted(distance_dict.items(), key=operator.itemgetter(1))
			
			sorted_distance_dict.pop(0) #keep the shortest, remove others
			for element in sorted_distance_dict:
				seqs2remove.append(str(element[0].name)) #want just the number, or SeqCode
			
	return seqs2remove
	
	
def remove_paralogs(results2keep, forconcat):
	
	for file in os.listdir(results2keep):
		if('gapTrim' in file and not file.endswith('.tre') and 'OG5_' in file):
			og = 'OG5_' + file.split('OG5_')[-1][:6]
			
			try:
				tre_f = [f for f in os.listdir(results2keep) if f.endswith('.tre') and og in f][0]
			except IndexError:
				print('Tree not found for OG ' + og)
				continue
				
			tree = Phylo.read(results2keep + tre_f, 'newick')
			tree.root_at_midpoint()
			
			ui_list = []
			paralogDict = { }
			for rec in SeqIO.parse(results2keep + file, 'fasta'):
				try:
					ui = rec.description[:10] #ui is MC_mc_code
					if(ui not in paralogDict):
						paralogDict.update({ ui : [] })
					paralogDict[ui].append(rec) # so len is # of paralogs per taxon
					if ui not in ui_list:
						ui_list.append(ui)
				except:
					pass
					
			seqs2remove = []
			print('Processing ' + og + '... ')
			for ui in ui_list:	
				if len(paralogDict[ui]) > 1:
					for sequence in paralogDict[ui]:
						clade = get_clade(tree, sequence)
						try:
							seq_clade = get_parent(tree, clade)
							if seq_clade == None:
								seq_clade = clade
						except:
							seq_clade = clade
																										
						if seq_clade != None:
							is_monophyletic(tree,seq_clade,0,"None",0,str(sequence.description),[],[])
						else:
							seqs2remove.append(sequence.description[:10])
						paralogDict[ui] = []
					seqs2remove += get_seq_to_keep(tree)
								
			with open(forconcat + og + '_paralog_purged.fasta', 'w') as o:
				for rec in SeqIO.parse(results2keep + file, 'fasta'):
					if(rec.description not in seqs2remove):
						o.write('>' + rec.description + '\n' + str(rec.seq) + '\n\n')
						
		