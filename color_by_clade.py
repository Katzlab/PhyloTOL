#Importing basic dependencies
import os, re
import sys

#This Biopython library is helpful for rooting & ladderizing trees
from Bio import Phylo

#Checking to see if p4 is installed -- if not, then not a problem, but the only rooting option will be midpoint
force_mid = False
try:
	from p4 import *
	var.doRepairDupedTaxonNames = 1
except:
	print('\nThe p4 library is not installed! Defaulting to midpoint-rooting\n')
	force_mid = True

#Defining some global variables to fill in later
dir_name = ''
file_names = []
minor_clades = []

#Defining the colors! This done using hex codes.
colors = ['[&!color=#0000ff]', '[&!color=#ab2121]', '[&!color=#7b25aa]', '[&!color=#12aaff]', '[&!color=#006300]', '[&!color=#ffa100]', '[&!color=#000000]', '[&!color=#ff6288]', '[&!color=#808080]']

taxon_colors = { 'ba' : '[&!color=#000000]', 'za' : '[&!color=#808080]', 'sr' : '[&!color=#7b2516]', 'op' : '[&!color=#12aaff]', 'pl' : '[&!color=#006300]', 'ex' : '[&!color=#ffa100]', 'ee' : '[&!color=#ff6288]', 'am' : '[&!color=#aa00ff]' }
	

#A function to remove dashes from a tree before parsing it, because p4 cannot read dashes (annoyingly)
def correct_tree(fname): 
	tree2correct = open(fname, 'r')
	tree2correct = tree2correct.readline()
	tree_corrected = open (fname.split('.tre')[0] + '_correcting.tre', 'w')
													
	if '-' in tree2correct:							
		tree2correct = re.sub('-', '', tree2correct) 
		
	tree_corrected.write(tree2correct)
	tree_corrected.close()	
	

#This function converts NEXUS files to Newick strings (in case your trees are already colored or have been saved in FigTree)
def nexus_to_newick(fname):

	nexus = False
	
	#Checking to see if it is a NEXUS file
	for line in open(fname):
		if('NEXUS' in line):
			nexus = True
	
	#If it is a NEXUS file, then rewrite it using the Newick string portion of the file
	if(nexus):
		newick = ''
		for line in open(fname):
			#line = line.split(' ')[-1]
			if(line.startswith('(') or line.startswith('tree1=')):
				newick = line.split('tree1=')[-1].replace("'", '').replace('\\', '')
				
		with open(fname, 'w') as o:
			o.write(newick)
	

#This function does most of the work, and reads in all of the
def read_newick(input_handle, use_minor_clades, bacterial_bin, directory, rooting):

	print('\nProcessing tree OG5_' + input_handle.split('OG5_')[1][:6] + '\n')
	
	#Formatting and reading in the tree to Phylo (Biopython) to ladderize and reroot it
	nexus_to_newick(input_handle)
	temp_handle = input_handle.replace('.tre', '_temp.tre')
	tree = Phylo.read(input_handle, 'newick')
	
	tree.ladderize(reverse = True)
	if(rooting == 1):
		tree.root_at_midpoint()
	Phylo.write(tree, temp_handle, 'nexus')

	#If coloring by minor clade, just filling a dictionary with the colors to use
	mcandcol = {}
	if(use_minor_clades == True):
		try:
			for m, mc in enumerate(minor_clades):	
				mcandcol.update({mc.lower() : colors[m]})
		except IndexError:
			bad_script_call()

	#Where the magic happens -- here, all of the tips in the tree are read and recorded, along with the color that they are supposed to be (note that all of this could be done in one line using p4, but this is very annoying to install)
	newick = ''
	taxa_and_colors = []
	#For all of the lines in the tree file
	for line in open(temp_handle, 'r'):
		temp = line.split(' ')[-1]
		
		#Once we get to the Newick string part of the file
		if(temp.startswith('tree1=')):
			line = line.split('tree1=')[1].replace("'", '').replace('\\', '')
			newick = temp.split('tree1=')[1].replace("'", '').replace('\\', '')
			
			#Splitting the Newick string by commas (each tip/node is separated by a comma)
			line = line.split(',')
			
			#For every node in the Newick string
			for chunk in line:
				#Reformatting the name
				chunk = chunk.split('(')[-1].split(')')[0]
				if("'" in chunk):
					chunk = chunk.split("'")[1]
				
				#Getting the tip name
				try:
					#og_idx = chunk.index('OG5_')
					#taxon = chunk[0:og_idx + 10]
					taxon = chunk.split(':')[0]
				except:
					taxon = chunk[0:10]
					
				#If not using minor clades									
				if(use_minor_clades == False):
					
					#If not coloring bacterial bins
					if(bacterial_bin == False):
						#Get the color for the taxon by the major clade (first two digits)
						color = taxon_colors[taxon[:2].lower()]
						#Recording the tip and the color in a list
						taxa_and_colors.append(taxon + color)
					#If coloring bacterial bins
					else:
						#If the tip is not a eukaryote, then do the same as above
						if(taxon[:2].lower() == 'ba' or taxon[:2].lower() == 'za'):
							color = taxon_colors[taxon[:2].lower()]
							taxa_and_colors.append(taxon + color)
						else:
							#If the tip is bacterial bin, give it the bacterial bin color (first in colors list)
							if(taxon[4] == 'b'):
								taxa_and_colors.append(taxon + colors[0])
							#If not, color the same way as above
							else:
								color = taxon_colors[taxon[:2].lower()]
								taxa_and_colors.append(taxon + color)
				#If coloring by minor clade -- do the same as above, but using the first 5 digits of the taxon instead of the first two
				else:
					try:
						color = mcandcol[taxon[:5].lower()]
						taxa_and_colors.append(taxon + color)
					except KeyError:
						#If the minor clade isn't being colored, then don't color it
						taxa_and_colors.append(taxon)
	
	#All of this is if you want to root by the largest bacterial or second most abundant clade in the tree. You need p4 for this, so you probably won't touch it.
	if(rooting == 0):
		p4_root = True
	
		with open(temp_handle, 'w') as o:
			o.write(newick)
			
		var.trees = []
			
		correct_tree(temp_handle)				
												
		tree_file = temp_handle.split('.tre')[0] + '_correcting.tre'
		read(tree_file) 
		tree = var.trees[0]
		
		os.remove(tree_file)
		
		sizes_cladesBaZa = {}			
		sizes_cladesOp = {}				
		sizes_cladesPl = {}				
		sizes_cladesAm = {}				
		sizes_cladesEx = {}				
		sizes_cladesSr = {}				
				
		for node in tree.iterNodesNoRoot():				
			allTaxa_node = tree.getAllLeafNames(node)	
			
			MC_list =[]

			for taxon in allTaxa_node:					
				MC_list.append(taxon[:2])	
																
			if MC_list.count('Ba') + MC_list.count('Za') == (len(MC_list) or len(MC_list)-1):
				sizes_cladesBaZa[tree.node(node)] = len(allTaxa_node)	
			if MC_list.count('Op') == (len(MC_list) or len(MC_list)-1):
				sizes_cladesOp[tree.node(node)] = len(allTaxa_node)
			if MC_list.count('Pl') == (len(MC_list) or len(MC_list)-1):
				sizes_cladesPl[tree.node(node)] = len(allTaxa_node)
			if MC_list.count('Am') == (len(MC_list) or len(MC_list)-1):
				sizes_cladesAm[tree.node(node)] = len(allTaxa_node)
			if MC_list.count('Ex') == (len(MC_list) or len(MC_list)-1):
				sizes_cladesEx[tree.node(node)] = len(allTaxa_node)
			if MC_list.count('Sr') == (len(MC_list) or len(MC_list)-1):
				sizes_cladesSr[tree.node(node)] = len(allTaxa_node)
						
		if sizes_cladesBaZa: 		
			biggestBaZa = max(sizes_cladesBaZa, key=sizes_cladesBaZa.get) 
			tree.reRoot(biggestBaZa)
		else:
			if sizes_cladesOp:
				biggestOp = max(sizes_cladesOp, key=sizes_cladesOp.get)
				tree.reRoot(biggestOp)
			else:
				if sizes_cladesPl:
					biggestPl = max(sizes_cladesPl, key=sizes_cladesPl.get)
					tree.reRoot(biggestPl)
				else:
					if sizes_cladesAm:
						biggestAm = max(sizes_cladesAm, key=sizes_cladesAm.get)
						tree.reRoot(biggestAm)
					else:
						if sizes_cladesEx:
							biggestEx = max(sizes_cladesEx, key=sizes_cladesEx.get)
							tree.reRoot(biggestEx)
						else:
							if sizes_cladesSr:
								biggestSr = max(sizes_cladesSr, key=sizes_cladesSr.get)
								tree.reRoot(biggestSr)
							else:
								print('Very few or no taxa found in tree')
								p4_root = False
		
		if(p4_root):				
			tree.writeNewick(temp_handle)					
			newick = open(temp_handle).read()
				
	os.remove(temp_handle)
					
	#Returning the newick string as well as the list of taxa and their corresponding colors. The hard part is done!
	return [newick, taxa_and_colors]

#Writing the output file with the NEXUS format
def write_lines(o, newick, taxa_and_colors, tree_font_size):
	ntax = str(len(taxa_and_colors))
			
	o.write('#NEXUS\n')	
	o.write('begin taxa;\n')
	o.write('\tdimensions ntax=' + ntax + ';\n')
	o.write('\ttaxlabels\n')
			
	for taxon in taxa_and_colors:
		o.write('\t' + taxon + '\n')
			
	o.write(';\nend;\n\n')
			
	o.write('begin trees;\n')
	o.write('\ttree tree_1 = [&R]\n')
	o.write(newick.replace('LKH', '-LKH').replace('--LKH', '-LKH'))
	o.write('end;\n\n')
			
	with open('figtree_format.txt', 'r') as ff:
		for line in ff:
			if('.fontSize' in line):
				o.write(line.replace('8', tree_font_size))
			else:
				o.write(line)

#A wrapper for writing a nexus file
def write_nexus(newick, taxa_and_colors, tree_font_size, dir_name = '', directory = False, name_idx = 0):

	if(directory == False):
		out_path = file_names[0][:-4] + '_colored.tre'
	else:
		if(not os.path.isdir(str(dir_name) + '_colored')):
			os.mkdir(str(dir_name) + '_colored')
			
		out_path = str(dir_name) + '_colored/' + str(file_names[name_idx][:-4]) + '_colored.tre'
		
	with open(out_path, 'w') as o:
		write_lines(o, newick, taxa_and_colors, tree_font_size)
					

#The main wrapper scripts to call all of the above functions
def main():
	
	#Reading in all of the arguments
	dth, tree_font_size = [sys.argv[1], sys.argv[2]]

	bacterial_bin = False
	directory = True
	rooting = 1
	use_minor_clades = False
	
	f = 0
	#For every file in the folder
	for file in os.listdir(dth):
		#Make sure that it is a tree
		if(file.endswith('.tre')):
			fname = dth + '/' + file
			#Recording the file name (I forget why)
			file_names.append(file)
			
			#Calling an above function to get the Newick string for the file and a list of all of the taxa and their corresponding colors
			newick, taxa_and_colors = read_newick(fname, use_minor_clades, bacterial_bin, directory, rooting)
			
			#Writing out the final NEXUS file
			write_nexus(newick, taxa_and_colors, tree_font_size, dth, directory, f)
			f += 1
	
#Calling the main function. Nothing will happen if you don't do this.
main()

