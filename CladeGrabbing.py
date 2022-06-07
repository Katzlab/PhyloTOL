#Dependencies
from p4 import *
import os, re
import sys
import csv
from Bio import SeqIO
import subprocess

var.doRepairDupedTaxonNames = 1
	
	
def nexus_to_newick(fname):

	nexus = False
	
	for line in open(fname):
		if('NEXUS' in line):
			nexus = True
			
	if(nexus):
		newick = ''
		for line in open(fname):
			line = line.split(' ')[-1]
			if(line.startswith('(') or line.startswith('tree1=')):
				newick = line.split('tree1=')[-1].replace("'", '').replace('\\', '')
				
		with open(fname, 'w') as o:
			o.write(newick)

	
def correct_tree(fname): 
	tree2correct = open(fname, 'r')
	tree2correct = tree2correct.readline()
	tree_corrected = open (fname.split('.tre')[0] + '_temp.tre', 'w')
													
	if '-' in tree2correct:							
		tree2correct = re.sub('-', '', tree2correct) 
		
	if(tree2correct[0].startswith('(') or tree2correct.startswith('\t(')):
		newick = tree2correct
		newick = newick.split(',')
		
		to_join = []
		for leaf in newick:
			if(leaf.startswith('(') and ')' in leaf):
				leaf = leaf.replace('(', '').split(')')
				del leaf[1]
				leaf = ')'.join(leaf)
				to_join.append(leaf)
			else:
				to_join.append(leaf)
				
		tree2correct = ','.join(to_join)
	
	tree_corrected.write(tree2correct)
	tree_corrected.close()	
	
	
def get_best_clade(PathtoFiles, tree, cont_num_contams, target_clade, nodes_to_exclude, min_presence, target_taxa_list, at_least_file, num_at_least):

	if(target_taxa_list.lower() != 'na' and target_taxa_list != ''):
		if(os.path.isfile(PathtoFiles + '/' + target_taxa_list)):
			target_taxa_list = [line.strip() for line in open(PathtoFiles + '/' + target_taxa_list)]
		else:
			print('\nClade-grabbing contamination loop error: it looks like you tried to input a target_taxa_file file, but it could not be found in the "DataFiles" folder. Fix the path here or make sure that this parameter is set to "NA".\n')
			quit()

	target_minor = target_clade[:4]

	if(at_least_file.lower() != 'na' and at_least_file != ''):
		if(os.path.isfile(PathtoFiles + '/' + at_least_file)):
			at_least_list = [line.strip() for line in open(PathtoFiles + '/' + at_least_file)]
		else:
			print('\nClade-grabbing contamination loop error: it looks like you tried to input an at_least_sisters_file file, but it could not be found in the "DataFiles" folder. Fix the path here or make sure that this parameter is set to "NA".\n')
			quit()
	else:
		at_least_list = []
		num_at_least = 0
		
	# for tax in at_least_list[::-1]:
	# 	if(tax[:10] in target_taxa_list or tax[:8] in target_taxa_list):
	# 		at_least_list.remove(tax)

	forbidden_nodes = [node for node in nodes_to_exclude]
	for node in nodes_to_exclude:
		for num in tree.getNodeNumsAbove(node):
			forbidden_nodes.append(tree.node(num))
	
	best_node = None
	best_size = 0
	for node in tree.iterNodesNoRoot():
		if(node not in forbidden_nodes):
			leaves = tree.getAllLeafNames(node)
							
			num = 0.0; dem = 0.0;
		
			non_minor = []
			for leaf in leaves:
				if(leaf[:2] != target_minor and leaf[:4] != target_minor):
					num += 1.0;
					non_minor.append(leaf[:10])
					
			if(target_taxa_list == 'na' or target_taxa_list == '' or target_taxa_list == 'NA'):
				n_targets = len(list(dict.fromkeys([tip[:10] for tip in leaves if(tip[:2] == target_clade or tip[:3] == target_clade or tip[:4] == target_clade or tip[:5] == target_clade or tip[:7] == target_clade or tip[:8] == target_clade)])))
			else:
				n_targets = len(list(dict.fromkeys([tip[:10] for tip in leaves if((tip[:2] == target_clade or tip[:3] == target_clade or tip[:4] == target_clade or tip[:5] == target_clade or tip[:7] == target_clade or tip[:8] == target_clade) and (tip[:10] in target_taxa_list or tip[:8] in target_taxa_list))])))
							
			at_least_taxa = len(list(dict.fromkeys([leaf[:10] for leaf in tree.getAllLeafNames(node) if leaf[:10] in at_least_list])))

			if(num <= cont_num_contams and n_targets > best_size and n_targets >= min_presence and at_least_taxa >= num_at_least):
				best_node = node
				best_size = n_targets
			
	return best_node
	
def get_subtree(PathtoFiles, fname, cont_num_contams, target_clade, min_presence, target_taxa_list, at_least_file, num_at_least, cont_filter_differential_coverage, cont_coverage_diff_OM, cont_coverage_diff_abs, return_cladegrabbing_subtrees):
		
	var.trees = []
	
	nexus_to_newick(fname)
														
	correct_tree(fname)							
												
	tree_file = fname.split('.tre')[0] + '_temp.tre'
	read(tree_file) 
	tree = var.trees[0] 						
												
	os.remove(tree_file)
		
	seen_clades = []
	
	best_node = get_best_clade(PathtoFiles, tree, cont_num_contams, target_clade, seen_clades, min_presence, target_taxa_list, at_least_file, num_at_least)
	if(best_node != None):
		seen_clades.append(best_node)
						
	while best_node != None:
		best_node = get_best_clade(PathtoFiles, tree, cont_num_contams, target_clade, seen_clades, min_presence, target_taxa_list, at_least_file, num_at_least)
		seen_clades.append(best_node)
		
	seqs2keep = []
	if(len(seen_clades) == 0):
		print('\nFound no clades that meet the contamination-loop selection criteria in the tree ' + fname.split('/')[-1] + '\n')
	else:
		print('\nFound ' + str(len(seen_clades) - 1) + ' clades that meet the contamination-loop selection criteria in the tree ' + fname.split('/')[-1] + '\n')
	for n, node in enumerate(seen_clades):
		if(node != None):
			if(cont_filter_differential_coverage):
				print('Filtering for differential coverage within clade...\n')
				all_taxa = [tip[:10] for tip in tree.getAllLeafNames(node)]
				tax_counts = { tax : all_taxa.count(tax) for tax in all_taxa}
				duped_taxa = [tax for tax in tax_counts if tax_counts[tax] > 1]
				tip_groups = [[tip for tip in tree.getAllLeafNames(node) if tax in tip] for tax in duped_taxa]

				tips_to_exclude = []
				for group in tip_groups:
					cov_dict = { }
					for tip in group:
						if('Cov' in tip[10:]):
							cov_dict.update({ tip : int(tip.split('Cov')[-1].split('_')[0].split('.')[0]) })

							for t2 in cov_dict:
								if(cont_coverage_diff_abs not in ['na', 'NA'] and cont_coverage_diff_OM not in ['na', 'NA']):
									if(max(cov_dict.values()) / cov_dict[t2] >= 10**int(cont_coverage_diff_OM) or max(cov_dict.values()) - cov_dict[t2] >= cont_coverage_diff_abs):
										if(t2 not in tips_to_exclude):
											tips_to_exclude.append(t2)
											print(t2 + ' is being excluded due to highly differential coverage with other target taxa in the clade:\n\tCoverage of this sequence: ' + str(cov_dict[t2]) + '\n\tBest coverage: ' + str(max(cov_dict.values())))
								elif(cont_coverage_diff_abs not in ['na', 'NA']):
									if(max(cov_dict.values()) - cov_dict[t2] >= cont_coverage_diff_abs):
										if(t2 not in tips_to_exclude):
											tips_to_exclude.append(t2)
											print(t2 + ' is being excluded due to highly differential coverage with other target taxa in the clade:\n\tCoverage of this sequence: ' + str(cov_dict[t2]) + '\n\tBest coverage: ' + str(max(cov_dict.values())))
								elif(cont_coverage_diff_OM not in ['na', 'NA']):
									if(max(cov_dict.values()) / cov_dict[t2] >= 10**int(cont_coverage_diff_OM)):
										if(t2 not in tips_to_exclude):
											tips_to_exclude.append(t2)
											print(t2 + ' is being excluded due to highly differential coverage with other target taxa in the clade:\n\tCoverage of this sequence: ' + str(cov_dict[t2]) + '\n\tBest coverage: ' + str(max(cov_dict.values())))
						else:
							break

				print('\n' + str(len(tips_to_exclude)) + ' tips excluded for highly differential coverage within the clade\n')

				seqs2keep += [tip for tip in tree.getAllLeafNames(node) if tip not in tips_to_exclude]
			else:	
				seqs2keep += [tip for tip in tree.getAllLeafNames(node)]

			if(return_cladegrabbing_subtrees == 'y'):
				clade_tree = tree.dupeSubTree(node, up = True)
				clade_tree.writeNewick('CladeGrabbingSubtrees/' + fname.split('.tre')[0].split('/')[-1] + '_best_clade_' + str(n) + '.tre')

	seqs2remove = [tip.replace('LKH', '-LKH').replace('--LKH', '-LKH') for tip in tree.getAllLeafNames(0) if tip not in seqs2keep and (tip[:2] == target_clade or tip[:3] == target_clade or tip[:4] == target_clade or tip[:5] == target_clade or tip[:7] == target_clade or tip[:8] == target_clade)]		
	
	return seqs2remove


def rewrite_preguidance(trees_handle, all_seqs2remove, include_outgroup, target_clade):

	all_seqs2remove = [seq.replace('-', '') for seq in all_seqs2remove]

	print('\nRewriting pre-Guidance files...\n')

	recs2keep_by_file = { }; files2delete = []
	for file in os.listdir(trees_handle):
		if('OG5_' in file and 'preguidance' in file):
			recs2keep_by_file.update({ file : [] })
			for rec in SeqIO.parse(trees_handle + '/' + file, 'fasta'):
				if(rec.description.replace('-', '') not in all_seqs2remove):
					if(include_outgroup.lower() == 'y' or include_outgroup.lower() == 'true'):
						recs2keep_by_file[file].append(rec)
					elif(include_outgroup.lower() != 'y' and include_outgroup.lower() != 'true' and target_clade in rec.description[:10]):
						recs2keep_by_file[file].append(rec)
			files2delete.append(file)

	for file in files2delete:
		os.remove(trees_handle + '/' + file)

	for file in recs2keep_by_file:
		with open(trees_handle + '/' + file, 'w') as o:
			for rec in recs2keep_by_file[file]:
				o.write('>' + rec.description + '\n' + str(rec.seq) + '\n\n')

		
def CladeGrabbingWrapper(PathtoFiles, trees_handle, cont_num_contams, target_clade, min_presence, at_least_file, num_at_least, target_taxa_list, include_outgroup, cont_filter_differential_coverage, cont_coverage_diff_OM, cont_coverage_diff_abs, return_cladegrabbing_subtrees):

	#script_call, PathtoFiles, trees_handle, cont_num_contams, target_clade, min_presence, at_least_file, num_at_least, target_taxa_list, include_outgroup, cont_filter_differential_coverage, cont_coverage_diff_OM, cont_coverage_diff_abs, return_cladegrabbing_subtrees = sys.argv
	
	try:
		min_presence = int(min_presence)
	except ValueError:
		print('\nInvalid input for contamination-loop parameter "min_target_presence". This required argument must be an integer value.\n')
		return 'FAILED'

	if(at_least_file.lower() != 'na'):
		try:
			num_at_least = int(num_at_least)
		except ValueError:
			print('\nInvalid input for contamination-loop parameter "at_least_sisters_num". This must be an integer value or NA.\n')
			return 'FAILED'
	else:
		print('\nWARNING: grabbing clades without an at-least list. The use of this argument is recommended\n')

	try:
		cont_num_contams = int(cont_num_contams)
	except ValueError:
		print('\nInvalid input for contamination-loop parameter "num_contams". This required arument must be an integer value.\n')
		return 'FAILED'

	if(cont_filter_differential_coverage == 'y'):
		if(cont_coverage_diff_OM.lower() == 'na' and cont_coverage_diff_abs.lower() == 'na'):
			print('\nIf you want to filter by differential coverage, you must either give an order of magnitude or absolute value by which your coverage values must differ\n')
		elif(cont_coverage_diff_OM.lower() != 'na'):
			try:
				cont_coverage_diff_OM = int(cont_coverage_diff_OM)
			except ValueError:
				print('\nInvalid input for contamination-loop parameter "cont_coverage_diff_OM". This must be an integer value.\n')
				return 'FAILED'
		elif(cont_coverage_diff_abs.lower() != 'na'):
			try:
				cont_coverage_diff_abs = int(cont_coverage_diff_abs)
			except ValueError:
				print('\nInvalid input for contamination-loop parameter "cont_coverage_diff_abs". This must be an integer value.\n')
				return 'FAILED'
	else:
		print('\nGrabbing clades without filtering by differential coverage\n')

	if(return_cladegrabbing_subtrees == 'y' and not os.path.isdir('CladeGrabbingSubtrees')):
		os.mkdir('CladeGrabbingSubtrees')

	seqs2remove_out = open('seqs2remove_out', 'w')
	seqs2remove_out_treesWcont = open('seqs2remove_out_treesWcont', 'w')
	
	f = 0
	all_seqs2remove = []
	for file in os.listdir(trees_handle):
		if('OG5_' in file and '.tre' in file):
			print(str(f + 1) + '. ' + file)
			f += 1
				
			seqs2remove = get_subtree(PathtoFiles, trees_handle + '/' + file, cont_num_contams, target_clade, min_presence, target_taxa_list, at_least_file, num_at_least, cont_filter_differential_coverage, cont_coverage_diff_OM, cont_coverage_diff_abs, return_cladegrabbing_subtrees)
			if(len(seqs2remove) > 0):
				seqs2remove_out_treesWcont.write('OG5_' + file.split('OG5_')[-1][:6] + '\n')
				for seq in seqs2remove:
					seqs2remove_out.write(seq + '\n')
					all_seqs2remove.append(seq)
					
	seqs2remove_out.close()	
	seqs2remove_out_treesWcont.close()	

	rewrite_preguidance(trees_handle, all_seqs2remove, include_outgroup, target_clade)										
		
