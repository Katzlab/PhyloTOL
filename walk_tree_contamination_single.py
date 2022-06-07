# Takes phylogenetic trees and list the sister clade of every taxon. 
# - When the sister clade if monophyletic, it reports the sister's taxa minor clade or 'same clade' if the sister taxa are the same than the compared taxon. 
# - If the sister clade is polyphyletic, it reports all minor clades of the sister taxa.
#
# Some chaneges were made to work with single cells:
# - In LKH data when sister has same first 8 digits in the name, it takes the whole clade and reports its sister. 
# - For all other taxa when the 10 digits codes is the same for the sister, it takes the whole clade and reports its sister.
# - Now it reports branch lengths

# python walk_tree_contamination.py treesFolder output
# python walk_tree_contamination.py ./trees ./output

#import dendropy
from p4 import *
import os, re
import csv
var.doRepairDupedTaxonNames = 1

if(len(sys.argv) == 8):
	script, treesFolder, out, query_clades, sister_clades, branch_length_filter, break_up_clades, single_sister_only = sys.argv
	summarize_sisters = True
else:
	summarize_sisters = False
	script, treesFolder, out = sys.argv


report = open(out, 'w')	# Here is the output

# In the next line you should list the weird taxa
#weirdtaxalist = ["Ex_pa_Tfoe", "Sr_st_tpse", "Sr_st_Csub", "Am_di_Naes", "Op_me_hsap", "Op_me_sman", "Op_me_cele", "Op_me_Dpul", "Op_me_Ctel", "Op_me_tadh", "Op_fu_Amac", "Ex_eu_Bsal", "Ex_he_Ngru", "EE_ap_Ttra", "Sr_st_Bpac", "Pl_gr_Pspg", "Op_ch_mbre", "EE_is_Drot", "Sr_st_Dspe", "EE_cr_Gcry", "Sr_st_Pinf", "Sr_di_Hsps", "Sr_rh_Sspa", "Sr_di_Aspi", "Sr_di_Gcat", "Sr_ci_Ptet", "Sr_ci_Scer", "Ex_pa_tvag", "Sr_ap_Cpar", "Sr_ap_pfal", "Ex_ma_Mjak", "Ex_is_Tpyr", "Pl_rh_Ccho", "Pl_rh_Rmar", "Pl_rh_Ccoe", "Pl_rh_Gsul", "Sr_rh_Lvor", "EE_ce_Rhet", "Am_di_Mspa", "Am_is_Fnol", "EE_ce_Chsp", "EE_ap_Mpla", "EE_ap_Rram", "Ex_ox_Mono", "Sr_rh_Bmot", "EE_ap_Nlon", "EE_is_Tmar", "Ex_eu_Egym", "Sr_ci_Slem", "Op_fu_Aalg", "Am_my_Dpur", "Sr_ch_Vbra", "EE_ap_Ftro", "Am_is_Fflu", "Sr_st_Esil", "Za_as_Heia", "Za_as_Thob", "Za_as_Loki", "Za_as_Odin", "Sr_pe_Perk", "Sr_rh_Cten", "Sr_st_Cfra", "Sr_st_Aana", "Am_ar_Enut", "Am_ar_Mbal", "Op_fu_Bden", "Pl_gr_Atri", "Sr_st_Ospa", "Op_ch_Sros", "Pl_gr_Cvar", "Op_me_Hvul", "Sr_rh_Bnat", "Sr_st_Goce", "Ex_ja_Rame", "Am_my_Ppol", "Pl_gr_Pcol", "Ex_eu_linf", "Ex_eu_tcon", "Pl_gl_Cpad", "Sr_rh_Astr", "Sr_rh_Erot", "Pl_gr_atha", "Pl_gr_Mpol", "Pl_gr_Tchu", "Pl_gr_crei", "Am_di_Acas", "EE_is_Tglo", "Am_tu_Nabe", "Op_ic_Sarc", "Op_is_Mvib", "Op_me_Cfol", "Op_ic_Cowc", "Op_me_Cpul", "Op_fu_Dspa", "Op_fu_Npat", "Op_fu_Ccor", "Op_fu_scer", "Op_fu_Lcor", "Pl_gl_Gnos", "Op_me_Skow", "Op_ic_Apar", "Sr_rh_Asco", "EE_ha_Ehux", "EE_ha_Igal", "EE_is_Tsub", "EE_br_Bant", "EE_ka_Rtru", "Sr_st_Ngad", "Sr_st_Bhom", "Sr_st_Espi", "Sr_st_Ptri", "Sr_st_Spus", "Sr_di_Omar", "Sr_st_Aman", "Sr_st_Croe", "Op_fu_Rall", "Op_me_Ppil", "Ex_fo_Sbar", "Sr_di_Smic", "Sr_ap_Gnip", "Sr_st_Ppar", "EE_br_Stet", "Ba_pg_Abau", "Ba_ac_Cdip", "Ba_pd_Daes", "Ba_cy_Pspb", "Ba_cy_Acyl", "Ba_cy_Onig", "Ba_ch_Caur", "Ba_di_Dtur", "Ba_fb_Gkau", "Ba_ni_Tyel", "Za_ko_ckor", "Ba_sp_Sple", "Ba_pl_Plim", "Ba_bc_Ctha", "Ba_fu_Fnuc", "Ba_te_Alai", "Ba_pg_ecol", "Ba_pb_bpse", "Ba_pa_Abra", "Ba_de_Trad", "Ba_th_tmar", "Ba_pb_Vpar", "Ba_ba_Bfra", "Ba_fc_Oval", "EE_ap_Asig", "Ba_cd_Cmur", "Za_eh_Haci", "Za_eh_Ngre", "Ba_aq_aaeo", "Za_et_Tkod", "Za_ec_Minf", "Za_ey_Mkan", "Za_eb_Mspa", "Za_em_Mhol", "Za_cr_Sisl", "Za_th_Csym", "Ex_eu_Dpap", "Za_ep_tvol", "Za_na_nequ", "Ba_cv_Amuc", "Ba_pa_rpro", "Za_pa_Maci", "Za_cr_Tneu", "Sr_ap_Bequ", "Za_ba_Crea", "Op_nu_Falb"]
weirdtaxalist = [] # Leave empty if you want to run script in all taxa included in your trees

print("\nwalk_tree_contamination_single.py: Detects sister taxa per clade in every tree ...")

if weirdtaxalist == []: 
	taxlist = 'no'
else:
	taxlist = 'yes'

def get_clades(taxon):								# Get different forms of the taxa/clade identifier
	major = taxon.split('_')[0]   					# The identifier as major clade (e.g., Am, Ex)
	minor = taxon.split('_')[1]	  					# The identifier as minor clade (e.g., di, eu)
	sp = taxon.split('_')[2]						# The identifier as 'species' (eg., Acas, Bsal)
	clade_name = major + '_' + minor 				# The identifier as Major_minor (eg., Am_di, Ex_eu)
	taxon_name = major + '_' + minor + '_' + sp 	# The identifier as Major_minor_sp (eg., Am_di_Acas) 
	return (major, minor, sp, clade_name, taxon_name)

def correct_tree(t): 								# Some trees have '-', which causes problems in p4.
	tree2correct = open ('%s%s' % (treesFolder, t), 'r')	# Open the file 't' and takes the information
	tree2correct = tree2correct.readline()			# Reads the info as tree for python 
	tree_corrected = open ('%stemporal_%s' % (treesFolder, t), 'w')
													# The line above creates a new file that will contain ...
													# ... the corrected tree. This is a TEMPORAL file
	if '-' in tree2correct:							# If the tree has '-', delete it
		tree2correct = re.sub('-', '', tree2correct) 
	
	tree_corrected.write(tree2correct)				# Write the corrected tree in the temporal file 
	tree_corrected.close()							# Finish the writing. 

n = 0	# This counter should be initizalize in 0.
		# This is going to show you how many trees have been analized


for t in os.listdir('%s' % treesFolder):		# Take each file in the folder that have the trees
	t = t.strip('\n')					
										
	if 'RAxML' in t and t.endswith('.tre'):						# Consider only the files that contain the word 'RAxML'.
		OG5 = 'OG5_' + t.split('OG5_')[-1][:6]	# Take the OG code from the file that contains the tree
										
		n = n+1							# Each time that you take a file with a tree. This counter add 1...
		print("walk_tree_contamination_single.py:\t%s\t%s" % (n,t))  		# ... and shows you in the terminal.
										# In this way you can track how many trees have been analyzed
		
		var.trees = []								# Here we are initializing the variable (for p4)...
													# ... that contains the tree
		correct_tree(t)								# Correcting the format of the tree, using the function...
													# correct_tree. Read the coments of "def correct_tree(t):"
		tree_file = '%stemporal_%s' % (treesFolder, t)	# Reading the CORRECTED file and taking the tree for python
		read(tree_file) 							# Reading the tree for p4 
		tree = var.trees[0] 						# The tree	itself in var 'tree' - going forward using p4
													# As the tree is now in 'tree', remove the temporal file
		os.remove('%stemporal_%s' % (treesFolder, t))		# After putting , remove that file
		
		taxon_names = tree.getAllLeafNames(0) 		# node 0 is the root, so this should list all taxon names
													# But, they will be like Sr_ap_pfal_PF110062_OG5_126572
													# You need a list of modified names like Sr_ap_pfal.
		taxon_names_mod = []						# This var accumulates the modified taxa names PER TREE
		taxon_names_mod2 = []						# Then, you should initialize this var here.
				
		for taxon in taxon_names:							# Take each taxon name
			taxon_mod = get_clades(taxon)[4]				# Modify the name using function 'get_clades' 
			taxon_names_mod.append(taxon_mod)				# Append the modified name in taxon_names_mod
			
			if taxlist == 'no':								# Next 3 lines are important for running without a taxa list
				if taxon_mod not in weirdtaxalist:
					weirdtaxalist.append(taxon_mod)
					
			taxon_names_mod2.append(get_clades(taxon)[3])	# Read at the top what 'get_clades' does 
															# append the modified name in taxon_names_mod
			
		taxon_names_mod2 = list(set(taxon_names_mod2))

		doBL = 'y'
		branches = {}
		numberNodes = 0
		numberLeaves = 0
		error_tree = ''	

		if len(taxon_names_mod2) > 1:					# We are only considering trees with at least 4 minor clades 
			for taxa in weirdtaxalist:				# For each weird taxon ...
				if taxa in taxon_names_mod:			# If the weird taxon is in taxon_names_mod...
													# As taxon_names_mod contains all taxa in the tree...
													# This conditional means that the weird taxon is in the tree
													# Create hashes for Ba/Za, Op, Pl, Am, Ex and Sr.
													# sister taxa and their size are appended in the hashes
													# This info changes per node. Hashes should be initialized 
					sizes_cladesBaZa = {}			
					sizes_cladesOp = {}				
					sizes_cladesPl = {}				
					sizes_cladesAm = {}				
					sizes_cladesEx = {}				
					sizes_cladesSr = {}				
				
					for node in tree.iterNodesNoRoot():				# For each node in the tree
						allTaxa_node = tree.getAllLeafNames(node)	# Get all taxa from leaves per node...
																	# ... and make a list.
						
						if doBL == 'y':
							numberNodes += 1
							branch = node.br.len											
							branches[str(node)] = float(branch)
							
						MC_list =[]

						for taxon in allTaxa_node:					# Take each taxon from the list above 
							MC_list.append(get_clades(taxon)[0])	# and get the major clade using 'get_clades'
																	# it will end up with ...
																	# list of MC of the taxa of each node
																	# The list above will be in MC_list.
							
						# In the next bunch of lines, the script determines if MC_list contains only Ba/za,
						# Op, Pl, Am, Ex or Sr. This means that the script search the nodes that only
						# contains an specific MC. Finally, the script appends the nodes of the clades 
						# and their sizes in the hashes. Also, the script allows one contaminant. For instance,
						# A clade with 9 Op leaves and 1 Pl leaf will be recorder as node_XXXXX -> 10.
						# This step is very important for the rooting.
							
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
						
						# In the next bunch of lines the tree is re-rooted by the biggest Ba/Za clade, or
						# the biggest Op clade, or the biggest Pl clade, or the biggest Am clade, or the
						# biggest Am clade, or the biggest Ex clade, or the biggest Ex clade.
						# In theory, only the first three options should occur. If the rooting fails,
						# print error tree in the terminal.
						
					if sizes_cladesBaZa: 		# If size_cladesBaZa is not 0
												# Take the node with the biggest number of taxa (leaves)
												# Re-root by this node
											
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
											error_tree = 'error tree'
				
					# The next loop iterates by the nodes of the tree
					# considering leaves as well as nodes
				
					if doBL == 'y':
						total_av = 0
					
						for i in branches: total_av = total_av + branches[i]
						averageBL = (total_av / numberNodes)				

					
					for node in tree.iterNodesNoRoot(): 				# but I really want to start with leaves!
						branchFinal = node.br.len
					
						if node.getNChildren() == 0: 			# If node is a leaf...						
							taxon_full = node.name					# Take the taxon that is in the leave 
							clade_name = get_clades(taxon_full)[4]	# Put taxon in the format MC_mc_sp
																# Ex: Sr_ap_pfal
																# save taxon in var clade_name
							taxa_interest = taxa				
							if taxa_interest in clade_name:		# If the leaf is the weird taxon...
																
								sisterTaxa = tree.getAllLeafNames(node.parent)	# Take all taxa that are in 
																				# the parental node of the leaf
																				# This list contains the weird
																				# taxon, and the sister taxa.
							
								sisterTaxa_real = []
								sisterTaxa_match = []
								taxon8times = 0      							# used later for checking monophyletic clades with DIFFERENT cells
								allowedStrains = ["1","2","3","4","5","6","7","8","9","0"]
							
								for taxon in sisterTaxa:						# used later for checking monophyletic clades SAME cells
									if taxa_interest in taxon: sisterTaxa_match.append(taxon)
							
								# In the next loop we will evaluate if it is monophyletic clades with DIFFERENT cells
							
								for taxon in sisterTaxa:								
									taxonfull = (get_clades(taxon)[4])
									if taxonfull[:-2] in taxa_interest:
										if taxonfull[-1] in allowedStrains:
											if taxonfull[8] in allowedStrains:
												taxon8times += 1									
											elif taxonfull[8] == taxa_interest[8]:
												taxon8times += 1
										
								# if monophyletic (with different cells) then we need to move to a deeper node
							
								if taxon8times == len(sisterTaxa):
									branch2report = (node.parent).br.len
									sisterNode = (node.parent).sibling
									if not sisterNode : sisterNode = (node.parent).leftSibling()
									sisterTaxa_real = tree.getAllLeafNames(sisterNode)
									if not sisterTaxa_real :   # Here we had to crete this conditional becase the step above didn't work when branches are 0																					
										sisterTaxa_real = tree.getAllLeafNames(sisterNode.parent)	
										for i in sisterTaxa : sisterTaxa_real.remove(i)
									
								else:
																
									if len(sisterTaxa_match) == len(sisterTaxa):	# So, if we also have a monophyletic clade (not from different cells), 
																					# then it should be processed as above.
									
										if not (node.parent).br : continue
										branch2report = (node.parent).br.len
										sisterNode = (node.parent).sibling
										if not sisterNode : sisterNode = (node.parent).leftSibling()
										sisterTaxa_real = tree.getAllLeafNames(sisterNode)
										if not sisterTaxa_real : 
											sisterTaxa_real = tree.getAllLeafNames(sisterNode.parent)	
											for i in sisterTaxa : sisterTaxa_real.remove(i)								
									
									else:										# when no mpnophyletic
										branch2report = branchFinal
										sisterTaxa_real = sisterTaxa
										sisterTaxa_real.remove(taxon_full)		# Remove leaf to take sister taxa
															
	#							sisterClades = []
								sisterMinors = []
								sisterSequences = ''								
								for taxon in sisterTaxa_real:		
	#								sisterClades.append(get_clades(taxon)[0])	# Modify the sister taxa
																				# from MC_mc_sp to MC
																				# Ex: Sr_ap_pfal to Sr
									sisterMinors.append(get_clades(taxon)[3])
									sisterSequences = sisterSequences + ',' + taxon
								
								# The next two lines take the list of MC of the sister taxa and remove duplicates
								# Ex: if the list were: [Sr, Sr, Sr, Am], it will produce [Sr, Am].
								# The example above is a non-monophyletic sister clade that produces 2 elements.
								# A monophyletic sister clade produces only one element in the list.
								# Ex: [Sr, Sr, Sr], produces [Sr]
								# The new list w/o duplicates is saved in var sisterClades
							
	#							sisterClades = set(sisterClades)	# 'set' makes an iterable but not-indexable object
	#							sisterClades = list(sisterClades)	# 'list' makes it indixeble
								sisterMinors = set(sisterMinors)
								sisterMinors = list(sisterMinors)

							
								minors = ''
								for minor in sisterMinors:
									minors = minors + ',' + minor
								
							
								if len(sisterMinors) == 1:							# If sisterMinors is 1 element
									minorClade = sisterMinors[0]					# get its minor clade

									if minorClade in taxa_interest:					# If weird taxon also has this MC
										result = 'same_minor'							# return 'same_MC'
										report.write (OG5 + '\t' + taxa + '\t' + taxon_full.replace('LKH', '-LKH') + '\t' +  result + '\tNA' + '\t' + str(branch2report) + '\t' + str(averageBL) + '\n')	
								
									else:											# If weird taxon hasn't this MC 	
										result = minorClade							# return the MC
																					# report result in output
																				
										if len(sisterTaxa_real) > 10 : sisterSequences = 'too-long'	
										report.write (OG5 + '\t' + taxa + '\t'+ taxon_full.replace('LKH', '-LKH') + '\t' +  result + '\t' + sisterSequences + '\t' + str(branch2report) + '\t' + str(averageBL) + '\n')	

								if len(sisterMinors) > 1:							# If sister clades has more than 1
									result = "non-monophyletic"						# retrieve 'non-monophyletic'
																					# report result in output
								
									if len(sisterMinors) > 20 : minors = 'too-long'
									report.write (OG5 + '\t' + taxa + '\t' + taxon_full.replace('LKH', '-LKH') + '\t' +  result + '\t'  + minors + '\t' + str(branch2report) + '\t' + str(averageBL) + '\n')
									
					doBL = 'n'

				else:
					result = 'no_taxaOFinterest' # If the tree does not have the weird taxon, return 'no_taxaOFinterest'
	#				report.write (OG5 + ',' + taxa + ',' +  result + '\n') # report result in output

			if 'error tree' in error_tree: # if the tree couldn't be re-rooted, retrieve 'OG cannot be annalized' 
				print("walk_tree_contamination.py: " + OG5 + ' cannot be annalized')
				error_tree = ''
		
		else:
			print("walk_tree_contamination.py: Tree is ignored because it contains fewer than 2 minor clades --> " + OG5)

report.close() 

sisters_summary = open('SisterSummary.csv', 'w')

#This summarizes the huge table initially generated, which isn't very human friendly
if(summarize_sisters):

	query_clades = [clade.strip() for clade in query_clades.split(',')]
	sister_clades = [clade.strip() for clade in sister_clades.split(',')]
	break_up_clades = [clade.strip() for clade in break_up_clades.split(',')]
	
	#Reading the output file from above
	report = open('sisterReport', 'r').readlines()

	print('report', report)
	
	#A list of all the trees and taxa in the set of trees
	trees = []
	taxa = []

	for i in report:
		trees.append(i.split("\t")[0])	
		taxa.append(i.split("\t")[1])
		
	trees = sorted(list(set(trees)))
	taxa = sorted(list(set(taxa)))
	
	print('trees', trees)
	print('taxa', taxa)
			
	#A list of all possible sister taxa to fill out later
	sisters_list = []
	#If list is constrained by user
	if(len([val for val in sister_clades if val.lower() != 'na']) > 0):
		for taxon in taxa:
			if(taxon[:5] in sister_clades or taxon[:2] in sister_clades or 
			taxon[:4] in sister_clades or taxon[:8] in sister_clades or taxon[:10] in sister_clades or taxon[:7] in sister_clades):
				sisters_list.append(taxon)
	else:
		sisters_list = [taxon[:10] for taxon in taxa]

	print('sisters', sisters_list)
				
	#Sorting sisters alphabetically (see output)
	sisters_list = sorted(sisters_list)
	
	sisters_template = { taxon[:10] : 0 for taxon in sisters_list }
		
	maj_clades = { 'Sr' : 0, 'Pl' : 0, 'Op' : 0, 'Am' : 0, 'Ex' : 0, 'EE' : 0, 'Ba' : 0, 'Za' : 0 }
	final_clades = { }; clades_in_order = []
	sisters_summary.write('Taxon,same,')
	for clade in break_up_clades:
		if(clade not in maj_clades):
			print('\nA clade that you input to be broken up (with the --break_up argument) is not valid. Make sure it is (at least) one of Sr, Pl, Op, Am, Ex, EE, Ba, or Za. Defaulting to summarizing at the major clade level.\n')
	for clade in maj_clades:
		if(clade in break_up_clades):
			minors = list(dict.fromkeys([taxon[:5] for taxon in taxa if taxon[:2] == clade]))	
			for minor in minors:
				clades_in_order.append(minor)
		else:
			clades_in_order.append(clade)
	sisters_summary.write(','.join(clades_in_order))
	sisters_summary.write(',,non-monophyletic,Total,Prop. Same,,' + ','.join(sisters_list) + '\n')

	print('qclades', query_clades)
	
	#If the user input a list of taxa/clade identifiers to constrain the tips for which we look for sisters
	if(len([val for val in query_clades if val.lower() != 'na']) > 0):
		for taxon in taxa[::-1]:
			if(taxon[:5] not in query_clades and taxon[:2] not in query_clades and 
			taxon[:4] not in query_clades and taxon[:8] not in query_clades and taxon[:10] not in query_clades and taxon[:7] not in query_clades):
				taxa.remove(taxon)

	print('taxa2', taxa)
			
	#For every taxon that appears across all trees (and not filtered out by user in taxon list or --query_clades args)
	for taxon in taxa:
		nm = 0; sm = 0; to = 0
		
		if(taxon not in final_clades):
			final_clades.update({ taxon : { clade : 0 for clade in clades_in_order } })
		
		sister_pa = sisters_template.copy()
		
		#For every line in the large non-human-friendly spreadsheet
		for line in report:
			values = line.split('\t')
			#If the line is about the current taxon of interest
			if values[1] == taxon:
			
				#Adding to the total number of appearances of the taxon
				to += 1
				
				unique_sisters = []
				if((branch_length_filter and float(values[5]) < float(values[6])) or (not branch_length_filter)):
					for sister in [sis[:10] for sis in line.split('\t')[4].split(',') if sis.strip() != '']:
						if(sister in sister_pa):
							unique_sisters.append(sister)
							
				if(single_sister_only and len(dict.fromkeys(unique_sisters)) > 1):
					continue
			
				#If it has sisters outside of its minor clade
				if values[3] != 'same_minor':
					#Apply branch length filter (or not) and then fill out values depending on the taxonomy of the sisters
					if((branch_length_filter and float(values[5]) < float(values[6])) or (not branch_length_filter)):
					
						if values[3] == 'non-monophyletic' : nm += 1
						
						for clade in final_clades[taxon]:
							if clade in values[3]: final_clades[taxon][clade] += 1
					
				#If sisters are the same minor clade
				else : sm += 1
				
				#Recording the sisters
				if((branch_length_filter and float(values[5]) < float(values[6])) or (not branch_length_filter)):
					for sister in [sis[:10] for sis in line.split('\t')[4].split(',') if sis.strip() != '']:
						if(sister in sister_pa):
							sister_pa[sister] += 1
					
		#Writing all requested data to the summary spreadsheet				
		sisters_summary.write(taxon + ',' + str(sm) + ',')
		for clade in clades_in_order:
			sisters_summary.write(str(final_clades[taxon][clade]) + ',')
		sisters_summary.write(',')
		sisters_summary.write(str(nm) + ',' + str(to) + ',' + str(round(float(sm) / float(to), 3)))
		sisters_summary.write(',,')
		for sister in sisters_list:
			if(sister_pa[sister] == 0):
				sisters_summary.write(',')
			else:
				sisters_summary.write(str(sister_pa[sister]) + ',')
		sisters_summary.write('\n')

















