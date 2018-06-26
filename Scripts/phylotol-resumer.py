import os, sys
import Utilities

arg = sys.argv

if 2 <= len(arg) <= 3:
	path = arg[1] 
	if len(arg) == 3:
		mode = (sys.argv[-1]).strip()
		if not mode == 'nr' : mode = 'df'
	else:
		mode = 'df'
else:
	print "Requires at least two arguments ...\n python phylotol-resumer.py path_to_pre-guidance_files mode(nr = no raxml, optional)"
	quit()

print "\n** mode -> %s **" % mode

# read parameters file and take the path to fasta files and the name of the OG list. 
infile = open('pipeline_parameter_file.txt','r').readlines()
for line in infile:
	if line[0] == '#':
		attribute = line.split()[0].strip('#')
		value = line.split()[2].strip()
		if attribute == 'PathtoFiles':
			PathtoFiles = value
		elif attribute == 'testPipelineList':
			oglistName = value
		elif attribute == 'guidanceIter':
			guidanceIter = value
		elif attribute == 'seqcutoff':
			seqcutoff = float(value)			
		elif attribute == 'colcutoff':
			colcutoff = float(value)
		elif attribute == 'rescutoff':
			rescutoff = float(value)

# make directories
PathtoOutput = os.path.abspath(path) + "/out_resume"
if os.path.exists(PathtoOutput):
	print "The folder " + PathtoOutput + " exists. Choose another path\n\n"
	quit()
else:
	os.system('mkdir ' + PathtoOutput)
	os.system('mkdir ' + PathtoOutput + "/Guidance")

oglist = open(PathtoFiles + '/' + oglistName,'r').readlines()  #list of ogs of interest
if oglist == []: 
	print 'terminating PhyloToL-resulmer: Your list of OGs is empty\n\n'
	quit()

for og in oglist:
	og = og.replace("\n", "")
	os.system("cp " + path + "/*" + og + "* " + PathtoOutput + '/Guidance/'+ og + 'forGuidance.fas')
	path2og = PathtoOutput + '/Guidance/'+ og + 'forGuidance.fas'
	og = og.replace("\n", "")
	Utilities.iterGuidance(oglistName, og, PathtoOutput, guidanceIter, seqcutoff, colcutoff, rescutoff, mode)

print "find your output here: " + PathtoOutput + "/" + oglistName + '_results2keep/'