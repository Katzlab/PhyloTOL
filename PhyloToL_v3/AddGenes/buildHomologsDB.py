import os,re

orthomclDB = open("OrthoMCL.fasta", "r").readlines()
ogs = open("oglistTest", "r").readlines()
os.system("mkdir -p Data/OGs")

tag = ""
seq = ""
line = ""
orthomclHash = {}

for line in orthomclDB:
	line = line.replace("\n", "")
	
	if ">" in line:
		print line
		if not seq == "" : orthomclHash[tag] = seq
		tag = line
		seq = ""
	
	if re.match("(^[A-z])|(^\*)", line):
		seq = seq + line

orthomclHash[tag] = seq

for og in ogs:
	og = og.replace("\n", "")
	ogFile = open("Data/OGs/" + og, "w")
	print og + "\n"
	
	for key in orthomclHash:
		if og in key:
			ogFile.write(key + "\n" + orthomclHash[key] + "\n")
	ogFile.close()