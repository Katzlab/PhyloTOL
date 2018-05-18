orthomclDB = File.open("OrthoMCL.fasta", "r").readlines()
ogs = File.open("oglistTest", "r").readlines()
`mkdir -p Data/OGs`

tag = ''
seq = ''
line = ''
orthomclHash = Hash.new

orthomclDB.each do |line|
	line = line.gsub("\n", "")
	
	if line =~ />/
		puts line
		if not seq == "" then orthomclHash[tag] = seq end
		tag = line
		seq = ""
	end
	
	if line =~ /(^[A-z])|(^\*)/
		seq = seq + line
	end
end 

orthomclHash[tag] = seq

ogs.each do |og|
	og = og.gsub("\n", "")
	ogFile = File.open("Data/OGs/" + og, "w")
	puts og + "\n"
	
	orthomclHash.each do |key, value|
		if key.include? og
			ogFile.write(key + "\n" + value + "\n")
		end
	end
	ogFile.close
end
