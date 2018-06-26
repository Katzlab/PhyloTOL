# This script only work if all taxa have different genus

phyloFile = File.open("/Users/marioceron/Documents/katzlab/duplications/speciesTree/ottt/allphylogenies6", "r").readlines()
#phyloFile = File.open("/Users/marioceron/Downloads/ott/taxonomy.tsv", "r").readlines()
phylogenies = Hash.new
phylogeniesNew = Hash.new
num_levels = 8

# collecting genera

nd = 0
phyloFile.each do |phylogeny|
	code = phylogeny.split("\t")[0]
	phylogeny = (phylogeny.split("\t")[1]).gsub("\n", "")
	phylogenies[phylogeny] = code 
end

# Sorting phylogenies by level

for i in 1..(num_levels	-1)
	puts "\n--- Sorting phylgenies --- Level: " + i.to_s
	phylogenies.each do |key, value|	
		clade = key.split(";")[1]
		levelId = clade.split("=")[0]
		acualLevel = clade.split("=")[1]
	
		if not acualLevel
			nd += 1
			clade = levelId + "=" + "ND" + nd.to_s
		end

		key = clade + ";" + ((key.split(";")[2..-1]).join(";"))

		if not phylogeniesNew[key]
			phylogeniesNew[key] = value
		else
			phylogeniesNew[key] = "1\t" + phylogeniesNew[key] + "," + value
		end
	end

	phylogeniesNew.each do |key, value|
		if value =~ /^1\t/
			value = value.gsub("1\t", "")
			value = "(" + value + ")"
			phylogeniesNew[key] = value
		end
		puts "value: " + value + " key: " + key
	end

	phylogenies = phylogeniesNew
	phylogeniesNew = Hash.new
end

tree = ''
phylogenies.each do |key, value|
	tree += value + ","
end

puts "(" + tree.sub(/,$/, "") + ");"
