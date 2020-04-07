# This script removes sequences according to specified rules in local databades (in each instance of the pipeline). 
# It also produces the list of removed sequences so that they can be used for removing in the ready to go folder. 

# input: 
# - Report of sister taxa
# - rules
# - folder of ncbiFiles
# - empty folder for new ncbiFiles

# Running:
# - Put all the imput files and folders in the same folder and run scrip with ruby: "ruby seqs2remove.rb"

# ruby seqs2remove.rb path2Files sisterReport rules seqs2remove_out nonHomologs
# ruby seqs2remove.rb ./ sisterReport rules seqs2remove_out nonHomologs

path = ARGV[0]
report_summary = File.open(ARGV[1], 'r').readlines()
rules = File.open(ARGV[2], 'r').readlines()
sequences_contamination = File.open(ARGV[3], 'w')
trees_contamination = File.open(ARGV[3] + '_treesWcont', 'w')
ncbiFiles = path + '/ncbiFiles/'
system "mkdir " + path + "/newncbi/"
newncbi = path + '/newncbi/'
seqs2remove = Array.new
nonhomologs = File.open(ARGV[4], 'r').readlines()
total2remove = nonhomologs


count = 0
treesWcontamination = Array.new
report_summary.each do |line|
	count += 1
	line = line.chomp
	line = line.split("\t")
	og = line[0]
	taxon = line[1]
	sequence = line[2]
	sister = line[3]

	rules.each do |rule|
		rule = rule.chomp
		rule = rule.split("\t")
		taxon_rule = rule[0]
		
		if taxon == taxon_rule
			contamination = rule[1..-1]	
			contamination.each do |taxon_contamination|

				if sister.include? taxon_contamination
					puts "seqs2remove.rb : contamination --> " + sequence
					seqs2remove << sequence
					treesWcontamination << og
					sequences_contamination.write(sequence + "\n")
				end
			end	
		end
	end
end

trees2reprocess = treesWcontamination.uniq
if treesWcontamination != []
	trees2reprocess.each do |tree|
		puts "seqs2remove.rb: tree with contamination --> " + tree
		trees_contamination.write(tree + "\n")
	end
end

(total2remove << seqs2remove).flatten!

puts "\nseqs2remove.rb: removing sequences from Databases:"
(Dir.open(ncbiFiles)).each do |ncbiFile|
	if ncbiFile.include? ".fasta"
		taxon = ncbiFile[0..9]
		to_remove = Array.new
		
		ncbisequences_raw = File.open(ncbiFiles + ncbiFile, "r").readlines()
		
		# ---- correcting .fasta format ----
		
		ncbisequences = Array.new
		seq = ''
		ncbisequences_raw.each do |line|
			if line =~ /^>/
				if seq != '' then ncbisequences << seq end	
				ncbisequences << line.gsub("\n", "")
				seq = ''
			end
			if line =~ /^([A-z]|\*)/ then seq = seq + line.gsub("\n", "") end
		end
		ncbisequences << seq
		
		#------------------------------------
		
		newncbiTags = Array.new
		newncbiFile = File.open(newncbi + ncbiFile, "w")

		puts ncbiFile
	
#		seqs2remove.each do |seq2remove|
		total2remove.each do |seq2remove|
			seq2remove = seq2remove.gsub(/\n/, "")
			if seq2remove.include? taxon
				to_remove << seq2remove
			end
		end
		
		index = 0
		ncbisequences.each do |ncbisequence|
			if ncbisequence =~ /^>/
				tag = ncbisequence.gsub(/>|\n/, "")
				unless to_remove.include? tag
					newncbiFile.write(ncbisequence + "\n" + ncbisequences[index + 1] + "\n")
				else
					puts "\nremoved: " + ncbisequence
					puts "removed: " + ncbisequences[index + 1] + "\n"
				end
			end
			index += 1
		end
		newncbiFile.close
	end
end

system "rm -r " + ncbiFiles
system "mv " + newncbi + "\t" + ncbiFiles
