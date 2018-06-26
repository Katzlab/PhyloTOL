#taxonomyFile  = File.open("taxonomyNOsps2.txt", 'r').readlines()
taxonomyFile = File.open("/Users/marioceron/Downloads/ott/taxonomy.tsv", "r").readlines()
#taxa2search = ["Entamoeba nuttalli", "Mastigamoeba balamuthi", "Acanthamoeba castellanii Neff", "Mayorella sp BSH02190019", "Neoparamoeba aestuarina", "Flamella fluviatilis", "Filamoeba nolandi ATCC50430", "Dictyostelium purpureum", "Physarum polycephalum", "Nolandella abertawensis", "Ancyromonas sigmoides", "Fabomonas tropica", "Mantamonas plastica", "Nutomonas longa", "Rigifila ramosa", "Thecamonas trahens", "Breviata anathema", "Subulatomonas tetraspora", "Choanocystis sp. FB-2015", "Raphidiophrys heterophryoidea", "Geminigera cryophila CCMP2564", "Emiliania huxleyi", "Isochrysis galbana", "Diphylleia rotans", "Tsukubamonas globosa", "Trimastix marina (NEW)", "Telonema subtile", "Roombia truncata", "Bodo saltans", "Diplonema papillatum [DP]", "Eutreptiella gymnastica NIES381", "Leishmania infantum JPCM5", "Trypanosoma congolense", "Spironucleus barkhanus", "Naegleria gruberi strain NEG-M", "Trimastix pyriformis ATCC 50562", "Reclinomonas americana", "Malawimonas jakobiformis", "Monocercomonoides sp.", "Tritrichomonas foetus", "Trichomonas vaginalis G3", "Monosiga brevicollis M", "Salpingoeca rosetta", "Anncaliia algerae PRA109", "Allomyces macrogynus GCA 000151295.1 A macrogynus V3", "Batrachochytrium dendrobatidis [BD]", "Conidiobolus coronatus", "Dacryopinax sp", "Lichtheimia corymbifera", "Neocallimastix patriciarum [NP]", "Rozella allomycis", "Saccharomyces cerevisiae S288c", "Amoebidium parasiticum", "Capsaspora owczarzaki", "Sphaeroforma arctica", "Ministeria vibrans", "Caenorhabditis elegans", "Carteriospongia foliascens", "Convoluta pulchra", "Capitella teleta", "Daphnia pulex", "Homo sapiens", "Hydra vulgaris", "Pleurobrachia pileus", "Saccoglossus kowalevskii", "Schistosoma mansoni", "Trichoplax adhaerens", "Fonticula alba", "Cyanophora paradoxa [Durnford lab] [CD]", "Glaucocystis nostochinearum", "Arabidopsis thaliana", "Amborella trichopoda", "Chlamydomonas reinhardtii", "Chlorella variabilis", "Marchantia polymorpha", "Prasinoderma coloniale CCMP1413", "Pterosperma sp CCMP1384", "Tetraselmis chuii", "Chondrus crispus", "Compsopogon coeruleus SAG 3694", "Galdieria sulphuraria", "Rhodosorus marinus CCMP 769", "Cryptosporidium parvum Iowa II", "Gregarina niphandrodes", "Plasmodium falciparum 3D7", "Vitrella brassicaformis CCMP3155", "Paramecium tetraurelia strain d4-2", "Stentor coeruleus", "Stylonychia lemnae", "Azadinium spinosum 3D9", "Gymnodinium catenatum GC744", "Hematodinium sp. SG 2012", "Oxyrrhis marinia", "Symbiodinium microadriaticum", "Perkinsus marinus", "Aulacantha scolymantha", "Astrolonche sp.", "Brevimastigomonas motovehiculus", "Bigelowiella natans", "Corallomyxa tenera", "Euglypha rotunda CCAP 1520/1", "Leptophrys vorax", "Sorites sp", "Aureococcus anophagefferens CCMP1850", "Aurantiochytrium mangrovei", "Blastocystis hominis", "Bolidomonas pacifica RCC208", "Chrysocystis fragilis CCMP3189", "Cafeteria roebergensis E410", "Chattonella subsalsa CCMP2191", "Dictyocha speculum CCMP1381", "Ectocarpus siliculosus", "Extubocellulus spinifer CCMP396", "Grammatophora oceanica CCMP 410", "Nannochloropsis gaditana", "Ochromonas sp CCMP1393", "Phytophthora infestans", "Phaeomonas parva CCMP2877", "Phaeodactylum tricornutum", "Synchroma pusillum CCMP3072 - MMETSP1452", "Thalassiosira pseudonana CCMP1335"]
#taxa2search = ["Corynebacterium diphtheriae 241","Aquifex aeolicus VF5","Bacteroides fragilis 638R","Chloroherpeton thalassium ATCC 35110","Chlamydia muridarum Nigg","Chloroflexus aurantiacus J 10 fl","Akkermansia muciniphila ATCC BAA 835","Anabaena cylindrica PCC 7122","Oscillatoria nigroviridis PCC 7112","Pseudanabaena sp PCC 7367","Truepera radiovictrix DSM 17093","Dictyoglomus turgidum DSM 6724","Geobacillus kaustophilus HTA426","Oscillibacter valericigenes","Fusobacterium nucleatum ATCC 25586","Thermodesulfovibrio yellowstonii DSM 11347","Azospirillum brasilense Sp245","Rickettsia prowazekii str. Madrid E","Burkholderia pseudomallei 1710b","Variovorax paradoxus EPS","Desulfovibrio aespoeensis Aspo 2","Acinetobacter baumannii 1656 2","Escherichia coli str. K-12 substr. W3110","Planctomyces limnophilus DSM 3776","Sphaerochaeta pleomorpha Grapes","Acholeplasma laidlawii PG 8A","Thermotoga maritima MSB8","Candidatus Heimdallarchaeota archaeon AB 125","Lokiarchaeum sp GC14 75","Candidatus Odinarchaeota archaeon LCB 4","Candidatus Thorarchaeota archaeon SMTZ1-45","miscellaneous Crenarchaeota group-1 archaeon SG8-32-1","Sulfolobus islandicus HVE10 4","Thermoproteus neutrophilus V24Sta","Methanobacterium SWAN 1","Methanocaldococcus infernus ME","Halarchaeum acidiphilum 489138","Natronobacterium gregoryi SP2 797304","Methanomethylovorans hollandica DSM 15978 867904","Thermoplasma volcanium GSS1","Thermococcus kodakarensis KOD1","Methanopyrus kandleri AV19","Candidatus Korarchaeum cryptofilum OPF8","Nanoarchaeum equitans Kin4-M","Candidatus Micrarchaeum acidiphilum ARMAN 2 425595","Cenarchaeum symbiosum A"]
taxa2search = ["Lokiarchaeum sp GC14 75","Sulfolobus islandicus HVE10 4","Thermoproteus neutrophilus V24Sta","Methanobacterium SWAN 1","Methanocaldococcus infernus ME","Halarchaeum acidiphilum 489138","Natronobacterium gregoryi SP2 797304","Methanomethylovorans hollandica DSM 15978 867904","Thermoplasma volcanium GSS1","Thermococcus kodakarensis KOD1","Methanopyrus kandleri AV19","Nanoarchaeum equitans Kin4-M","Cenarchaeum symbiosum A"]

puts "\nCollecting genera ...\n"

taxonomyNG = Array.new
correctedTaxa = Array.new
taxonomyStart = Hash.new
taxonomies = Array.new
genera = Array.new

taxa2search.each do |taxon2search|
	taxon2search = taxon2search.gsub("candidatus ", "")
	taxon2search = taxon2search.gsub("_", " ")
	taxon2search = taxon2search.gsub(/ .*$/, "")
	correctedTaxa << taxon2search
end

taxonomyFile.each do |line|
	lineSplitted = line.split("|")
	currentNode = lineSplitted[0].gsub("\t", "")
	parentNode = lineSplitted[1].gsub("\t", "")
	taxonomy = lineSplitted[2].gsub("\t", "")
	level = lineSplitted[3].gsub("\t", "")

	if level == "genus"
		if correctedTaxa.include? taxonomy
			if taxonomyStart[taxonomy]
				taxonomyStart[taxonomy] += "#"+ parentNode
				puts taxonomyStart[taxonomy]
			else
				taxonomyStart[taxonomy] = level + '=' + taxonomy + "#"+ parentNode
				puts taxonomyStart[taxonomy]
			end
		end
		
		genera << taxonomy + "," + currentNode + "," + parentNode
	else taxonomyNG << currentNode + "|" + parentNode + "|" + taxonomy + "|" + level end
end

puts "\nCollecting full taxonomies ...\n"
	
correctedTaxa.each do |taxon2search|
	nextNodes = taxonomyStart[taxon2search].split("#")[1..-1]
	nextNodes = nextNodes.uniq
	
	count = 0
	nextNodes.each do |nextNode|
		count += 1
		fullTaxonomy = taxonomyStart[taxon2search].split("#")[0] + ";"
		stop = 'no'
		status = 'running'

		while stop == 'no'
			findNode = 'no'
		
			taxonomyNG.each do |line|
				lineSplitted = line.split("|")
				currentNode = lineSplitted[0]
				parentNode = lineSplitted[1]
				taxonomy = lineSplitted[2]
				level = lineSplitted[3]
			
				if currentNode == nextNode
					findNode = 'yes'
					fullTaxonomy = fullTaxonomy + level + '=' + taxonomy + ";"
					nextNode = parentNode
				
					if nextNode == '' then status = 'finished' end
				
				end
			
				if status == 'finished' then stop = 'yes' end
			
			end
		
			if findNode == 'no'
				genera.each do |genus|
					genusCurrent = genus.split(",")[1]
					genusParent = genus.split(",")[2]
					genus = genus.split(",")[0]
				
					if genusCurrent == nextNode
						findNode = 'yes'
						fullTaxonomy = fullTaxonomy + 'genus' + '=' + genus + ";"
						nextNode = genusParent
					end
				end
			end
		
			if findNode == 'no' then stop = 'yes' end
		
		end

		taxonomies << taxon2search + " --> 4" + fullTaxonomy
		puts fullTaxonomy	

	end	
end

puts "\n****  Full taxonomies  ****\n"

taxonomies.each do |taxonomy| puts taxonomy end
