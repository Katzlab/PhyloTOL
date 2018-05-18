use strict;
use warnings;

my $orthomclFile = open(ODB, "OrthoMCL.fasta") or die "the file OrthoMCL.fasta doesn't exist\n";
my @orthomclDB = <ODB>;
my $path2data = "/data1/mario/addGenes/Data/";
my $path2temp = "$path2data/Temp/";
system("mkdir -p $path2data/OGs");
system("mkdir -p $path2data/OGs_filtered");
system("mkdir -p $path2temp/OGs4guidance");
system("mkdir -p $path2temp/names");
system("mkdir -p $path2temp/OGs_filteredNoName");

my $oglistFile = open(OGL, "oglist") or die "the file oglist does not exist\n";
my @ogs = <OGL>;

my $tag = '';
my $seq = '';
my $line = '';
my %orthomclHash;
my %names;

foreach $line (@orthomclDB){
	chomp $line;
	
	if ($line =~ />/){
		unless ($seq eq ''){$orthomclHash{$tag} = $seq}
		$tag = $line;
		$seq = '';
	} 
	
	if ($line =~ /(^[A-z])|(^\*)/) {
		$seq = $seq . $line;
	}
}

$orthomclHash{$tag} = $seq;

foreach my $og (@ogs){
	chomp $og;
	my $count = 0;
	my $ogFile = open(OG, ">>$path2data/OGs/$og") or die "the file $og doesn't exist\n";
	my $ogFile2 = open(OG2, ">>$path2temp/OGs4guidance/$og") or die "the file $og doesn't exist\n";
	my $name = open(NAME, ">>$path2temp/names/$og") or die "the file $og doesn't exist\n";
	print $og . "\n";
	
	foreach my $key (keys %orthomclHash) {
		if ($key =~ $og) {
			$count += 1;
			print OG "$key\n$orthomclHash{$key}\n";
			print OG2 ">$count\n$orthomclHash{$key}\n";
			print NAME "$key\t>$count\n";
			$names{">$count"} = $key;
		}
	}
	
	close OG;
	close OG2;
	close NAME;
	
	system("bash exeGuidance.sh -i $path2data/Temp/OGs4guidance/$og -o $path2data/Temp/ -t 7 -c ./01-Guidance_out -s 0.4");
	
	my $toRenameFile = open(TRF, "$path2temp/OGs_filteredNoName/$og");
	my @toRename = <TRF>;
	my $renamedFile = open(RFL, ">>$path2data/OGs_filtered/$og");
	
	foreach my $toRenameLine (@toRename){
		if ($toRenameLine =~ /^>/) {
			chomp $toRenameLine;
			my $tag4rename = $names{$toRenameLine};
			print RFL "$names{$toRenameLine}\n$orthomclHash{$tag4rename}\n";
		}
	}
	close RFL;
}
