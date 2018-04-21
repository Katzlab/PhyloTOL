use strict;
use warnings;

my $orthomclFile = open(ODB, "OrthoMCL.fasta") or die "the file OrthoMCL.fasta doesn't exist\n";
my @orthomclDB = <ODB>;

my $og = 'OG5_181689';

my $ogFile = open(OG, ">>$og") or die "the file $og doesn't exist\n";

my $copy = 'n';
my $tag = '';
my $seq = '';
my $line = '';
my %orthomclHash;

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

#close OG;

foreach my $key (keys %orthomclHash) {
	print $key . "\n";
	print $orthomclHash{$key} . "\n";
}

