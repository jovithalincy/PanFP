#!/usr/bin/perl

# usage:
#
# perl bin/lineage_16Scopynum.pl -i <full_path_to/otu_table_modified.txt> -d <full_path_to/org_chr_16SrRNA_KO.TH30.txt> -o <full_path_to/lineage_16Scopynum.txt>
#
# written by: Se-Ran Jun

use warnings;
use strict;

use Getopt::Std;
use Statistics::Basic qw(:all);

my $undefined = 'undef';

my %args;

# read arguments ...
getopt('i:d:o', \%args);

my $otu_table = $args{i}; # input: an otu table with modified lineage
my $db_file = $args{d}; # input: a list of organisms with functional annotations stored in data directory
my $copynum_dist = $args{o}; # output: a list of distributions of 16S copy numbers for lineages in the otu table

die "usage: bin/lineage_16Scopynum.pl -i <full_path_to/otu_table_modified.txt> -d <full_path_to/org_chr_16SrRNA_KO.TH30.txt> -o <full_path_to/lineage_16Scopynum.txt>\n"
unless defined ($otu_table) && defined ($db_file) && defined ($copynum_dist);

# read organisms with 16S copy number and full lineage for a functional database of interest
my %organism_copynum = (); my %organism_lineage = ();
open (IN, "<$db_file") || die "$!"; 
while (<IN>) {

	chomp $_;
	next if($_ =~ /^#/);
	$_ = trimwhitespace($_);
	my @line = (); @line = split('\t', $_);

	my $organism = $line[0];
	my $copynum = $line[2];
	my $lineage = $line[-1]; # it should have all ranks: k:;p:;c:;o:;f:;g:;s:

	# filter out "undef"
	my @rank_taxo = (); @rank_taxo = split(';', $lineage);
	if ( $rank_taxo[$#rank_taxo] =~ /$undefined/ ) {
		my $idx = -1;
		for ( my $i = $#rank_taxo ; $i > 0 ; $i-- ) {
			if ( ( $rank_taxo[$i] =~ /$undefined/ ) && ( $rank_taxo[$i-1] !~ /$undefined/ ) ) {
				$idx = $i-1;
				last;
			}
		}

		$lineage = join(';', @rank_taxo[0 .. $idx]);
	}

	$organism_copynum{$organism} = $copynum;
	$organism_lineage{$organism} = $lineage;
}
close IN;

# read taxons present in a modified otu-table
my %otu_lineages = ();
open (IN, "<$otu_table") || die "$!"; 
while(<IN>) {

	chomp $_;
	next if($_ =~ /^#/);
	$_ = trimwhitespace($_);
	my @line = (); @line = split('\t', $_);

	my $lineage = $line[-1];
	$otu_lineages{$lineage} = 1;
}
close IN;

# make distributions of 16S copy numbers of taxons from organisms belonging to taxons
my %lineage_copynum = ();
foreach my $lineage ( keys %otu_lineages) {

	@{$lineage_copynum{$lineage}} = ();

	foreach my $organism ( keys %organism_lineage ) {
		if ( $organism_lineage{$organism} =~ /$lineage/ ) {
			push( @{$lineage_copynum{$lineage}}, $organism_copynum{$organism} );
		}
	}
}

# print distribution of 16S copy numbers of lineages present in a modified otu-table
open( OUT, ">$copynum_dist" ) || die "$!";
print OUT "#lineage\tmedian16Scopynum\tavg16Scopynum\t16Scopynums\n";
foreach my $lineage ( keys %otu_lineages ) {

	my $median16S = median(@{$lineage_copynum{$lineage}});
	my $mean16S = mean(@{$lineage_copynum{$lineage}});

	printf (OUT "%s\t%d\t%.3f", $lineage, $median16S, $mean16S);

	for my $copynum ( @{$lineage_copynum{$lineage}} ) {
		print OUT "\t$copynum";
	}
	print OUT "\n";
}
close OUT;

exit 0;

###############################################################################################
# SUBROUTINE
###############################################################################################

# trimwhitespace : remove whitespace from the start and end of the string
sub trimwhitespace{

        use strict;
        use warnings;

        my $string=shift;

        $string=~s/^\s+//;
        $string=~s/\s+$//;

        return $string;
}
