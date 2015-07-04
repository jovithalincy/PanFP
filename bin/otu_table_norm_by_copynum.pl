#!/usr/bin/env perl

# usage:
#
# perl bin/otu_table_norm_by_copynum.pl -i <full_path_to/otu_table_modified.txt> -c <full_path_to/lineage_16Scopynum.txt> -o <full_path_to/otu_table_modified_norm_by_copynum.txt>
#
# written by: Se-Ran Jun

use warnings;
use strict;

use Getopt::Std;

my %args;

# read arguments ...
getopt('i:c:o', \%args);

my $otu_table = $args{i}; # input: an otu table with modified lineage
my $copynum_file = $args{c}; # input: a list of lineages present in the otu table with their 16S copy number distributions
my $normalized_otu_table = $args{o}; # output: an otu table normalized by the median 16S copy number of OTUs' lineages

die "usage: bin/otu_table_norm_by_copynum.pl -i <full_path_to/otu_table_modified.txt> -c <full_path_to/lineage_16Scopynum.txt> -o <full_path_to/otu_table_modified_norm_by_copynum.txt>\n"
unless defined ($otu_table) && defined ($copynum_file) && defined ($normalized_otu_table);

# read a modified otu-table
my @samples = (); my @otus = (); my %otu_sample = (); my %otu_lineage = ();
open (IN, "<$otu_table") || die "$!";
while (<IN>) {

     chomp $_;
	$_ = trimwhitespace($_);
	my @line = (); @line = split('\t', $_);

     if($_ =~ /^\#/){
     	@samples = @line[1 .. ($#line - 1)];
     } else {
		my $otu = $line[0];
		my $lineage = $line[-1];

		push(@otus, $otu);
		$otu_lineage{$otu} = $lineage;

		for my $i ( 1 .. ($#line - 1) ) {
			my $sample = $samples[$i - 1];
			$otu_sample{$otu}{$sample} = $line[$i];
		}
	}
}
close IN;
my $notu = scalar @otus;
my $nsample = scalar @samples;

print "the number of otus = $notu\n";
print "the number of samples = $nsample\n";

# read median 16S copy numbers for lineages present in the otu-table
my %lineage_copynum = ();
open ( IN, "<$copynum_file" ) || die "$!";
while (<IN>) {
	
	chomp $_;
	next if ($_ =~ /^#/);
	$_ = trimwhitespace($_);
	my @line = (); @line = split('\t', $_);

	my $lineage = $line[0];
	my $copynum = $line[1]; # take medium copynumber

	$lineage_copynum{$lineage} = $copynum;
}
close IN;

# normalize the otu-table by median 16S copy number of otus' lineages
my %otu_sample_norm_by_copynum = ();
for my $otu ( @otus ) {

	my $lineage = $otu_lineage{$otu};

	for my $sample ( @samples ) {
     	my $normalized_otu_sample = $otu_sample{$otu}{$sample}/$lineage_copynum{$lineage};
		$otu_sample_norm_by_copynum{$otu}{$sample} = $normalized_otu_sample;
	}
}

# print a modified otu-table normalized by median 16S copy number of otus' lineages
open ( OUT, ">$normalized_otu_table" ) || die "$!";
print OUT "#OTU";
for my $sample ( @samples ) {
        print OUT "\t$sample";
}
print OUT "\tmodifiedConsensusLineage\n";
for my $otu ( @otus ) {
	print OUT "$otu";
	for my $sample ( @samples ) {
     	printf(OUT "\t%.10f", $otu_sample_norm_by_copynum{$otu}{$sample});
     }
	print OUT "\t$otu_lineage{$otu}\n";
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
