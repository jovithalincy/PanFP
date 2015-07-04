#!/usr/bin/env perl

# usage:
#
# perl bin/lineage_func.pl -i <full_path_to/otu_table_modified.txt> -d <full_path_to/org_chr_16SrRNA_KO.TH30.txt> -a <annotation/KO> -o <full_path_to/lineage_func>
#
# written by: Se-Ran Jun

use warnings;
use strict;

use Getopt::Std;

my $undefined = 'undef';
my $funcDB = 'KO'; # functional annotation of interest

my %args;

# read arguments ...
getopt('i:f:d:a:o', \%args);

my $otu_table = $args{i}; # input: an otu table with modified lineage
my $db_file = $args{d}; # input: a list of organisms with functional annotations
my $annotation_dir = $args{a}; # input: a directory which includes functional annotations of organisms
my $output_dir = $args{o}; # output directory

die "usage: bin/lineage_func.pl -i <full_path_to/otu_table_modified.txt> -d <full_path_to/org_chr_16SrRNA_KO.TH30.txt> -a <annotation/KO> -o <full_path_to/lineage_func>\n"
unless defined ($otu_table) && defined ($db_file) && defined ($annotation_dir) && defined ($output_dir);

# check whether outpur directory exists. If not, create the directory
if ( !(-d $output_dir) ) {
	system "mkdir $output_dir";
}

# read a list of lineages present in a modified otu-table
my %lineages = ();
open (IN, "<$otu_table") || die "$!";
while (<IN>) {
	
	chomp $_;
	next if ($_ =~ /^#/);
	$_ = trimwhitespace($_);
	my @line = (); @line = split('\t', $_);

	my $lineage = $line[-1];
	$lineages{$lineage} = 1;
}
close IN;
my $nlineage = scalar keys %lineages;
print "the number of lineages = $nlineage\n";

# read organisms with 16S copy number and full lineage for a functional database of interest
my %organism_lineage = ();
open ( IN, "<$db_file" ) || die "$!";
while (<IN>) {

	chomp $_;
	next if ($_ =~ /^#/);
	$_ = trimwhitespace($_);
	my @line = (); @line = split('\t', $_);

	my $organism = $line[0];
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

	$organism_lineage{$organism} = $lineage;
}
close IN;

# make a functional profile for each lineage by collecting functional profiles of organisms belonging to the lineage
foreach my $taxon ( keys %lineages ) {

	my %lineage_func = (); my $lineage_size = 0;
	foreach my $organism ( keys %organism_lineage ) {

		my $lineage = $organism_lineage{$organism};
		if ( $lineage =~ /$taxon/ ) {

			$lineage_size++; # the number of organisms pooled for the taxon

			my $infile = "$organism.$funcDB";
			open (IN, "<$annotation_dir/$infile") || die "$!";
			while (<IN>) {

				chomp $_;
				next if ($_ =~ /^#/);
				$_ = trimwhitespace($_);
				my @line = (); @line = split('\t', $_);

				for my $i ( 1 .. $#line ) { # more than one functional term can be assigned
					my $func = $line[$i];
					if(defined $lineage_func{$func}){
						$lineage_func{$func}++;
					} else {
						$lineage_func{$func} = 1;
					}
				}
			}
			close IN;
		}
	}

	my $retaxon = $taxon;
	$retaxon =~ s/:/__/g;
     $retaxon =~ s/;/_/g;
	my $outfile = "$retaxon.$funcDB";

	# print a functional profile of taxon present in a modified otu-table
	open ( OUT, ">$output_dir/$outfile" ) || die "$!";
	print OUT "#$taxon\t$lineage_size\n";
	foreach my $func ( keys %lineage_func ) {
		print OUT "$func\t$lineage_func{$func}\n";
	}
	close OUT;
}

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
