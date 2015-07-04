#!/usr/bin/env perl

# usage:
#
# perl bin/modify_otu_table.pl -i <full_path_to/otu_table.txt> -d <full_path_to/org_chr_16SrRNA_KO.TH30.txt> -o <full_path_to/otu_table_modified.txt>
#
# written by: Se-Ran Jun

use strict;
use warnings;

use Getopt::Std;

#------------------------------------------------------------------
# make sure that the following words are used or something else
#------------------------------------------------------------------
my $nolineage = '__'; # if an otu has a lineage assigned, then the lineage has at least '__' in the string
my $undefined = 'undef';
my $candidatus = 'Candidatus';
#------------------------------------------------------------------

my %args;

# read arguments ...
getopt('i:d:o', \%args);

my $otu_table = $args{i}; # input: an otu table: tab delimited QIIME format
my $db_file = $args{d}; # input: a list of organisms with functional annotations stored in data directory
my $modified_otu_table = $args{o}; # output: an otu table with modified lineage

die "usage: bin/modify_otu_table.pl -i <full_path_to/otu_table.txt> -d <full_path_to/org_chr_16SrRNA_KO.TH30.txt> -o <full_path_to/otu_table_modified.txt>\n"
unless defined ($otu_table) && defined ($db_file) && defined ($modified_otu_table);

# read an otu-table
my @samples = (); my @otus = (); my %otu_sample = (); my %otu_lineage = ();
open (IN, "<$otu_table") || die "$!";
while (<IN>) {

	chomp $_;
	$_ = trimwhitespace($_);
	next if ( $_ =~ /# Constructed from biom file/);
	my @line = (); @line = split('\t', $_);

	my $lineage = $line[-1];

	if ($_ =~ /^#/) {
		for my $i ( 1 .. ($#line - 1) ) {
			my $sample = $line[$i];
			push(@samples, $sample);
		}
	} else {
		if( $lineage =~ /$nolineage/ ) { # exclude otu with no lineage information

			my $otu = $line[0];
			push (@otus, $otu);

			for my $i ( 1 .. ($#line - 1) ) {
				my $sample = $samples[$i - 1];
				$otu_sample{$otu}{$sample} = $line[$i];
			}

			my @rank_taxo = (); @rank_taxo = split('; ', $lineage);

			my $lineage = '';
			for my $i ( 0 .. $#rank_taxo ) {

				$rank_taxo[$i] = trimwhitespace($rank_taxo[$i]);

				( my $rank, my $taxo ) = split('__', $rank_taxo[$i]);

				# for example: [Clostridium]
				$taxo = $1 if($taxo =~ /\[(\w+)\]/); # parse out in [ ]

				# delete "Candidatus" since "Candidatus" was deleted in a db-file
				if ( $taxo =~ /^$candidatus/ ) {
					$taxo = substr( $taxo, length($candidatus), length($taxo) - length($candidatus) );
				}

				if( length($taxo) == 0 ) {
					$taxo = $undefined; # unify into 'undef'
				}

				# it should be up to s: for all otus
				$lineage = $lineage."$rank:$taxo;";
			}

			$lineage = substr($lineage, 0, length($lineage) - 1); # delete ';' at the last lineage
			$otu_lineage{$otu} = $lineage;
		} else {
			print "no lineage information : $_\n";
		}
	}
}
close IN;
my $notu = scalar @otus;
my $nsample = scalar @samples;

print "the number of otus = $notu\n";
print "the number of samples = $nsample\n";

# read organisms with functional annotations
my %db_lineages = ();
open (IN, "<$db_file") || die "$!";
while (<IN>) {

	chomp $_;
	next if ($_ =~ /^#/);
	$_ = trimwhitespace($_);
	my @line = (); @line = split('\t', $_);

	my $organism = $line[0];
	my $lineage = $line[-1];

	# it should be up to species
	$db_lineages{$lineage} = 1;
}
close IN;

# modify and trim otu's lineage that the resulting lineage agrees with NCBI taxonomy nomenclature
my %otu_modified_lineage = ();
foreach my $otu (keys %otu_lineage) {

	# all otus have the same length of lineage: k__;p__;c__;o__;f__;g__;s__
	my @rank_taxo = (); @rank_taxo = split(';', $otu_lineage{$otu});

	# check if lineage information at the taxonomic rank exists in NCBI taxonomy nomenclature
     my @rank_taxo_info = ();
	for my $i (0 .. $#rank_taxo) {
		if($rank_taxo[$i] !~ /$undefined/) {
			my $chk = 0;
			foreach my $lineage ( keys %db_lineages ) {
				if($lineage =~ /$rank_taxo[$i]/) {
					$chk = 1;
					last;
				}
			}
              	# 0: undefined in NCBI taxonomy nomenclature, 1: defined in NCBI taxonomy nomenclature
              	$rank_taxo_info[$i] = $chk;
		} else {
          	# 0: undefined in NCBI taxonomy nomenclature
               $rank_taxo_info[$i] = 0;
		}
	}

	# keep the last level (which is not 'undef') identified in taxonomy in db-file
     my $modified_lineage = ''; my $idx1 = -1;
     for (my $i = $#rank_taxo ; $i >= 0 ; $i-- ) {
     	if ( $rank_taxo_info[$i] == 1 ) {
			foreach my $db_lineage ( keys %db_lineages ) {
				if ( $db_lineage =~ /$rank_taxo[$i]/ ) {
					my $idx2 = index ($db_lineage, $rank_taxo[$i]);
					$modified_lineage = substr($db_lineage, 0, $idx2);
					last;
				}
			}
			$idx1 = $i;
               last;
		}
     }
	$modified_lineage = $modified_lineage."$rank_taxo[$idx1]";

	my $ext_modified_lineage = $modified_lineage;
	if ( $idx1 == 0 ) {
		$ext_modified_lineage = "$modified_lineage;p:undef;c:undef;o:undef;f:undef;g:undef;s:undef";
	} elsif ( $idx1 == 1 ) {
		$ext_modified_lineage = "$modified_lineage;c:undef;o:undef;f:undef;g:undef;s:undef";
	} elsif ( $idx1 == 2 ) {
		$ext_modified_lineage = "$modified_lineage;o:undef;f:undef;g:undef;s:undef";
	} elsif ( $idx1 == 3 ) {
		$ext_modified_lineage = "$modified_lineage;f:undef;g:undef;s:undef";
	} elsif ( $idx1 == 4 ) {
		$ext_modified_lineage = "$modified_lineage;g:undef;s:undef";
	} elsif ( $idx1 == 5 ) {
		$ext_modified_lineage = "$modified_lineage;s:undef";
	}

     # choose only otu with at least kingdom and phylum information
     if ( ($modified_lineage =~ /k:/ ) && ($modified_lineage =~ /p:/) ) {
		$otu_modified_lineage{$otu} = $modified_lineage;
     }

	if ( $otu_lineage{$otu} ne $ext_modified_lineage ) {
		#print "$otu_lineage{$otu}\t$ext_modified_lineage\n";
	}
}

print "the number of otus in an otu-table = ", scalar keys %otu_lineage,"\n";
print "the number of otus with functional annotation after modifying lineage = ", scalar keys %otu_modified_lineage,"\n";

# print an otu-table with modified lineages
open (OUT, ">$modified_otu_table") || die "$!";
print OUT "#OTU";
for my $sample ( @samples ) {
	print OUT "\t$sample";
}
print OUT "\tmodifiedConsensusLineage\n";
for my $otu ( @otus ) {
	if ( defined $otu_modified_lineage{$otu} ) {

		print OUT "$otu";
		for my $sample ( @samples ) {
			print OUT "\t$otu_sample{$otu}{$sample}";
		}
		print OUT "\t$otu_modified_lineage{$otu}\n";
	}
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
