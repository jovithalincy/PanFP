#!/usr/bin/env perl

# usage:
#
# perl bin/otu_table_norm_by_samplesize.pl -i <full_path_to/otu_table_modified_norm_by_copynum.txt> -o <full_path_to/otu_table_modified_norm_by_copynum_samplesize.txt>
#
# written by: Se-Ran Jun

use warnings;
use strict;

use Getopt::Std;

my %args;

# read arguments ...
getopt('i:o', \%args);

my $otu_table = $args{i}; # input: an otu table normalized by 16S copy number
my $normalized_otu_table = $args{o}; # output: an otu table normalized by sample size

die "usage: bin/otu_table_norm_by_samplesize.pl -i <full_path_to/otu_table_modified_norm_by_copynum.txt> -o <full_path_to/otu_table_modified_norm_by_copynum_samplesize.txt>\n"
unless defined ($otu_table) && defined ($normalized_otu_table);

# read a modified otu-table normalized by median 16S copy number of otus' lineage
my @samples = (); my @otus = (); my %otu_sample = (); my %otu_lineage = ();
open ( IN, "<$otu_table" ) || die "$!";
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

# calculate sample size
my %sample_size = ();
for my $sample ( @samples ) {
	$sample_size{$sample} = 0.0;
	for my $otu ( @otus ){
     	$sample_size{$sample} = $sample_size{$sample} + $otu_sample{$otu}{$sample};
     }
}

# print an otu-table normalized by sample size
open ( OUT, ">$normalized_otu_table" ) || die "$!";
print OUT "#OTU";
for my $sample ( @samples ) {
	print OUT "\t$sample";
}
print OUT "\tmodifiedConsensusLineage\n";
for my $otu ( @otus ) {

	print OUT "$otu";
     for my $sample ( @samples ) {
		my $normalized_otu_sample = $otu_sample{$otu}{$sample}/$sample_size{$sample};
          printf(OUT "\t%.10f", $normalized_otu_sample);
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
