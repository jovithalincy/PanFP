#!/usr/bin/env perl

# usage:
#
# perl bin/func_sample.pl -i <full_path_to/lineage_sample.txt> -l <full_path_to_lineage_func> -o <full_path_to/func_sample.txt>
#
# written by: Se-Ran Jun

use warnings;
use strict;

use Getopt::Std;

my %args;
my $funcDB = 'KO'; # functional annotation of interest

# read arguments ...
getopt('i:f:l:o', \%args);

my $lineage_sample_table = $args{i}; # input: a lineage sample table
my $lineage_func_dir = $args{l}; # input: a directory which includes lineages' functional profiles
my $func_sample_table = $args{o}; # output: a function sample table

die "usage: bin/func_sample.pl -i <full_path_to/lineage_sample.txt> -l <full_path_to_lineage_func> -o <full_path_to/func_sample.txt>\n"
unless defined ($lineage_sample_table) && defined ($lineage_func_dir) && defined ($func_sample_table);

# read a lineage-sample table
my @samples = (); my @lineages = (); my %lineage_sample = ();
open (IN, "<$lineage_sample_table" ) || die "$!";
while (<IN>) {

	chomp $_;
	$_ = trimwhitespace($_);
	my @line = (); @line = split('\t', $_);

     if ($_ =~ /^#/ ) {
     	@samples = @line[ 1 .. $#line ];
     } else {
     	my $lineage = $line[0];
          push(@lineages, $lineage);

          for my $i ( 1 .. $#line ) {
          	my $sample = $samples[$i - 1];
               $lineage_sample{$lineage}{$sample} = $line[$i];
       	}
	}
}
close IN;
my $nsample = scalar @samples;
my $nlineage = scalar @lineages;

print "the number of samples = $nsample\n";
print "the number of lineages = $nlineage\n";

# read functional profiles for lineages
my %lineage_func = (); my %funcs = ();
for my $lineage ( @lineages ) {

	my $relineage = $lineage;
	$relineage =~ s/:/__/g;
  	$relineage =~ s/;/_/g;	

	my $infile = "$relineage.$funcDB"; my $lineage_size;
     open (IN, "<$lineage_func_dir/$infile") || die "$!";
     while (<IN>) {

		chomp $_;
		$_ = trimwhitespace($_);

		if($_ =~ /^#/){
			(my $lineage, $lineage_size) = split('\t', $_);
			$lineage = substr($lineage, 1, length($lineage) - 1);
		} else {
			(my $func, my $func_fq) = split('\t', $_);
			$lineage_func{$lineage}{$func} = $func_fq/$lineage_size;

			if ( defined $funcs{$func} ) {
				$funcs{$func}++;
			} else {
				$funcs{$func} = 1;
			}
		}
	}
	close IN;
}
my @funcs_ordered = (); @funcs_ordered = sort{$funcs{$b} <=> $funcs{$a}} keys %funcs;
my $nfunc = scalar @funcs_ordered;

print "the number of functions = $nfunc\n";

# make a function sample table by plug-in lineages' functional profiles after normalizing
# the occurrence of functional terms by the number of organisms belonging to lineages
my %func_sample = ();
for my $func ( @funcs_ordered ) {
	for my $sample ( @samples ) { 

		$func_sample{$func}{$sample} = 0.0;

		for my $lineage ( @lineages ) { 

			if ( defined $lineage_func{$lineage}{$func} ) {
				$func_sample{$func}{$sample} = $func_sample{$func}{$sample} + $lineage_sample{$lineage}{$sample} * $lineage_func{$lineage}{$func};
			}
		}
	}
}

# print a function sample table
open (OUT, ">$func_sample_table") || die "$!";
print OUT "#$funcDB"."terms";
for my $sample ( @samples ) {
	print OUT "\t$sample";
}
print OUT "\n";
for my $func ( @funcs_ordered ) {

	print OUT "$func";	

	for my $sample ( @samples ) {
		printf(OUT "\t%.10f", $func_sample{$func}{$sample});
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
