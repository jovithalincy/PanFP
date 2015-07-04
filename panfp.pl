#!/usr/bin/env perl

# usage:
#
# perl panfp.pl -i otu_sample.txt
#
# written by: Se-Ran Jun

use strict;
use warnings;
use Cwd;

use Getopt::Std;

# setup directory
my $panfp_dir = $dir; # should be fixed if necessary
my $dir = getcwd;

# read arguments ...
my %args;
getopt('i', \%args);
my $otu_sample_table = $args{i}; # input: an otu-sample table

die "usage: panfp.pl -i <otu_sample.txt>\n"
unless defined ($otu_sample_table);

system "$panfp_dir/bin/modify_otu_table.pl -i $dir/$otu_sample_table -d $panfp_dir/data/org_chr_16SrRNA_KO.TH30.txt -o $dir/otu_table_modified.txt";

system "$panfp_dir/bin/lineage_16Scopynum.pl -i $dir/otu_table_modified.txt -d $panfp_dir/data/org_chr_16SrRNA_KO.TH30.txt -o $dir/lineage_16Scopynum.txt";

system "$panfp_dir/bin/lineage_func.pl -i $dir/otu_table_modified.txt -d $panfp_dir/data/org_chr_16SrRNA_KO.TH30.txt -a $panfp_dir/annotation/KO -o $dir/lineage_func";

system "$panfp_dir/bin/otu_table_norm_by_copynum.pl -i $dir/otu_table_modified.txt -c $dir/lineage_16Scopynum.txt -o $dir/otu_table_modified_norm_by_copynum.txt";

system "$panfp_dir/bin/otu_table_norm_by_samplesize.pl -i $dir/otu_table_modified_norm_by_copynum.txt -o $dir/otu_table_modified_norm_by_copynum_samplesize.txt";

system "$panfp_dir/bin/lineage_sample.pl -i $dir/otu_table_modified_norm_by_copynum_samplesize.txt -o $dir/lineage_sample.txt";

system "$panfp_dir/bin/func_sample.pl -i $dir/lineage_sample.txt -l $dir/lineage_func -o $dir/func_sample.txt";

system "rm -f $dir/otu_table_modified.txt";
system "rm -f $dir/lineage_16Scopynum.txt";
system "rm -f $dir/otu_table_modified_norm_by_copynum.txt";
system "rm -f $dir/otu_table_modified_norm_by_copynum_samplesize.txt";
system "rm -f $dir/lineage_sample.txt";
system "rm -rf $dir/lineage_func";
