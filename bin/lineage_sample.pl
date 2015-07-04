#!/usr/bin/env perl

# usage:
#
# perl bin/lineage_sample.pl -i <full_path_to/otu_table_modified_norm_by_copynum_samplesize.txt> -o <full_path_to/lineage_sample.txt>
#
# written by: Se-Ran Jun

use warnings;
use strict;

use Getopt::Std;

my %args;

# read arguments ...
getopt('i:o', \%args);

my $otu_table = $args{i}; # input: an otu table normalized by 16S copy number and sample size
my $lineage_sample_table = $args{o}; # output: a lineage smaple table 

die "usage: bin/lineage_sample.pl -i <full_path_to/otu_table_modified_norm_by_copynum_samplesize.txt> -o <full_path_to/lineage_sample.txt>\n"
unless defined ($otu_table) && defined ($lineage_sample_table);

# read an otu-table
my @samples = (); my @otus = (); my %otu_sample = (); my %otu_lineage = (); my %lineages = ();
open ( IN, "<$otu_table" ) || die "$!";
while (<IN>) {

	chomp $_;
	$_ = trimwhitespace($_);
	my @line = (); @line = split('\t', $_);

     if ( $_ =~ /^#/ ) {
     	@samples = @line[1 ..  ($#line - 1) ];
     } else {
		my $otu = $line[0];
		my $lineage = $line[-1];

		push(@otus, $otu);
		$otu_lineage{$otu} = $lineage;

		if (defined $lineages{$lineage} ) {
			$lineages{$lineage}++;
		} else {
			$lineages{$lineage} = 1;
		}

		for my $i ( 1 .. ($#line  - 1) ) {
			my $sample = $samples[$i - 1];
			$otu_sample{$otu}{$sample} = $line[$i];
		}
	}
}
close IN;
my $notu = scalar @otus;
my $nsample = scalar @samples;
my $nlineage = scalar keys %lineages;

print "the number of otus = $notu\n";
print "the number of samples = $nsample\n";
print "the number of lineages = $nlineage\n";

# for printing lineages in order by the taxonomic levels that this step is not necessary
my %kp = (); my %pc = (); my %co = (); my %of = (); my %fg = (); my %gs = ();
foreach my $lineage ( keys %lineages ) {

	my @level = (); @level = split(';', $lineage);

	for my $i ( 0 .. ($#level - 1) ) {
		if (($level[$i] =~ /k:/) && ($level[$i+1] =~ /p:/)) {
          	$kp{$level[$i]}{$level[$i+1]} = 1;
		} elsif (($level[$i] =~ /p:/) && ($level[$i+1] =~ /c:/)) {
          	$pc{$level[$i]}{$level[$i+1]} = 1;
		} elsif (($level[$i] =~ /c:/) && ($level[$i+1] =~ /o:/)) {
          	$co{$level[$i]}{$level[$i+1]} = 1;
		} elsif (($level[$i] =~ /o:/) && ($level[$i+1] =~ /f:/)) {
          	$of{$level[$i]}{$level[$i+1]} = 1;
		} elsif (($level[$i] =~ /f:/) && ($level[$i+1] =~ /g:/)) {
         		$fg{$level[$i]}{$level[$i+1]} = 1;
		} elsif (($level[$i] =~ /g:/) && ($level[$i+1] =~ /s:/)) {
          	$gs{$level[$i]}{$level[$i+1]} = 1;
          }
	}
}

# thip step is also not necessary
my @lineages_ordered = ();
foreach my $kingdom (keys %kp) {
	foreach my $phylum (keys %{$kp{$kingdom}}) {

     	my $lineage = "$kingdom;$phylum";

          if (defined $lineages{$lineage}) {
			push (@lineages_ordered, $lineage);
          }

          foreach my $class (keys %{$pc{$phylum}}){

          	my $lineage = "$kingdom;$phylum;$class";

               if (defined $lineages{$lineage}) {
				push(@lineages_ordered, $lineage);
            	}

              	foreach my $order (keys %{$co{$class}}){

               	my $lineage = "$kingdom;$phylum;$class;$order";

                 	if (defined $lineages{$lineage}) {
					push(@lineages_ordered, $lineage);
                	}

                   	foreach my $family (keys %{$of{$order}}) {

                    	my $lineage = "$kingdom;$phylum;$class;$order;$family";

                         if (defined $lineages{$lineage}) {
						push(@lineages_ordered, $lineage);
                       	}

                       	foreach my $genus (keys %{$fg{$family}}) {

                         	my $lineage = "$kingdom;$phylum;$class;$order;$family;$genus";

                              if (defined $lineages{$lineage}) {
							push(@lineages_ordered, $lineage);
                             	}

						foreach my $species ( keys %{$gs{$genus}} ) {

							my $lineage = "$kingdom;$phylum;$class;$order;$family;$genus;$species";

                                  	if (defined $lineages{$lineage}) {
                                  		push(@lineages_ordered, $lineage);
                                  	}
						}
                       	}
              		}
           	}
    		}
	}   
}
print "the number of lineages ordered = ",scalar @lineages_ordered,"\n";

# convert an otu-table into a lineage-sample table
my %lineage_sample = ();
for my $lineage ( @lineages_ordered ) { # lineages in order by the taxonomic levels
	for my $sample ( @samples ) {

		$lineage_sample{$lineage}{$sample} = 0.0;

		for my $otu ( @otus ) {
			if($otu_lineage{$otu} eq $lineage){
				$lineage_sample{$lineage}{$sample} = $lineage_sample{$lineage}{$sample} + $otu_sample{$otu}{$sample};
			}
		}
	}
}

# print a lineage-sample table
open (OUT, ">$lineage_sample_table" ) || die "$!";
print OUT "#Lineage";
for my $sample ( @samples ) {
	print OUT "\t$sample";
}
print OUT "\n";
for my $lineage ( @lineages_ordered ) {

	print OUT "$lineage";

	for my $sample ( @samples ) {
     	printf(OUT "\t%.10f", $lineage_sample{$lineage}{$sample});
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
