## PanFP: Pangenome-based functional profiles for microbial communities

The software called PanFP (Pangenome-based functional profile) is a tool to predict functional 
composition of microbial communities based on communities survey data which builds a pangenome for 
each OTU and assigns the pangenome's functional profile to the OTU. For microbial communities for 
which metagenomes are not available, the PanFP would provide new ways to study complex ecosystems
comparative functional metatgnomes and metadata analysis statistically highlighting the most
overpresented biological annotation and which functions are critical between environmental samples.

The PanFP contains the following three directories:

###1. data
The directory includes a file, **org_chr_16SrRNA_KO.TH30.txt**, which is a list of orgranisms whose 
(1) chromosomal complete proteomes are available from **NCBI**, (2) 16S rRNA seqeunces are identified,
and (3) functional coverages by KEGG Orthology (KO) are at least 30%.

The file contains the following information in a tab delimited format: 

  * _the name of organism_,
  * _the number of chromosomes_,
  * _the number of 16S copy number_,
  * _the number of proteins_,
  * _the number of proteins with functional annotation by Pfam_,
  * _the full lineage obtained from **NCBI**_.

###2. annotation
The direcotory contains functional annotations by KEGG Orthology for organisms in a file, 
**org_chr_16SrRNA_KO.TH30.txt**. 

For example, a file named **Yersinia_pestis_biovar_Medievalis_Harbin_35_uid158537.GO** contains proteins
with **NCBI** reference ids along with assigned KO terms, which are ordered by locus on a chromosome(s) for
an organism named **Yersinia_pestis_biovar_Medievalis_Harbin_35_uid158537** by **NCBI** in a tab delimited 
format. In the file, the first line contains the organism's full lineage, the number of proteins, the number of 
proteins with functional annotation by KEGG Orthology. The functional annotation was done based on cross 
reference between **NCBI** and **Uniprot**.

###3. bin
The directory contains seven Perl scripts: These scripts have been tested on Perl 5.16 and use a collection
of basic statistics moduels **Statistics::Basic**.

  * modify_otu_table.pl
  * lineage_16Scopynum.pl
  * lineage_func.pl
  * otu_table_norm_by_copynum.pl
  * otu_table_norm_by_samplesize.pl
  * lineage_sample.pl
  * func_sample.pl

##Usage

The PanFP takes an OTU table of microbial communities as an input file, and produces the predicted
functional profiles of the microbial communities by performing the following seven steps:

####STEP 1:####
```
bin/modify_otu_table.pl -i <full_path_to/otu_table.txt> -d <full_path_to/org_chr_16SrRNA_KO.TH30.txt> -o <full_path_to/modified_otu_table.txt>
```

The script (1) modifies OTU's lineage based on **NCBI** taxonomy nomenclature, (2) modifies OTU's lineage 
that the OTU can have at least one organism whose taxonomy belongs to its lineage, (3) discards OTUs with 
no lineage information, and (4) discards OTUs if they have only kingdom information.

input:
  * an otu table
  * a list of organisms with functional annotations stored in a directory, **PanFP/data**

output:
  * an otu table with modified lineage

An output, otu table should be in a tab delimited format as follows:

|\#OTU|sample1|...|sampleN|ConsensusLineage|
|:----|:------|:--|:------|:---------------|
|4479946|0.0|...|2.0|k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__MND1; f__; g__; s__|
|64356|4.0|...|24.0|k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__; f__; g__; s__|
|...|...|...|...|...|
|247460|0.0|...|2.0|k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacteriales; f__Enterobacteriaceae; g__Gluconacetobacter; s__liquefaciens|

where the first column represents OTU ids, numbers represent raw frequency of 16S rRNA, and the last column
represents lineage of OTUs.

An output modified otu-table is in a tab delimited format as follows:

|\#OTU|sample1|...|sampleN|modifiedConsensusLineage|
|:----|:------|:--|:------|:-----------------------|
|4479946|0.0|...|2.0|k:Bacteria;p:Proteobacteria;c:Betaproteobacteria|
|64356|4.0|...|24.0|k:Bacteria;p:Proteobacteria;c:Betaproteobacteria|
|...|...|...|...|...|
|247460|0.0|...|2.0|k:Bacteria;p:Proteobacteria;c:Gammaproteobacteria;o:Enterobacteriales;f:Enterobacteriaceae;g:Gluconacetobacter;s:liquefaciens|

where the first column represents OTU ids, numbers represent raw frequency of 16S rRNA, and the last column
represents lineage of OTUs.

####STEP 2:####

```
bin/lineage_16Scopynum.pl -i <full_path_to/otu_table_modified.txt> -d <full_path_to/org_chr_16SrRNA_KO.TH30.txt> -o <full_path_to/lineage_16Scopynum.txt>
```

The script produces distributions of 16S copy numbers for lineages. For a given taxon, the script collects
organisms belonging to the taxon, extracts their 16S copy numbers, and calcualtes an average and the median
16S copy number for each OTU through its lineage.

input:
  * an otu table with modified lineage
  * a list of organisms for a functional database of interest stored in a directory, **PanFP/data**

ouput:
  * distributions of 16S copy numbers for lineages

An output 16Scopynumber distribution is in a tab delimited format as follows:

|#lineage|median|avg|16Scopynums|
|:-------|:-----|:--|:----------|
|k:Bacteria;p:Chlorobi|2|1.727|2  1   1   2   2   1   2   1   2   3   2|
|k:Bacteria;p:Firmicutes;c:Bacilli;o:Bacillales;f:Alicyclobacillaceae;g:Alicyclobacillus|6|6.000|6   6|
|...|...|...|...|
|k:Bacteria;p:Proteobacteria;c:Betaproteobacteria;o:Burkholderiales;f:undef;g:Leptothrix|2|2.000|2|

where the first column represents lineage present in a modified otu-table, the second column represents
the median 16S copy number, third column represents an average 16S copy number, and the rest of columns 
represent 16S copy numbers of organisms belong to a given lineage.

####STEP 3:####

```
bin/lineage_func.pl -i <full_path_to/otu_table_modified.txt> -d <full_path_to/org_chr_16SrRNA_KO.TH30.txt> -a <annotation/KO> -o <full_path_to/lineage_func>
```

The script collects organisms belonging to the taxon, and then calculates the occurrence of functional 
terms in the set of organisms.

input:
  * an otu table with modified lineage
  * a list of organisms with functional annotations stored in a directory, **PanFP/data**
  * functional annotations of organisms stored in a directory, **PanFP/annotation/KO**

output:
  * functional profiles of lineages present in a modified otu-table

For a lineage **k:Bacteria;p:Verrucomicrobia**, an output file named **k__Bacteria_p__Verrucomicrobia.KO**
is formatted as follows:

|#k:Bacteria;p:Proteobacteria;c:Gammaproteobacteria;o:Chromatiales|12   |
|:-----------------------------|:---|
|K01698|12|
|K05602|1|
|......|...|
|K01154|27|
|K02026|2|

where the first line represents a taxon present in a modified otu-table and the number of ogranisms belonging 
to the taxon. From the second line, the first column represents functional term occurred at least one time
in a set of organisms and the second column represents their frequency among 4 organisms for this case.

####STEP 4:####

```
bin/otu_table_norm_by_copynum.pl -i <full_path_to/otu_table_modified.txt> -c <full_path_to/lineage_16Scopynum.txt> -o <full_path_to/otu_table_modified_norm_by_copynum.txt>
```

A script normalizes an otu-table by the median 16S copy number of OTUs' lineages.

input:
  * an otu table with modified lineage
  * a list of lineages present in a modified otu-table with their 16S copy number distributions

output:
  * an otu table normalized by the median 16S copy number of OTUs' lineages

An output, otu table normalized by the median 16S copy number of OTUs' lineage is in a tab delimited
format as follows:

|\#OTU|sample1|...|sampleN|modifiedConsensusLineage|
|:----|:------|:--|:------|:-----------------------|
|4479946|0.0|...|0.33|k:Bacteria;p:Proteobacteria;c:Betaproteobacteria|
|64356|2.0|...|12.0|k:Bacteria;p:Proteobacteria;c:Betaproteobacteria|
|.....|...|...|....|...|
|247460|0.0|...|0.5|k:Bacteria;p:Proteobacteria;c:Gammaproteobacteria;o:Enterobacteriales;f:Enterobacteriaceae;g:Gluconacetobacter;s:liquefaciens|

where the first column represents OTU ids, the numbers are normalized 16S abundance which is more likely
an approximation of organismal abundances, and the last column represents OTU's modified lineage.

####STEP 5:####

```
bin/otu_table_norm_by_samplesize.pl -i <full_path_to/otu_table_modified_norm_by_copynum.txt> -o <full_path_to/otu_table_modified_norm_by_copynum_samplesize.txt>
```

The script normalizes an otu-table by the sample size which is defined as the sum of organismal abundances
that the final abundances represent approximations of organismal abundances which are comparable across
different samples.

input:
  * an otu table normalized by the median 16S copy number of OTUs' lineages

output:
  * an otu table normalized by the sample size

An output, otu table normalized by the sample size is in a tab delimited format as follows:

|\#OTU|sample1|...|sampleN|modifiedConsensusLineage|
|:----|:------|:--|:------|:-----------------------|
|4479946|0.0|...|0.00345|k:Bacteria;p:Proteobacteria;c:Betaproteobacteria|
|64356|0.00221|...|0.07892|k:Bacteria;p:Proteobacteria;c:Betaproteobacteria|
|...|...|...|...|...|
|247460|0.0|...|0.01234|k:Bacteria;p:Proteobacteria;c:Gammaproteobacteria;o:Enterobacteriales;f:Enterobacteriaceae;g:Gluconacetobacter;s:liquefaciens|

where the first column represents OTU id, the numbers represent organismal abundances comparable across different
smaples, and the last column represents OTU's modified lineage.

####STEP 6:####

```
bin/lineage_sample.pl -i <full_path_to/otu_table_modified_norm_by_copynum_samplesize.txt> -o <full_path_to/lineage_sample.txt>
```

The script converts a otu-table normalized by the median 16S copy number and the sample size into
a lineage sample table by summing abundances of organisms belong to a given lineage.

input:
  * an otu table normalized by the median 16S copy number and the sample size 

output:
  * a lineage sample table

An output lineage-sample table is in a tab delimited format as follows:

|\#Lineage|sample1|...|sampleN|
|:--------|:------|:--|:------|
|k:Archaea;p:Nanoarchaeota|0.0002448750|...|0.0000000000|
|k:Bacteria;p:Synergistetes;c:Synergistia;o:Synergistales;f:Synergistaceae;g:Aminobacter|0.0000544167|...|0.0000411905|
|...|...|...|...|
|k:Bacteria;p:Actinobacteria;c:Actinobacteria|0.0000272083|...|0.0002780357|

####STEP 7:####

```
bin/func_sample.pl -i <full_path_to/lineage_sample.txt> -l <full_path_to_lineage_func> -o <full_path_to/func_sample.txt>
```

The script converts an lineage-sample table into a function-sample table by plug-in lineages' functional
profiles after normalizing the occurrence of functional terms by the number of organisms belonging to
lineages that the frequency less likely depends on the lineage size (the number of organisms pooled within
the taxon which varies greatly across different taxons)

input:
  * a lineage sample table
  * a directory which includes lineages' functional profiles

output:
  * a function sample table

An output function sample table is in a tab delimited format as follows:

|#KOs   |sample1|...|sampleN|
|:------|:------|:--|:------|
|K02946|0.9951432689|...|1.0065902626|
|K01265|1.6652932004|...|1.3552325900|
|...|...|...|...|
|K02867|1.0209092188|...|1.0006068186|

where the first line represents sample ids and the rest of lines represents functional terms with their frequency 
in samples.

### For a quick start:
A script, panfp.pl takes an otu-table as an input, calls Perl scripts described above in order based on
the algorithm, and produces a final result, a function-sample table.

```
panfp.pl -i otu_table.txt
```

#### Author

Direct questions and comments to the author, Se-Ran Jun, at https://github.com/srjun.
