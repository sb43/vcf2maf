=======================================================================
BCCA GENOME SCIENCES CENTRE TCGA EXOME DATA FILE HEADER INFORMATION
=======================================================================
DCC Level 2 data:
-----------------

The Strelka VCF format is explained here:
https://sites.google.com/site/strelkasomaticvariantcaller/home/somatic-variant-output

1) .snv.vcf
A tab-delimited text file containing putative single nucleotide variants identified from exome data in the VCF format.

2) .indel.vcf
A tab-delimited text file containing putative indels identified from exome data in the VCF format.


=======================================================================
BCCA GENOME SCIENCES CENTRE TCGA EXOME ANALYIS ALGORITHM DESCRIPTION
=======================================================================

Strelka [1] version 1.0.6 was used to identify somatic single nucleotide variants and short insertions and deletions from the TCGA MESO exome dataset. All parameters were set to default with the exception of "isSkipDepthFilters" which was set to 1 in order to skip depth filtration. 83 pairs of libraries were analyzed. If blood sample was available, it served as the matched normal. Otherwise, matched normal tissue was used.


[1] Saunders CT., Wong WS., Swamy S., Becq J., Murray LJ., Cheetham RK (2012) Strelka: accurate somatic small-variant calling from sequenced tumor-normal sample pairs. Bioinformatics, 28, 1811-7. [PMID: 22581179]