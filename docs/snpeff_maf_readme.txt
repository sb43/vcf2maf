The first 34 columns are standard MAF format, and described here:
https://wiki.nci.nih.gov/x/eJaPAQ

The subsequent 10 columns are relevant to most analyses.

35. HGVSc - the coding sequence of the variant in HGVS recommended format
36. HGVSp - the protein sequence of the variant in HGVS recommended format
37. HGVSp_Short - Same as HGVSp, but using 1-letter amino-acid codes
38. Transcript_ID - transcript onto which the consequence of the variant has been mapped
39. Exon_Number - the exon number (out of total number)
40. t_depth - read depth across this locus in tumor BAM
41. t_ref_count - read depth supporting the reference allele in tumor BAM
42. t_alt_count - read depth supporting the variant allele in tumor BAM
43. n_depth - read depth across this locus in normal BAM
44. n_ref_count - read depth supporting the reference allele in normal BAM
45. n_alt_count - read depth supporting the variant allele in normal BAM

The next column is relevant to analyses that consider the effect of the variant on all alternate
isoforms of the gene, or on non-coding/regulatory transcripts. The effects are sorted first by
transcript biotype priority, then by effect severity, and finally by decreasing order of transcript
length. Each effect in the list is in the format [Gene_Name,Effect,HGVSp,Transcript_ID,].

46. all_effects - a semicolon delimited list of all possible variant effects, sorted by priority

All remaining columns are straight out of the snpEff annotator, as described here:
http://snpeff.sourceforge.net/SnpEff_manual.html#output

47. Effect - Effect of variant. Details here: http://snpeff.sourceforge.net/SnpEff_manual.html#eff
48. Effect_Impact - Effect impact {High, Moderate, Low, Modifier}
49. Functional_Class - Functional class {NONE, SILENT, MISSENSE, NONSENSE}
50. Codon_Change - Old_codon/new_codon OR distance to transcript (for non-coding variants)
51. Amino_Acid_Change - Amino acid change: old_AA/AA_position/new_AA (e.g. 'E30K')
52. Amino_Acid_Length - Length of protein in amino acids (is actually transcription length divided by 3)
53. Gene_Name - Gene name, as per the transcript database used
54. Transcript_BioType - Transcript bioType, if available
55. Gene_Coding - [CODING | NON_CODING]. This field is 'CODING' if any transcript of the gene is marked as protein coding
56. Exon_Rank - Exon rank or Intron rank (e.g. '1' for the first exon, '2' for the second exon, etc.)
57. Genotype_Number - Genotype number corresponding to this effect (e.g. '2' if the effect corresponds to the second ALT)
58. ERRORS - Any errors (not shown if empty)
59. WARNINGS - Any warnings (not shown if empty)
