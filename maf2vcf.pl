#!/usr/bin/env perl

# maf2vcf - Reformat variants in a given MAF into generic VCFs with GT:AD:DP data if available

use strict;
use warnings;
use IO::File;
use Getopt::Long qw( GetOptions );
use Pod::Usage qw( pod2usage );

# Set any default paths and constants
my $ref_fasta = "$ENV{HOME}/.vep/homo_sapiens/78_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa";
my ( $tum_depth_col, $tum_rad_col, $tum_vad_col ) = qw( t_depth t_ref_count t_alt_count );
my ( $nrm_depth_col, $nrm_rad_col, $nrm_vad_col ) = qw( n_depth n_ref_count n_alt_count );

# Find out where samtools is installed, and warn the user if it's not
my $samtools = ( -e "/opt/bin/samtools" ? "/opt/bin/samtools" : "/usr/bin/samtools" );
$samtools = `which samtools` unless( -e $samtools );
chomp( $samtools );
( $samtools and -e $samtools ) or die "ERROR: Please install samtools, and make sure it's in your PATH\n";

# Check for missing or crappy arguments
unless( @ARGV and $ARGV[0]=~m/^-/ ) {
    pod2usage( -verbose => 0, -message => "$0: Missing or invalid arguments!\n", -exitval => 2 );
}

# Parse options and print usage syntax on a syntax error, or if help was explicitly requested
my ( $man, $help ) = ( 0, 0 );
my ( $input_maf, $output_dir );
GetOptions(
    'help!' => \$help,
    'man!' => \$man,
    'input-maf=s' => \$input_maf,
    'output-dir=s' => \$output_dir,
    'ref-fasta=s' => \$ref_fasta,
    'tum-depth-col=s' => \$tum_depth_col,
    'tum-rad-col=s' => \$tum_rad_col,
    'tum-vad-col=s' => \$tum_vad_col,
    'nrm-depth-col=s' => \$nrm_depth_col,
    'nrm-rad-col=s' => \$nrm_rad_col,
    'nrm-vad-col=s' => \$nrm_vad_col
) or pod2usage( -verbose => 1, -input => \*DATA, -exitval => 2 );
pod2usage( -verbose => 1, -input => \*DATA, -exitval => 0 ) if( $help );
pod2usage( -verbose => 2, -input => \*DATA, -exitval => 0 ) if( $man );

# Check that the minimum inputs needed 
( defined $input_maf and defined $output_dir ) or die "ERROR: At least input-maf and output-dir must be defined!\n";

# Parse through each variant in the MAF, and fill up the respective VCFs
my $maf_fh = IO::File->new( $input_maf ) or die "ERROR: Couldn't open file: $input_maf\n";
my $line_count = 0;
my %col_idx = (); # Hash to map column names to column indexes
while( my $line = $maf_fh->getline ) {

    # Skip comment lines
    next if( $line =~ m/^#/ );

    # Instead of a chomp, do a thorough removal of carriage returns, line feeds, and prefixed/suffixed whitespace
    my @cols = map{s/^\s+|\s+$|\r|\n//g; $_} split( /\t/, $line );

    # Parse the header line to map column names to their indexes
    if( $line =~ m/^(Hugo_Symbol|Chromosome)/ ) {
        my $idx = 0;

        # Fetch the column names and do some sanity checks (don't be case-sensitive)
        map{ my $c = lc; $col_idx{$c} = $idx; ++$idx; } @cols;
        map{ my $c = lc; ( defined $col_idx{$c} ) or die "ERROR: $_ is a required MAF column!\n" } qw( Chromosome Start_Position Reference_Allele Tumor_Sample_Barcode );
        ( defined $col_idx{tumor_seq_allele1} or defined $col_idx{tumor_seq_allele2} ) or die "ERROR: At least one MAF column for Tumor_Seq_Allele must be defined!\n";

        # Fetch all tumor-normal paired IDs from the MAF, doing some whitespace cleanup in the same step
        my $tn_idx = $col_idx{tumor_sample_barcode} + 1;
        $tn_idx .= ( "," . ( $col_idx{matched_norm_sample_barcode} + 1 )) if( defined $col_idx{matched_norm_sample_barcode} );
        my @tn_pair = map{s/^\s+|\s+$|\r|\n//g; s/\s*\t\s*/\t/; $_}`egrep -v "^#|^Hugo_Symbol|^Chromosome" $input_maf | cut -f $tn_idx | sort -u`;

        # For each TN-pair in the MAF, initialize blank VCFs with proper VCF headers in output directory
        unless( -e $output_dir ) { mkdir $output_dir or die "ERROR: Couldn't create directory $output_dir! $!"; }
        foreach my $pair ( @tn_pair ) {
            my ( $t_id, $n_id ) = split( /\t/, $pair );
            $n_id = "NORMAL" unless( defined $n_id ); # Use a placeholder name for normal if its undefined
            my $vcf_file = "$output_dir/$t_id\_vs_$n_id.vcf";
            my $vcf_fh = IO::File->new( $vcf_file, ">" );
            $vcf_fh->print( "##fileformat=VCFv4.2\n" );
            $vcf_fh->print( "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" );
            $vcf_fh->print( "##FORMAT=<ID=AD,Number=G,Type=Integer,Description=\"Allelic Depths of REF and ALT(s) in the order listed\">\n" );
            $vcf_fh->print( "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n" );
            $vcf_fh->print( "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$t_id\t$n_id\n" );
            $vcf_fh->close;
        }
        next;
    }

    # Print an error if we got to this point without parsing a header line, and increment a counter for all non-header lines
    ( %col_idx ) or die "ERROR: Couldn't find a header line in the MAF: $input_maf";
    $line_count++;

    # For a variant in the MAF, parse out the bare minimum data needed by a VCF
    my ( $chr, $pos, $ref, $al1, $al2, $t_id, $n_id, $n_al1, $n_al2 ) = map{ my $c = lc; ( defined $col_idx{$c} ? $cols[$col_idx{$c}] : undef )} qw( Chromosome Start_Position Reference_Allele Tumor_Seq_Allele1 Tumor_Seq_Allele2 Tumor_Sample_Barcode Matched_Norm_Sample_Barcode Match_Norm_Seq_Allele1 Match_Norm_Seq_Allele2 );

    # Parse out read counts for ref/var alleles, if available
    my ( $t_dp, $t_rad, $t_vad, $n_dp, $n_rad, $n_vad ) = map{my $c = lc; (( defined $col_idx{$c} and defined $cols[$col_idx{$c}] and $cols[$col_idx{$c}] =~ m/^\d+/ ) ? sprintf( "%.0f", $cols[$col_idx{$c}] ) : '.' )} ( $tum_depth_col, $tum_rad_col, $tum_vad_col, $nrm_depth_col, $nrm_rad_col, $nrm_vad_col );

    # Normal sample ID could be undefined for legit reasons, but we need a placeholder name
    $n_id = "NORMAL" unless( $n_id );

    # If normal alleles are unset in the MAF (quite common), assume homozygous reference
    $n_al1 = $ref unless( $n_al1 );
    $n_al2 = $ref unless( $n_al2 );

    # Make sure we have at least 1 variant allele. If 1 is unset, set it to the reference allele
    if( !$al1 and !$al2 ) { 
        warn "WARNING: Skipping variant at $chr:$pos without any variant alleles specified!\n";
        next;
    }
    $al1 = $ref unless( $al1 );
    $al2 = $ref unless( $al2 );

    # To represent indels in VCF format, we need to fetch the preceding bp from a reference FASTA
    my ( $ref_len, $al1_len, $al2_len ) = map{( $_=~m/^(\?|-|0)$/ ? 0 : length( $_ )) } ( $ref, $al1, $al2 );
    if( $ref_len == 0 or $al1_len == 0 or $al2_len == 0 ) {
        --$pos if( $ref_len > $al1_len or $ref_len > $al2_len ); # Decrement POS for deletions only
        my $prefix_bp = `$samtools faidx $ref_fasta $chr:$pos-$pos | grep -v ^\\>`;
        chomp( $prefix_bp );
        $prefix_bp = uc( $prefix_bp );
        ( $prefix_bp =~ m/^[ACGTN]$/ ) or die "ERROR: Cannot retreive bp at $chr:$pos! Please specify --ref-fasta appropriately\n";
        # Blank out the dashes (or other weird chars) used with indels, and prefix the fetched bp
        ( $ref, $al1, $al2, $n_al1, $n_al2 ) = map{s/^(\?|-|0)$//; $_=$prefix_bp.$_} ( $ref, $al1, $al2, $n_al1, $n_al2 );
    }

    # To simplify setting tumor genotype later, ensure that $al2 is always non-REF
    ( $al1, $al2 ) = ( $al2, $al1 ) if( $al2 eq $ref );
    # Do the same for the normal alleles, though it makes no difference if both are REF
    ( $n_al1, $n_al2 ) = ( $n_al2, $n_al1 ) if( $n_al2 eq $ref );

    # Fill an array with all unique REF/ALT alleles, and set their 0-based indexes like in a VCF
    # Notice how we ensure that $alleles[0] is REF and #alleles[1] is the major ALT allele in tumor
    my ( @alleles, %al_idx );
    my $idx = 0;
    foreach my $al ( $ref, $al2, $al1, $n_al2, $n_al1 ) {
        unless( defined $al_idx{$al} ) {
            push( @alleles, $al );
            $al_idx{$al} = $idx++;
        }
    }

    # Set tumor and normal genotypes (FORMAT tag GT in VCF)
    my $t_gt = join( "/", $al_idx{$al2}, $al_idx{$al1} );
    my $n_gt = join( "/", $al_idx{$n_al2}, $al_idx{$n_al1} );

    # Create the VCF's comma-delimited ALT field that must list all non-REF (variant) alleles
    my $alt = join( ",", @alleles[1..$#alleles] );

    # If there are >1 variant alleles, assume that depths in $t_vad and $n_vad are for $al2
    if( scalar( @alleles ) > 2 ) {
        $t_vad = join( ",", $t_vad, map{"."}@alleles[2..$#alleles] );
        $n_vad = join( ",", $n_vad, map{"."}@alleles[2..$#alleles] );
    }

    # Construct genotype fields for FORMAT tags GT:AD:DP
    my $t_fmt = "$t_gt:$t_rad,$t_vad:$t_dp";
    my $n_fmt = "$n_gt:$n_rad,$n_vad:$n_dp";

    # Contruct a VCF formatted line and append it to the respective VCF
    my $vcf_file = "$output_dir/$t_id\_vs_$n_id.vcf";
    my $vcf_line = join( "\t", $chr, $pos, ".", $ref, $alt, qw( . . . ), "GT:AD:DP", $t_fmt, $n_fmt );
    my $vcf_fh = IO::File->new( $vcf_file, ">>" );
    $vcf_fh->print( "$vcf_line\n" );
    $vcf_fh->close;
}
$maf_fh->close;

# Make sure that we handled a positive non-zero number of lines in the MAF
( $line_count > 0 ) or die "ERROR: No variant lines in the input MAF!\n";

__DATA__

=head1 NAME

 maf2vcf.pl - Reformat variants in a given MAF into generic VCFs with GT:AD:DP data if available

=head1 SYNOPSIS

 perl maf2vcf.pl --help
 perl maf2vcf.pl --input-maf test.maf --output-dir vcfs

=head1 OPTIONS

 --input-maf      Path to input file in MAF format
 --output-dir     Path to output directory where VCFs will be stored, one per TN-pair
 --ref-fasta      Path to reference Fasta file [~/.vep/homo_sapiens/78_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa]
 --tum-depth-col  Name of MAF column for read depth in tumor BAM [t_depth]
 --tum-rad-col    Name of MAF column for reference allele depth in tumor BAM [t_ref_count]
 --tum-vad-col    Name of MAF column for variant allele depth in tumor BAM [t_alt_count]
 --nrm-depth-col  Name of MAF column for read depth in normal BAM [n_depth]
 --nrm-rad-col    Name of MAF column for reference allele depth in normal BAM [n_ref_count]
 --nrm-vad-col    Name of MAF column for variant allele depth in normal BAM [n_alt_count]
 --help           Print a brief help message and quit
 --man            Print the detailed manual

=head1 DESCRIPTION

This script breaks down variants in a MAF into VCFs for each tumor-normal pair, in preparation for annotation with vcf2maf

=head2 Relevant links:

 Homepage: https://github.com/ckandoth/vcf2maf
 VCF format: http://samtools.github.io/hts-specs/
 MAF format: https://wiki.nci.nih.gov/x/eJaPAQ

=head1 AUTHORS

 Cyriac Kandoth (ckandoth@gmail.com)

=head1 LICENSE

 Apache-2.0 | Apache License, Version 2.0 | https://www.apache.org/licenses/LICENSE-2.0

=cut
