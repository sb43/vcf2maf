#!/usr/bin/env perl

# vcf2maf - Convert a VCF into a MAF by mapping each variant to only one of all possible gene isoforms

use strict;
use Data::Dumper;
use warnings;
use Capture::Tiny qw(:all);
use Try::Tiny qw(try catch finally);
use Pod::Usage qw(pod2usage);
use Getopt::Long;
use File::Path qw(mkpath);

use Vcf;







my $file_path=$ARGV[0];
my $file_extension = $ARGV[1];


try{

my ($options) = option_builder();

my $file_path=$options->{'d'};
my $file_extension = $options->{'e'};
my $disease = $options->{'p'};
my $outdir = $options->{'o'};

if (! -e $outdir) {
	mkpath($outdir);
}

my $type='somatic';

open(my $cgp_fh, '<', 'project_data/all_samples_sent.txt') or die "Could not open ref file\n";

my $sample_names=undef;

	while(<$cgp_fh>) {
			chomp;
			my ($ref_UUID, $ref_PD_ID) = split /\t/, $_;
			$sample_names->{uc($ref_UUID)} = $ref_PD_ID
	}
	close $cgp_fh;

	open(my $cghub_fh, '<', 'project_data/TCGA_MESO_Data_Summary.txt') or die "Could not open ref2 file";

	while(<$cghub_fh>) {
			chomp;
			my($REF_barcode,$REF_Ali) = (split /\t/, $_)[1,17];
			if ( exists $sample_names->{uc($REF_Ali)} ) {
					$sample_names->{$sample_names->{uc($REF_Ali)}}=[uc($REF_Ali),$REF_barcode];															  
			}
	}
	close $cghub_fh;
	my $ext="pindel_somatic";
	foreach my $sample (keys %$sample_names) {
		if( ref($sample_names->{$sample}) eq 'ARRAY' ) {
			my $input_vcf=$file_path.$sample.$file_extension;
			#my($dir,$file,$ext)=File::parse();
			if( -s $input_vcf.'.gz' ) {
					my $vcf = Vcf->new(file => $input_vcf.'.gz');
					$vcf->parse_header();	
					$vcf->recalc_ac_an(0);
					my($normal_id)=get_normal_id($vcf);
					#print "$sample---$normal_id------$sample_names->{$sample}[0]-----$sample_names->{$sample}[1]------$sample_names->{$normal_id}[0]-----$sample_names->{$normal_id}[1]\n";
					my $cmd = "bgzip -f -c -d $input_vcf.gz >$input_vcf";
					#run_cmd($cmd);
					#$input_vcf=~s/\.gz//g;
					$cmd = "bsub -o $sample\_$ext.o -e $sample\_$ext.e \'perl vcf2maf.pl --input-vagrent $input_vcf ".
								" --tumor-id $sample_names->{$sample}[1] --normal-id $sample_names->{$normal_id}[1] ".
								" --vcf-tumor-id TUMOUR --vcf-normal-id NORMAL ".
								" --tumor-uuid $sample_names->{$sample}[0]  --normal-uuid $sample_names->{$normal_id}[0] ".
								" --maf-center WTSI --seq-source WXS ".
								" --output-maf $outdir/$sample\_$ext.maf ".
								" --dbsnp /nfs/users/nfs_s/sb43/scripts/vcf2maf/project_data/dbsnp_subset/0000.vcf.gz ".
								" --process-flag PASS ".
								" --somatic-maf 1 \' ";
								#" --dbsnp /nfs/users/nfs_s/sb43/scratch_tmp_storage_not_backed_up/vcf2MAF/00-All.vcf.gz \n";
								#" --process-flag PASS";
				   
				      
					     run_cmd($cmd);
				  #print "$cmd\n";
					#exit;
				#print "$sample\n";
			}		
		}
	}
	
	print "completed analysis -- time to merge";
	#awk '($1!~/Hugo_Symbol/ && $1!~/#version/)' *.maf >../pindel.maf
	#sanger.ac.uk_".$disease.'IlluminaHiSeq_DNASeq.Level_2.1.somatic.maf'
	

}
catch{
	print "Error occured $_\n";
};


sub get_normal_id {
	my ($vcf)=@_;
	my $sample_info=$vcf->get_header_line(key=>'SAMPLE', ID=>'NORMAL' );
	foreach my $sample_line ( @$sample_info) {
		foreach my $key (keys %$sample_line) {
			if(defined $sample_line->{$key} and $key eq "SampleName") {
				return $sample_line->{$key};
			}
		}
	}
}



=head2 _run_cmd
runs external command
Inputs
=over 2
=item cmd - command to run
=back
=cut

sub run_cmd {
	my($cmd)=@_;
	my ($out,$stderr,$exit)=capture{system($cmd)};
	if($exit) {
			print "Failed to run <<<<<<< \n $cmd  <<<<<< \n with status <<<<<< \n OUT:\n $out  :ERR:\n $stderr EXIT:\n $exit \n <<<<<<< \n";
	}
	else {
		print "\ncommand <<<<<< \n $cmd \nrun successfully <<<<<<<<< ";
	}
	return $out;
}


sub option_builder {
	my ($factory) = @_;
	my %opts;
	&GetOptions (
					'h|help'    => \$opts{'h'},
					'd|vcfDir=s' => \$opts{'d'},
					'e|vcfExt=s' => \$opts{'e'},
					'p|project=s' => \$opts{'p'},
					'o|outdir=s'  => \$opts{'o'},
					'v|version'  => \$opts{'v'},
	);

  pod2usage(-verbose => 1) if(defined $opts{'h'});

	pod2usage(q{'-d' vcf dir path must be specified.}) unless(defined $opts{'d'}) ;
	pod2usage(q{'-e' vcf file extension must be specified.}) unless(defined $opts{'e'}) ;
	pod2usage(q{'-p' project name must be specified.}) unless(defined $opts{'p'}) ;
	pod2usage(q{'-o' output location must be specified.}) unless(defined $opts{'o'});

	return \%opts;
}

__END__

=head1 NAME

runVCF2MAF.pl - run vcf2maf conversion 

=head1 SYNOPSIS

runVCF2MAF.pl  -d -e -p -o [ -h -v ]

Required Options (bam and bed interval files must be defined):

  --vcfDir        (-d) vcf dir path [ path to input directory ]
  --vcfExt        (-d) vcf file extension [ e.g., .caveman_c.annot.vcf]
  --project       (-p) project name [e.g MESO ]
  --outdir        (-o) outdir [ Path to output directory ]
  
Optional :
  --help             (-h)  This message
  --version          (-v) displays version number of this software

  Example:
      perl runVCF2MAF.pl -d test caveman_c.vcf -e .caveman_c.vcf  -p MESO -o testdir
=cut




