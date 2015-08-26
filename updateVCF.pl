#!/usr/bin/env perl

# updted VCF headers - Convert a VCF into a MAF by mapping each variant to only one of all possible gene isoforms

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
	my $i=0;
	foreach my $sample (keys %$sample_names) {
		if( ref($sample_names->{$sample}) eq 'ARRAY' ) {
			my $input_vcf=$file_path.$sample.$file_extension;
			if( -s $input_vcf.'.gz' ) {
					my $vcf_old = Vcf->new(file => $input_vcf.'.gz');
					$vcf_old->parse_header();	
					$vcf_old->recalc_ac_an(0);
					my($normal_id)=get_normal_id($vcf_old);
					$vcf_old->close();
					#print "$sample---$normal_id------$sample_names->{$sample}[0]-----$sample_names->{$sample}[1]------$sample_names->{$normal_id}[0]-----$sample_names->{$normal_id}[1]\n";
					#my $cmd = "bgzip -f -c -d $input_vcf.gz >$input_vcf";
					#run_cmd($cmd);
					my $sample_set=$sample;
	#*****				
					my $var_type='Indel'; # Sub/Indel
					
					my $tcga_vcf_file="$outdir/$sample_names->{$sample}[1]\_vs_$sample_names->{$normal_id}[1].flagged.$var_type.annot.vcf";
					
					
				 #my $cmd="bgzip $tcga_vcf_file";
				 #run_cmd($cmd);
				 #$cmd="tabix -p vcf $tcga_vcf_file.gz";
				 #run_cmd($cmd);
				 
					#next;
					
					open (my $tcga_vcf_file_fh, '>', $tcga_vcf_file) or warn "Unable to open file $!";
					
					$sample_set=~s/(.*)\w{1}/$1/g;
	#*****
					#caveman
					#my $vcf_file_tmp="/lustre/scratch112/sanger/cgppipe/nst_pipe/data/export/request_225/results/$sample_set/tumour/caveman/$sample_names->{$sample}[1]/$sample_names->{$sample}[1]\_vs_$sample_names->{$normal_id}[1].flagged.muts.annot.vcf.gz";
					#pindel
					my $vcf_file_tmp="/lustre/scratch112/sanger/cgppipe/nst_pipe/data/export/request_225/results/$sample_set/tumour/pindel/$sample_names->{$sample}[1]/$sample_names->{$sample}[1]\_vs_$sample_names->{$normal_id}[1].flagged.annot.vcf.gz";
					#CaVEMan, Pindel
	#*****
					my $software_name='<pindel>';
					#caveman :1.7.0;  pindel: 1.5.2 ,
	#*****
					my $software_version='<1.5.2>';
					
					my $vcf = Vcf->new(file => $vcf_file_tmp);					
					my $header=$vcf->parse_header();
										
					#print Dumper $vcf->get_header_line(key=>'vcfProcessLog_20150815.1');
					# pindel version 1.5.2 , caveman 1.7.0
					$vcf->add_header_line({key =>'tcgaversion', value=>'1.2'});
					$vcf->remove_header_line(key=>'reference');
					$vcf->add_header_line({	key => 'reference' , value => 'ftp://ftp.sanger.ac.uk/pub/cancer/support-files/reference/GRCh37d5.fa'});
					
					$vcf->add_header_line({	key=>'SAMPLE', ID=>'NORMAL',
						SampleTCGABarcode => $sample_names->{$normal_id}[1],
						SampleUUID 				=> lc($sample_names->{$normal_id}[0]),
						Description 			=> 'NORMAL',
						softwareName 			=> $software_name,
						softwareVer 			=> $software_version,
						Platform					=> "Illumina",
						SequenceSource    =>  "WXS",
						Source        		=> 'CGHub',
						SampleName        => $sample_names->{$normal_id}[1],
						File							=> $sample_names->{$normal_id}[1].'.bam',
					});
					
					
					$vcf->add_header_line({	key=>'SAMPLE', ID=>'TUMOUR',
						SampleTCGABarcode => $sample_names->{$sample}[1],
						SampleUUID 				=> lc($sample_names->{$sample}[0]),
						Description 			=> 'TUMOUR',
						softwareName 			=> $software_name,
						softwareVer 			=>  $software_version,
						Platform					=> "Illumina",
						SequenceSource    =>  "WXS",
						Source        		=> 'CGHub',
						SampleName        => $sample_names->{$sample}[1],
						File							=> $sample_names->{$sample}[1].'.bam'		
					});
					
				  print $tcga_vcf_file_fh $vcf->format_header();
					while (my $x=$vcf->next_line()) { 
								my ($ref,$alt)=(split /\t/ ,$x)[3,4];
								next if ($ref =~ m/N/ || $alt =~ m/N/);
								print $tcga_vcf_file_fh $x;
        		
    			}

			}		
		}
	}




#"bsub -oo logs/caveaman_pass%I.log  -q normal -J 'caveman_pass[1-83]' -P analysis-cgp -n 1 -R'select[mem>=500] span[hosts=1] rusage[mem=500]' -M 500  'perl farm_idx_exec.pl caveman_pass.cmd $LSB_JOBINDEX'"

#print "bsub -oo logs/test%I.log  -q normal -J \'UNM[1-$i]\' -P analysis-cgp -n 1 -R\'select[cgp_nfs>=13] span[hosts=1] rusage[cgp_nfs=13]\' -M 500  \'perl farm_idx_exec.pl UNM.cmd \$LSB_JOBINDEX\'";
	
	#print "completed analysis -- time to merge";
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
  --vcfExt        (-e) vcf file extension [ e.g., .caveman_c.annot.vcf]
  --project       (-p) project name [e.g MESO ]
  --outdir        (-o) outdir [ Path to output directory ]
  
Optional :
  --help             (-h)  This message
  --version          (-v) displays version number of this software

  Example:
      perl runVCF2MAF.pl -d test caveman_c.vcf -e .caveman_c.vcf  -p MESO -o testdir
=cut




