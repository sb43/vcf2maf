#!/usr/bin/env perl

# vcf2maf - Convert a VCF into a MAF by mapping each variant to only one of all possible gene isoforms

use strict;
use Data::Dumper;
use warnings;
use Capture::Tiny qw(:all);

use Vcf;


my $file_path=$ARGV[0];
my $file_extension = $ARGV[1];

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



foreach my $sample (keys %$sample_names) {
	if( ref($sample_names->{$sample}) eq 'ARRAY' ) {
		my $input_vcf=$file_path.'/'.$sample.$file_extension;
		if( -s $input_vcf ) {
				my $vcf = Vcf->new(file => $input_vcf);
				$vcf->parse_header();	
			  $vcf->recalc_ac_an(0);
			  my($normal_id)=get_normal_id($vcf);
			  
			  print "$sample---$normal_id------$sample_names->{$sample}[0]-----$sample_names->{$sample}[1]------$sample_names->{$normal_id}[0]-----$sample_names->{$normal_id}[1]\n";
			  
			  my $cmd = "bgzip -d $input_vcf";
			  run_cmd($cmd);
			  $input_vcf=~s/\.gz//g;
			  $cmd = "perl vcf2maf.pl --input-vagrent $input_vcf ".
			  			" --tumor-id $sample_names->{$sample}[1] --normal-id $sample_names->{$normal_id}[1] 
			  			--vcf-tumor-id TUMOUR --vcf-norma-id NORMAL --tumor-uuid $sample_names->{$sample}[0]  --normal-uuid $sample_names->{$normal_id}[0]".
			  			"--maf-center WTSI ";
			  
			  print $cmd;
			  
			  exit;
			
			 
			
		}
				
	}
	
}




sub get_normal_id {
	my ($vcf)=@_;
	my $sample_info=$vcf->get_header_line(key=>'SAMPLE', ID=>'NORMAL' , );
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
