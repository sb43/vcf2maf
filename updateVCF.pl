#!/usr/bin/env perl

# updted VCF headers - Convert a VCF into a MAF by mapping each variant to only one of all possible gene isoforms

use strict;
use Data::Dumper;
use warnings;
use Capture::Tiny qw(:all);
use List::Util qw(first);

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
	#*****	Change var type 			
					my $var_type='indel'; # snv/indel
					
					#TCGA-YS-AA4M-01A-11D-A39R-32 
					my $tcga_id=join('-',(split '-', $sample_names->{$sample}[1])[0,1,2]);
					my $tcga_file_vcf_name="$tcga_id.sanger.ac.uk.$var_type.vcf";
					my $tcga_vcf_file="$outdir/$tcga_file_vcf_name";;
					
					
				 #my $cmd="bgzip $tcga_vcf_file";
				 #run_cmd($cmd);
				 #$cmd="tabix -p vcf $tcga_vcf_file.gz";
				 #run_cmd($cmd);
				 
					#next;
			
					open (my $tcga_vcf_file_fh, '>', $tcga_vcf_file) or warn "Unable to open file $!";
					
					$sample_set=~s/(.*)\w{1}/$1/g;
	#***** provide algorithm specific file names
					#caveman
					#my $vcf_file_tmp="/lustre/scratch112/sanger/cgppipe/nst_pipe/data/export/request_225/results/$sample_set/tumour/caveman/$sample_names->{$sample}[1]/$sample_names->{$sample}[1]\_vs_$sample_names->{$normal_id}[1].flagged.muts.annot.vcf.gz";
					#pindel
					my $vcf_file_tmp="/lustre/scratch112/sanger/cgppipe/nst_pipe/data/export/request_225/results/$sample_set/tumour/pindel/$sample_names->{$sample}[1]/$sample_names->{$sample}[1]\_vs_$sample_names->{$normal_id}[1].flagged.annot.vcf.gz";
					
	#***** change  #caveman, pindel 
	
					my $algo="pindel";
					
					#caveman :1.7.0;  pindel: 1.5.2 ,
	#*****
					my $software_version='<1.5.2>';
					
					my $vcf = Vcf->new(file => $vcf_file_tmp);	
					$vcf->recalc_ac_an(0);				
					$vcf->parse_header();
					my($process_logs)=get_process_log($vcf->format_header());	
					foreach my $proc_log(@$process_logs) {
						$vcf->remove_header_line(key=>$proc_log);
					}	
					if ($algo eq 'pindel') {
						
					
						my $prm={unmatched=>'pindel_np.gff3.gz',simrep=>'simpleRepeats.bed.gz',rules=>'pulldownRules.lst',annot=>'codingexon_regions.indel.bed.gz'};
						$vcf->add_header_line({key=>'vcfProcessLog',
						InputVCFSource => $tcga_file_vcf_name,
						InputVCFVer => $software_version,
						InputVCFParam => $prm});
					}
			 
				  #print	$vcf->format_header();
					#exit;			
					#print Dumper $vcf->get_header_line(key=>'vcfProcessLog_20150815.1');
					# pindel version 1.5.2 , caveman 1.7.0
					$vcf->add_header_line({key =>'tcgaversion', value=>'1.2'});
					$vcf->add_header_line({key =>'center', value=>'WTSI'});
					$vcf->add_header_line({key =>'PEDIGREE', Name_0=>$sample_names->{$sample}[1],Name_1=>'NA'});
					$vcf->add_header_line({key =>'phasing', value=>'none'});
					$vcf->add_header_line({key =>'assembly', value=> 'ftp://ftp.sanger.ac.uk/pub/cancer/support-files/reference/GRCh37d5.fa'});
					$vcf->add_header_line({key =>'geneAnno', value=> 'ftp://ftp.sanger.ac.uk/pub/cancer/support-files/VAGrENT/reference_data/Homo_sapiens.GRCh37.Ensembl75.VAGrENT_REFERENCE.tar.gz'});
					
					$vcf->remove_header_line(key=>'reference');
					$vcf->add_header_line({	key => 'reference' , value => 'ftp://ftp.sanger.ac.uk/pub/cancer/support-files/reference/GRCh37d5.fa'});
					
					$vcf->add_header_line({	key=>'SAMPLE', ID=>'NORMAL',
						SampleTCGABarcode => $sample_names->{$normal_id}[1],
						SampleUUID 				=> lc($sample_names->{$normal_id}[0]),
						Description 			=> 'NORMAL',
						softwareName 			=> "<$algo>",
						softwareVer 			=> $software_version,
						Platform					=> "Illumina",
						SequenceSource    =>  "WXS",
						Source        		=> 'CGHub',
						SampleName        => $sample_names->{$normal_id}[1],
						Individual        => lc($sample_names->{$normal_id}[0]),
						Accession         => 'NA',
						File							=> $sample_names->{$normal_id}[1].'.bam',
					});
					
					
					$vcf->add_header_line({	key=>'SAMPLE', ID=>'TUMOUR',
						SampleTCGABarcode => $sample_names->{$sample}[1],
						SampleUUID 				=> lc($sample_names->{$sample}[0]),
						Description 			=> 'TUMOUR',
						softwareName 			=>  "<$algo>",
						softwareVer 			=>  $software_version,
						Platform					=> "Illumina",
						SequenceSource    =>  "WXS",
						Source        		=> 'CGHub',
						SampleName        => $sample_names->{$sample}[1],
						Individual        => lc($sample_names->{$sample}[0]),
						Accession         => 'NA',
						File							=> $sample_names->{$sample}[1].'.bam'		
					});
					
					
					$vcf->remove_header_line(key=>'INFO', ID=>'VT');
					$vcf->remove_header_line(key=>'INFO', ID=>'RE');
					$vcf->remove_header_line(key=>'INFO', ID=>'DP');
					
					$vcf->add_header_line({key=>'INFO', ID=>"REn", Number=> "1", Type=> "Integer", Description=>"Range end"});
					$vcf->add_header_line({key=>'INFO', ID=>'VGT', Number => "1", Type=>"String", Description => "Variant type based on the Vagrent Default Annotation" });
					$vcf->add_header_line({key=>'INFO', ID=>'VLSC',Number=>"1",Type=>"Integer", Description=>"Final somatic score between 0 and 255 when multiple lines of evidence are available"});
					$vcf->add_header_line({key=>'INFO', ID=>'VT', Number => "1", Type=>"String", Description => "Variant type, can be SNP, INS , DEL, DNP, TNP, ONP or Consolidated" });
				  $vcf->add_header_line({key=>'INFO',ID=>'VLS', Number=>"1",Type=>'Integer',Description=>"Final validation status relative to non-adjacent Normal,0=wildtype,1=germline,2=somatic,3=LOH,4=post transcriptional modification,5=unknown"});
					$vcf->add_header_line({key=>'INFO',ID=>'GENE',Number=>'.',Type=>'String',Description=>"HUGO gene symbol or Unknown"});
				  $vcf->add_header_line({key=>'INFO',ID=>'RGN', Number=>'.',Type=>'String',Description=>"Region where nucleotide variant occurs in relation to a gene"});
				  $vcf->add_header_line({key=>'INFO',ID=>'RE',Number=>0,Type=>'Flag',Description=>"Position known to have RNA-edits to occur"});
				  $vcf->add_header_line({key=>'INFO',ID=>"DP", Number=> "1", Type=> "Integer", Description=>"Total Depth across samples"});
				  
				  $vcf->add_header_line({key=>'FORMAT',ID=>'DP', Number=>"1",Type=>"Integer", Description=>"Read depth at this position in the sample"});
					$vcf->add_header_line({key=>'FORMAT',ID=>'SS', Number=>"1",Type=>'Integer',Description=>"Variant status relative to non-adjacent Normal,0=wildtype,1=germline,2=somatic,3=LOH,4=post-transcriptional modification,5=unknown"});
					$vcf->add_header_line({key=>'FORMAT',ID=>'BQ', Number=>'.', Type=>'Integer', Description=>"Average base quality for reads supporting alleles"});
					$vcf->add_header_line({key=>'FORMAT',ID=>'AD',Number=>'.',Type=>'Integer', Description=>"Depth of reads supporting alleles 0/1/2/3..." });
					
				
				  my $tmp=$vcf->format_header();
				  my @header_data=split('\n',$tmp);
				  foreach my $line (@header_data) {
						if($line=~/vcfProcessLog/) {	
							$line=~s/<|>//g;	
							$line=~s/(vcfProcessLog.*?=)/$1</g;
							$line=~s/InputVCFParam=/InputVCFParam="/g;
						  print "$line\">\n";
						  print $tcga_vcf_file_fh "$line\">\n";
						}
						else {
							print $tcga_vcf_file_fh "$line\n";
						}
					}
				  
				  #exit;
				  
				  my $ss_val=5;
				  
				  while (my $x = $vcf->next_data_hash()) {
				    my $ref=$x->{'REF'};
					  my $alt=$x->{'ALT'}[0];
				  	next if ($ref =~ m/N/ || $alt =~ m/N/);
				  	
				    my $info_hash=$x->{INFO};
				    
				   	$info_hash->{VGT}=$info_hash->{VT};
				   	
				   	delete $info_hash->{VT};
				   	
				   	
				   			    
						$vcf->add_format_field($x,'DP');
						$vcf->add_format_field($x,'SS');
						$vcf->add_format_field($x,'BQ');
						$vcf->add_format_field($x,'AD');
									
						if(first {$_ eq 'PASS'} @{$x->{FILTER}}) {
              	$ss_val='2';
            }
						
					  if(first {$_ eq 'PR'} @{$x->{FORMAT}}) { # Pindel VCF
					  
					    my ( $ref_length, $var_length ) = ( length( $ref ), length( $alt ));
					    my @alleles = ( $ref, split( /,/, $alt ));
							# Remove the prefixed reference bp from all alleles, using "-" for simple indels
								( $ref, $alt, @alleles ) = map{$_ = substr( $_, 1 ); ( $_ ? $_ : "-" )} ( $ref, $alt, @alleles );
								--$ref_length; --$var_length;
								if( $ref_length < $var_length ) { # Handle insertions, and the special case for complex ones
									$var_type = "INS";
								}
								else { # Handle deletions
									$var_type = "DEL";
								}
							
					     $info_hash->{'VT'}=$var_type;
					    
					    $info_hash->{'REn'}=$info_hash->{'RE'};
					    delete $info_hash->{'RE'};
					    
					  	$$x{gtypes}{'TUMOUR'}{'AD'}=( ($$x{gtypes}{'TUMOUR'}{'PR'} + $$x{gtypes}{'TUMOUR'}{'NR'}) - ($$x{gtypes}{'TUMOUR'}{'NP'} + $$x{gtypes}{'TUMOUR'}{'PP'}).','.($$x{gtypes}{'TUMOUR'}{'NP'} + $$x{gtypes}{'TUMOUR'}{'PP'})  );
   					  $$x{gtypes}{'TUMOUR'}{'DP'}=( ($$x{gtypes}{'TUMOUR'}{'PR'} + $$x{gtypes}{'TUMOUR'}{'NR'}) );
   					  
   					  $$x{gtypes}{'NORMAL'}{'AD'}=( ($$x{gtypes}{'NORMAL'}{'PR'} + $$x{gtypes}{'NORMAL'}{'NR'}) - ($$x{gtypes}{'NORMAL'}{'NP'} + $$x{gtypes}{'NORMAL'}{'PP'}).','.($$x{gtypes}{'NORMAL'}{'NP'} + $$x{gtypes}{'NORMAL'}{'PP'})  );
   					  $$x{gtypes}{'NORMAL'}{'DP'}=( ($$x{gtypes}{'NORMAL'}{'PR'} + $$x{gtypes}{'NORMAL'}{'NR'}) );
   					  
   					  $$x{gtypes}{'TUMOUR'}{'SS'}=$ss_val;
   					  $$x{gtypes}{'NORMAL'}{'SS'}=5;
   					  $$x{gtypes}{'TUMOUR'}{'BQ'}='.';
   					  $$x{gtypes}{'NORMAL'}{'BQ'}='.';
   					  $info_hash->{'DP'}=$$x{gtypes}{'NORMAL'}{'DP'}+$$x{gtypes}{'TUMOUR'}{'DP'};
   					  
					  }
						
					  if(first {$_ eq 'FAZ'} @{$x->{FORMAT}}) { # Caveman VCF
					    $info_hash->{'VT'}='SNP';
					    #delete $info_hash->{'DP'};
					  	$$x{gtypes}{'TUMOUR'}{'AD'}=( ($$x{gtypes}{'TUMOUR'}{'F'.$ref.'Z'} + $$x{gtypes}{'TUMOUR'}{'R'.$ref.'Z'}).','.($$x{gtypes}{'TUMOUR'}{'F'.$alt.'Z'} + $$x{gtypes}{'TUMOUR'}{'R'.$alt.'Z'})  );
              $$x{gtypes}{'TUMOUR'}{'DP'}= $$x{gtypes}{'TUMOUR'}{'FAZ'} + $$x{gtypes}{'TUMOUR'}{'FCZ'} + $$x{gtypes}{'TUMOUR'}{'FGZ'} + $$x{gtypes}{'TUMOUR'}{'FTZ'} + $$x{gtypes}{'TUMOUR'}{'RAZ'} + $$x{gtypes}{'TUMOUR'}{'RCZ'} + $$x{gtypes}{'TUMOUR'}{'RGZ'} + $$x{gtypes}{'TUMOUR'}{'RTZ'};
							
							$$x{gtypes}{'NORMAL'}{'AD'}=( ($$x{gtypes}{'NORMAL'}{'F'.$ref.'Z'} + $$x{gtypes}{'NORMAL'}{'R'.$ref.'Z'}).','.($$x{gtypes}{'NORMAL'}{'F'.$alt.'Z'} + $$x{gtypes}{'NORMAL'}{'R'.$alt.'Z'})  );
              $$x{gtypes}{'NORMAL'}{'DP'}= $$x{gtypes}{'NORMAL'}{'FAZ'} + $$x{gtypes}{'NORMAL'}{'FCZ'} + $$x{gtypes}{'NORMAL'}{'FGZ'} + $$x{gtypes}{'NORMAL'}{'FTZ'} + $$x{gtypes}{'NORMAL'}{'RAZ'} + $$x{gtypes}{'NORMAL'}{'RCZ'} + $$x{gtypes}{'NORMAL'}{'RGZ'} + $$x{gtypes}{'NORMAL'}{'RTZ'};
              
              $$x{gtypes}{'TUMOUR'}{'SS'}=$ss_val;
   					  $$x{gtypes}{'NORMAL'}{'SS'}=5;
   					  $$x{gtypes}{'TUMOUR'}{'BQ'}='.';
   					  $$x{gtypes}{'NORMAL'}{'BQ'}='.';
   					  
					  }
            $x->{INFO}=$info_hash;
						#print Dumper $x;
						print $tcga_vcf_file_fh $vcf->format_line($x);       	          	
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


sub get_process_log {
	my ($header)=@_;
	my @data=split('\n',$header);
	my $process_logs;
	foreach my $line (@data) {
		if($line=~/vcfProcessLog_/) {	
		  $line=~s/##//g;
			push(@$process_logs, (split '=',$line)[0]);
		}
	}
	
	return $process_logs;
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
  --vcfExt        (-e) vcf file extension [ e.g., .caveman_c.annot.vcf, .pindel.annot.vcf]
  --project       (-p) project name [e.g MESO ]
  --outdir        (-o) outdir [ Path to output directory ]
  
Optional :
  --help             (-h)  This message
  --version          (-v) displays version number of this software

  Example:
      perl runVCF2MAF.pl -d test caveman_c.vcf -e .caveman_c.vcf  -p MESO -o testdir
=cut




