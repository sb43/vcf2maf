#!/usr/bin/perl

use strict;

use Data::Dumper;

my($tcga_file,$type)=@ARGV;

if(!defined $type || !defined $tcga_file) {
	warn "Please provide TCGA summary file and file type to generate[protected/somatic]";
	exit;
}
open my $tcga_fh, '<', $tcga_file || warn "unable to open input file $!";

# generated following header using data from : and command:https://wiki.nci.nih.gov/display/TCGA/Sample+and+Data+Relationship+Format
#~/scripts/vcf2maf/project_data$ cut -f 2 sdrf.headers |xargs -I {} echo \"{}\"\,|xargs
my @header=(
#Material
'Extract Name','Comment [TCGA Barcode]','Comment [is tumor]', 'Material Type', 'Annotation REF','Comment [TCGA Genome Reference]',
#Lib preparation
'Protocol REF','Parameter Value [Vendor]','Parameter Value [Catalog Name]','Parameter Value [Catalog Number]','Parameter Value [Annotation URL]','Parameter Value [Product URL]','Parameter Value [Target File URL]','Parameter Value [Target File Format]','Parameter Value [Target File Format Version]','Parameter Value [Probe File URL]','Parameter Value [Probe File Format]','Parameter Value [Probe File Format Version]','Parameter Value [Target Reference Accession]',
#sequencing
'Protocol REF',
#mapping
'Protocol REF','Comment [Derived Data File REF]','Comment [TCGA CGHub ID]','Comment [TCGA CGHub metadata URL]','Comment [TCGA Include for Analysis]',
'Derived Data File','Comment [TCGA Include for Analysis]','Comment [TCGA Data Type]','Comment [TCGA Data Level]','Comment [TCGA Archive Name]',
'Parameter Value [Protocol Min Base Quality]', 'Parameter Value [Protocol Min Map Quality]','Parameter Value [Protocol Min Tumor Coverage]','Parameter Value [Protocol Min Normal Coverage]',
#variant calling
'Protocol REF','Derived Data File','Comment [TCGA Spec Version]','Comment [TCGA Include for Analysis]','Comment [TCGA Data Type]','Comment [TCGA Data Level]','Comment [TCGA Archive Name]',
#MAF generation
'Protocol REF','Derived Data File','Comment [TCGA Spec Version]','Comment [TCGA Include for Analysis]','Comment [TCGA Data Type]','Comment [TCGA Data Level]','Comment [TCGA Archive Name]',
# Mutation validation
'Protocol REF','Derived Data File','Comment [TCGA Spec Version]','Comment [TCGA Include for Analysis]','Comment [TCGA Data Type]','Comment [TCGA Data Level]','Comment [TCGA Archive Name]',
);

# for reference only
my @tcga_header=qw/study	barcode	disease	disease_name	sample_type	sample_type_name	analyte_type	library_type	center	center_name	platform[10]	platform_name	assembly	filename	files_size	checksum	analysis_id[16]	aliquot_id	participant_id[18]	sample_id	tss_id	sample_accession	published	uploaded	modified	state	reason/;
print join("\t",@header)."\n";

my $count=0;
while (<$tcga_fh>) {
	chomp;
	$count++;
	next if $count < 2;


	# Material	
	my($Extract_Name,$Comment_tcga_barcode,$is_tumour,$Material_Type,$Genome_reference, )=(split "\t", $_)[17,1,4,6,12];
	my $Annotation_REF = "->";
	if($is_tumour eq 'TP') {$is_tumour='yes';}
	else{$is_tumour='no';}

	# Library
	
	my $LibProtocol_REF='sanger.ac.uk:library_preparation:Illumina:01';
	my $Vendor='Agilent';
	my $Catalog_Name='SureSelect Human All Exon V3';
	my $Catalog_Number='->';
	my $Annotation_URL='->';
	my $Product_URL='->';
	my $Target_File_URL='->';
	my $Target_File_Format='->';
	my $Target_File_Format_Version='->';
	my $Probe_File_URL='->';
	my $Probe_File_Format='->';
	my $Probe_File_Format_Version='->';
	my $Target_Reference_Accession='->';

	#sequencing
	my $SEQ_Protocol_REF='sanger.ac.uk:DNA_sequencing:Illumina:01';

	#Mapping	
	my $MAP_Protocol_REF='sanger.ac.uk:alignment:bwa-mem:01';
	my $MAP_Derived_Data_File_REF="$Comment_tcga_barcode.bam";
	my $TCGA_CGHub_ID='->';
	my $TCGA_CGHub_metadata_URL='->';
	my $MAP_TCGA_Include_for_Analysis='yes';
	my $MAP_Derived_Data_File='->';
	my $MAP_TCGA_Include_for_Analysis_BIGWIG='->';
	my $MAP_TCGA_Data_Type='->';
	my $MAP_TCGA_Data_Level='->';
	my $MAP_TCGA_Archive_Name='->';
	my $Min_Base_Quality='->';
	my $Min_Map_Quality='->';
	my $Protocol_Min_Tumor_Coverage='->';
	my $Protocol_Min_Normal_Coverage='->';



	#VCF	

	my $tcga_id=join('-',(split '-', $Comment_tcga_barcode)[0,1,2]);
	my $tcga_file_vcf_name="$tcga_id.sanger.ac.uk.\%s.vcf";

	my $VCF_Protocol_REF='->';
	my $VCF_Derived_Data_File='->';
	my $VCF_TCGA_Spec_Version='->';
	my $VCF_TCGA_Include_for_Analysis='->';
	my $VCF_TCGA_Data_Type='->';
	my $VCF_TCGA_Data_Level='->';
	my $VCF_TCGA_Archive_Name='->';

	if($type eq 'protected') {
		$VCF_Protocol_REF='sanger.ac.uk:variant_calling:caveman_and_pindel:01';
		$VCF_Derived_Data_File=$tcga_file_vcf_name if ($is_tumour eq 'yes');
		$VCF_TCGA_Spec_Version='4.1';
		$VCF_TCGA_Archive_Name='sanger.ac.uk_MESO.IlluminaHiSeq_DNASeq_automated_Cont.Level_2.1.1.0';
		$VCF_TCGA_Include_for_Analysis='yes';
		$VCF_TCGA_Data_Type='Mutations';
		$VCF_TCGA_Data_Level='Level 2';
	}
	#MAF

	my $MAF_Protocol_REF='sanger.ac.uk:vcf2maf:data_consolidation:01';
	my $MAF_Derived_Data_File='->';
	if(($type eq 'protected') && ($is_tumour eq 'yes' )) {
		 $MAF_Derived_Data_File="sanger.ac.uk_MESO.IlluminaHiSeq_DNASeq_automated_Cont.Level_2.1.1.0.$type.maf";
	}
	elsif($is_tumour eq 'yes' ) {
		$MAF_Derived_Data_File="sanger.ac.uk_MESO.IlluminaHiSeq_DNASeq_automated.Level_2.1.1.0.$type.maf";
	}

	my $MAF_TCGA_Spec_Version='2.3';
	my $MAF_TCGA_Include_for_Analysis='yes';
	my $MAF_TCGA_Data_Type='Mutations';
	my $MAF_TCGA_Data_Level='Level 2';

	my $MAF_TCGA_Archive_Name='sanger.ac.uk_MESO.IlluminaHiSeq_DNASeq_automated.Level_2.1.1.0';
	if($type eq 'protected') {
	 $MAF_TCGA_Archive_Name='sanger.ac.uk_MESO.IlluminaHiSeq_DNASeq_automated_Cont.Level_2.1.1.0';
	}
	#VALIDATION
	my $MUT_Protocol_REF='->';
	my $MUT_Derived_Data_File='->';
	my $MUT_TCGA_Spec_Version='->';
	my $MUT_TCGA_Include_for_Analysis='->';
	my $MUT_TCGA_Data_Type='->';
	my $MUT_TCGA_Data_Level='->';
	my $MUT_TCGA_Archive_Name='->';

	my $data={ 
	#Material	
	1=>$Extract_Name, 2=>$Comment_tcga_barcode, 3=>$is_tumour, 4=>$Material_Type, 5=>$Annotation_REF, 6=>$Genome_reference,

	#Library protocol
	7=>$LibProtocol_REF , 8=>$Vendor , 10=>$Catalog_Name,
	11=>$Catalog_Number,
	12=>$Annotation_URL,
	13=>$Product_URL,
	14=>$Target_File_URL,
	15=>$Target_File_Format,
	16=>$Target_File_Format_Version,
	17=>$Probe_File_URL,
	18=>$Probe_File_Format,
	19=>$Probe_File_Format_Version,
	20=>$Target_Reference_Accession,
	# sequencing
	31=>$SEQ_Protocol_REF,
	# mapping
	614=>$MAP_Protocol_REF,
	615=>$MAP_Derived_Data_File_REF,
	616=>$TCGA_CGHub_ID,
	617=>$TCGA_CGHub_metadata_URL,
	618=>$MAP_TCGA_Include_for_Analysis,
	619=>$MAP_Derived_Data_File,
	620=>$MAP_TCGA_Include_for_Analysis_BIGWIG,
	621=>$MAP_TCGA_Data_Type,
	622=>$MAP_TCGA_Data_Level,
	623=>$MAP_TCGA_Archive_Name,
	624=>$Min_Base_Quality,
	625=>$Min_Map_Quality,
	626=>$Protocol_Min_Tumor_Coverage,
	627=>$Protocol_Min_Normal_Coverage,

	#VCF
	711=>$VCF_Protocol_REF,
	712=>$VCF_Derived_Data_File,
	713=>$VCF_TCGA_Spec_Version,
	714=>$VCF_TCGA_Include_for_Analysis,
	715=>$VCF_TCGA_Data_Type,
	716=>$VCF_TCGA_Data_Level,
	717=>$VCF_TCGA_Archive_Name,

	#MAF
	811=>$MAF_Protocol_REF,
	812=>$MAF_Derived_Data_File,
	813=>$MAF_TCGA_Spec_Version,
	814=>$MAF_TCGA_Include_for_Analysis,
	815=>$MAF_TCGA_Data_Type,
	816=>$MAF_TCGA_Data_Level,
	817=>$MAF_TCGA_Archive_Name,
	# Mutation_validation
	911=>$MUT_Protocol_REF,
	912=>$MUT_Derived_Data_File,
	913=>$MUT_TCGA_Spec_Version,
	914=>$MUT_TCGA_Include_for_Analysis,
	915=>$MUT_TCGA_Data_Type,
	916=>$MUT_TCGA_Data_Level,
	917=>$MUT_TCGA_Archive_Name,

	};


	my $line=undef;
	foreach my $key (sort {$a<=>$b} keys %$data) {
		
				if($data->{$key}){
					$line.="$data->{$key}\t";
				}
				else{
					$line.='->'."\t";	
				}
			}
	$line=~s/\t$//g;
	if($is_tumour eq 'yes') {
		
		my $tmp=$line;
		
		my $line_snv=sprintf($line,'snv');
		my $line_indel=sprintf($tmp,'indel');
		print "$line_snv\n$line_indel\n";
	}
	else {
		print "$line\n";
	}
 


}



