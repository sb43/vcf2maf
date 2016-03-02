perl create_SDRF.pl project_data/TCGA_MESO_Data_Summary.txt protected >sanger.ac.uk_MESO.IlluminaHiSeq_DNASeq_Cont.1.1.0.sdrf.txt
perl create_SDRF.pl project_data/TCGA_MESO_Data_Summary.txt somatic >sanger.ac.uk_MESO.IlluminaHiSeq_DNASeq.1.1.0.sdrf.txt
cp sanger.ac.uk_MESO.IlluminaHiSeq_DNASeq.1.1.0.sdrf.txt ~/copy_mac/magetab/sanger.ac.uk_MESO.IlluminaHiSeq_DNASeq.mage-tab.1.1.0/
cp sanger.ac.uk_MESO.IlluminaHiSeq_DNASeq_Cont.1.1.0.sdrf.txt ~/copy_mac/magetab/sanger.ac.uk_MESO.IlluminaHiSeq_DNASeq_Cont.mage-tab.1.1.0/
