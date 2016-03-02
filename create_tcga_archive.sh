mv *.vcf sanger.ac.uk_MESO.IlluminaHiSeq_DNASeq.Level_2.1.0.0/
cd sanger.ac.uk_MESO.IlluminaHiSeq_DNASeq.Level_2.1.0.0;

for i in `ls`
	do
		echo `pwd`
		bgzip $i
		tabix -p vcf $i.gz
		md5sum *.gz *.tbi >MANIFEST.txt
 done

cd ..

echo `pwd`


tar -cvzf sanger.ac.uk_MESO.IlluminaHiSeq_DNASeq.Level_2.1.0.0.tar.gz sanger.ac.uk_MESO.IlluminaHiSeq_DNASeq.Level_2.1.0.0
md5sum *.gz >sanger.ac.uk_MESO.IlluminaHiSeq_DNASeq.Level_2.1.0.0.tar.gz.md5 
