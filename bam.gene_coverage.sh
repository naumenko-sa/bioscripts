#get coverage for DMD exons

samtools view -bh DMD.bam X > DMD.X.bam &
bedtools coverage -a ~/Desktop/reference_tables/cheo.muscular_exons.bed.DMD_exons -b DMD.X.bam -mean > DMD.X.bam.coverage
