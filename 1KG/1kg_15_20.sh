#!/bin/bash

# download large genotype files and split them by sample, extract first 100 samples

# benchmark: chr22: 26G download - 1h

date
# download the reference - needed for normalization
wget -c https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
wget -c https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai

for chr in `seq 15 20`
do
#chr=21

    echo `date` $chr

    if [ ! -f chr${chr}.clean.vcf.gz ];then
        wget -c https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.recalibrated_variants.vcf.gz -O chr${chr}.vcf.gz
        wget -c https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.recalibrated_variants.vcf.gz.tbi -O chr${chr}.vcf.gz.tbi

        # remove not PASS
        # remove annotation to make the file smaller (does not make a file lot smaller maybe 2G, 
        # but otherwise they will go into every individual sample)
        # split multiallelic and normalise
        # piped together with -Ou it is actually very fast
        bcftools view -f "PASS,." -Ou chr${chr}.vcf.gz | bcftools annotate -x INFO/AC,INFO/AF,INFO/AN,INFO/AN_EUR,INFO/AN_EAS,INFO/AN_AMR,INFO/AN_SAS,INFO/AN_AFR,INFO/AF_EUR,INFO/AF_EAS,INFO/AF_AMR,INFO/AF_SAS,INFO/AF_AFR,INFO/AC_EUR,INFO/AC_EAS,INFO/AC_AMR,INFO/AC_SAS,INFO/AC_AFR,INFO/AC_Het_EUR,INFO/AC_Het_EAS,INFO/AC_Het_AMR,INFO/AC_Het_SAS,INFO/AC_Het_AFR,INFO/AC_Het,INFO/AC_Hom_EUR,INFO/AC_Hom_EAS,INFO/AC_Hom_AMR,INFO/AC_Hom_SAS,INFO/AC_Hom_AFR,INFO/AC_Hom,INFO/HWE_EUR,INFO/ExcHet_EUR,INFO/HWE_EAS,INFO/ExcHet_EAS,INFO/HWE_AMR,INFO/ExcHet_AMR,INFO/HWE_SAS,INFO/ExcHet_SAS,INFO/HWE_AFR,INFO/ExcHet_AFR,INFO/AN_EUR_unrel,INFO/AN_EAS_unrel,INFO/AN_AMR_unrel,INFO/AN_SAS_unrel,INFO/AN_AFR_unrel,INFO/AF_EUR_unrel,INFO/AF_EAS_unrel,INFO/AF_AMR_unrel,INFO/AF_SAS_unrel,INFO/AF_AFR_unrel,INFO/AC_EUR_unrel,INFO/AC_EAS_unrel,INFO/AC_AMR_unrel,INFO/AC_SAS_unrel,INFO/AC_AFR_unrel,INFO/AC_Het_EUR_unrel,INFO/AC_Het_EAS_unrel,INFO/AC_Het_AMR_unrel,INFO/AC_Het_SAS_unrel,INFO/AC_Het_AFR_unrel,INFO/AC_Hom_EUR_unrel,INFO/AC_Hom_EAS_unrel,INFO/AC_Hom_AMR_unrel,INFO/AC_Hom_SAS_unrel,INFO/AC_Hom_AFR_unrel -Ou - | \
             bcftools norm -m -both -f GRCh38_full_analysis_set_plus_decoy_hla.fa -Ov - | bgzip -c > chr${chr}.clean.vcf.gz
        tabix chr${chr}.clean.vcf.gz
    fi

    bcftools query -l chr${chr}.clean.vcf.gz | head -n100 > samples100.txt
    bcftools +split -S samples100.txt -o . -i'GT="alt"' chr${chr}.clean.vcf.gz

    for f in `cat samples100.txt`;do mv $f.vcf $f.chr${chr}.vcf;bgzip $f.chr${chr}.vcf;tabix $f.chr${chr}.vcf.gz;done;

    rm chr$chr.clean.vcf.gz chr$chr.clean.vcf.gz.tbi chr$chr.vcf.gz chr$chr.vcf.gz.tbi
done

date
