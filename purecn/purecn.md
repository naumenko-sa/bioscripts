# PureCN

## 1. Create SNV PON with Mutect2

- This step could be done once for a sequencing panel with > 40 normal samples.
- see [gatk documentation]
- https://github.com/lima1/PureCN/issues/117
- https://gatk.broadinstitute.org/hc/en-us/articles/360037058172-CreateSomaticPanelOfNormals-BETA-
- https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2)
- https://gatkforums.broadinstitute.org/gatk/discussion/11136/how-to-call-somatic-mutations-using-gatk4-mutect2#2
- progress with gatk4 in purecn: https://github.com/lima1/PureCN/issues/6
- https://github.com/lima1/PureCN/issues/117

- Use *normal* samples only, call variants in *tumor-only* mode with mutect2
- bcbio does not support yet running mutect2 in tumor-only mode without a PON, use commands from GATK site.

### 1.1 Create interval_list file:
```
gatk BedToIntervalList \
-I panel.bed \
-O panel.interval_list \
-SD /bcbio/genomes/Hsapiens/hg38/seq/hg38.dict
```

### 1.2 Download gnomad-af only vcf:
```
wget -c https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz
wget -c https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi
```

### 1.3 Call variants for each normal bam in tumor-only mode
```
gatk Mutect2 \
-R /bcbio/genomes/Hsapiens/hg38/seq/hg38.fa \
-I S1_N.bam \
-O S1_N_for_pon.vcf.gz \
-tumor S1_N \
--max-mnp-distance 0 \
--intervals panel.interval_list \
--interval-padding 50 \
--germline-resoure af-only-gnomad.hg38.vcf.gz \
--genotype-germline-sites

tabix S1_N.vcf.gz
```

### 1.4 Merge files with gatk3
Alternative steps 1.4.1-1.4.2 are not working for PureCN,
because CreateSomaticPanelOfNormals does not produce AD field in the vcf.

```
vcf_files=""
for f in *.for_pon.vcf.gz
do
    vcf_files="$vcf_files -V $f"
done

gatk3 -Xmx12g \
-T CombineVariants \
--minimumN 3 \
-R /data/genomes/Hsapiens/hg38/seq/hg38.fa \
-o snv_pon.vcf \
$vcf_files

bgzip snv_pon.vcf
tabix snv_pon.vcf.gz
```

### [future] 1.4.1 Create genomicsdb
```
ls -1 *.for_pon.vcf.gz | awk -F "." '{print $1"\t"$0}' > sample_list.tsv

gatk GenomicsDBImport \
-R /data/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa \
--intervals $1 \
--sample-name-map sample_list.tsv \
--genomicsdb-workspace-path pon_db
```

### [future] 1.4.2 Create snv_pon.vcf.gz
```
gatk CreateSomaticPanelOfNormals \
-R /data/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa \
--germline-resource af-only-gnomad.vcf.gz \
-V gendb://pon_db \
-O snv_pon.vcf.gz
```

## 2 Do segmentation with gatkcnv

call tumor-only variants with Mutect2```
gatk Mutect2 \
-R /bcbio/genomes/Hsapiens/hg38/seq/hg38.fa \
-I S1_N.bam \
-O S1_N.vcf.gz \
--max-mnp-distance 0 \
?use --intervals panel.interval_list \
?    --intterval-padding 50 \
--genotype-germline-sites

## 3. Run pureCN with snv_pon and segments from gatkcnv.

3.1 Prepare intervals
