# PureCN

## 1. Create SNV PON with Mutect2

- This step could be done once for a sequencing with > 40 normal samples.
- see [gatk manual](https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2)
- Use *normal* samples only, call variants in *tumor-only* mode with mutect2
- bcbio does not support yet running mutect2 in tumor-only mode without a PON, use commands from GATK site.

### 1.1 Create interval_list file:
```
gatk BedToIntervalList \
-I panel.bed \
-O panel.interval_list \
-SD /bcbio/genomes/Hsapiens/hg38/seq/hg38.dict
```

### 1.2 Call variants for each normal bam
- using interval padding
- output germline variants - PureCN needs them
```
gatk Mutect2 \
-R /bcbio/genomes/Hsapiens/hg38/seq/hg38.fa \
-I S1_N.bam \
-O S1_N.vcf.gz \
--max-mnp-distance 0 \
```

?use --intervals panel.interval_list \
?    --intterval-padding 50 \
--genotype-germline-sites

### 1.3 Create genomics.db
```
gatk GenomicsDBImport\
-R /bcbio/genomes/Hsapiens/hg38/seq/hg38.fa \
-L intervals.interval_list \
        --genomicsdb-workspace-path pon_db \
        -V normal1.vcf.gz \
        -V normal2.vcf.gz \
        -V normal3.vcf.gz
```

### 1.4 Download gnomad-af only vcf:
```
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz
```

### 1.5. combine calls into PON.vcf
```
gatk CreateSomaticPanelOfNormals -R reference.fasta \
      --germline-resource af-only-gnomad.vcf.gz \
      -V gendb://pon_db \
      -O pon.vcf.gz
      ```



result: snv_pon.vcf.gz


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
