# PureCN step by step

## 0. Clean up a bed file - remove intervals on ALT chromosomes
```
cat coverage.bed | grep -v alt > panel.bed
```

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
# -tumor is deprecated

gatk Mutect2 \
-R /bcbio/genomes/Hsapiens/hg38/seq/hg38.fa \
-I S1_N.bam \
-O S1_N.vcf.gz \
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
substitute with bcftools?

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

## 2. Create CNV pon
See https://bcbio-nextgen.readthedocs.io/en/latest/contents/somatic_variants.html#workflow3-copy-number-variants

### 2.1 preprocess intervals
```
# $1 = panel.bed or could be interval_list file
# $2 = hg38.fa

bname=`basename $1 .bed`

gatk PreprocessIntervals \
-R /data/genomes/Hsapiens/hg38/seq/hg38.fa \
--interval-merging-rule OVERLAPPING_ONLY \
-O $bname.padding250.interval_list \
-L $1 \
--bin-length 0 \
--padding 250
```

### 2.2 GC annotation
```
# $1 = panel.padding250.interval_list
# $2 = hg38.fa
# output = panel.gcannotated.tsv

bname=`basename $1 .interval_list`

gatk \
--java-options '-Xms500m -Xmx131859m -XX:+UseSerialGC -Djava.io.tmpdir=.' \
AnnotateIntervals \
-R /data/genomes/Hsapiens/hg38/seq/hg38.fa \
-L $1 \
--interval-merging-rule OVERLAPPING_ONLY \
-O $bname.gcannotated.tsv
```

### 2.3 CollectReadCounts
```
# $1 = sort.bam
# $2 = interval.list - not gc annotated

bname=`basename $1 .bam`

gatk --java-options '-Xms500m -Xmx131859m -XX:+UseSerialGC -Djava.io.tmpdir=.' \
CollectReadCounts \
-I $1 \
-L $2 \
--interval-merging-rule OVERLAPPING_ONLY \
-O $bname.counts.hdf5
```

### 2.4 gather a cnv.pon
```
# $1 = panel.gcannotated.tsv

hdf5_files=""
for f in *.hdf5
do
    hdf5_files="$hdf5_files -I $f"
done

unset JAVA_HOME && \
export PATH=/bcbio/anaconda/bin:"$PATH" && \
gatk --java-options '-Xms500m -Xmx131859m -XX:+UseSerialGC -Djava.io.tmpdir=.' \
CreateReadCountPanelOfNormals \
-O cnv.pon.hdf5 \
--annotated-intervals $1 \
$hdf5_files \
--maximum-zeros-in-sample-percentage 100
```

## 3. Call SNV a in a tumor-only sample
```
# $1 = sample_N.bam
# $2 = panel.interval_list

bname=`basename $1 .bam`

gatk Mutect2 \
-R /data/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa \
-I $1 \
-O $bname.vcf.gz \
--max-mnp-distance 0 \
--intervals $2 \
--interval-padding 50 \
--germline-resource af-only-gnomad.hg38.vcf.gz \
--genotype-germline-sites \
--native-pair-hmm-threads 16 \
--panel-of-normals snv.pon.vcf.gz

tabix -f $bname.vcf.gz
```

## 4. Do segmentation with gatkcnv and cnv.pon for a tumor only sample

### 4.1 CollectReadCounts
```
# $1 = sort.bam
# $2 = interval.list - not gc annotated

bname=`basename $1 .bam`

gatk --java-options '-Xms500m -Xmx131859m -XX:+UseSerialGC -Djava.io.tmpdir=.' \
CollectReadCounts \
-I $1 \
-L $2 \
--interval-merging-rule OVERLAPPING_ONLY \
-O $bname.counts.hdf5
```

### 4.2 Denoise read counts
```
# $1 = input bam.counts.hdf5
# $2 = pon.hdf5
# $3 = panel.gcannotated.tsv - necessary for PureCN

bname=`basename $1 .counts.hdf5`

gatk --java-options "-Xmx12g" \
DenoiseReadCounts \
-I $1 \
--count-panel-of-normals $2 \
--standardized-copy-ratios $bname.standardizedCR.tsv \
--denoised-copy-ratios $bname.denoisedCR.tsv \
--annotated-intervals $3
```

### 4.3 Model segments
```
# $1 = T.denoisedCR.tsv

bname=`echo $1 | awk -F "." '{print $1}'`

gatk --java-options "-Xmx4g" ModelSegments \
--denoised-copy-ratios $1 \
--output . \
--output-prefix $bname
```

### 4.4 Call copy ratio
```
# #1 = T.cr.seg

bname=`basename $1 .cr.seg`

gatk CallCopyRatioSegments \
--input $1 \
--output $bname.called.seg
```

### 4.5 Plot modelled segments
```
# $1 = T.denoisedCR.tsv
# $2 hg38.dict

bname=`basename $1 .counts.denoisedCR.tsv`

gatk PlotModeledSegments \
--denoised-copy-ratios $1 \
--segments $bname.modelFinal.seg \
--sequence-dictionary $2 \
--minimum-contig-length 46709983 \
--output segment_plots \
--output-prefix $bname
```

## 5. Run pureCN with snv_pon and segments from gatkcnv

needs a bit of installation tweaking
https://bioconductor.org/packages/release/bioc/vignettes/PureCN/inst/doc/Quick.html#5_run_purecn_with_third-party_segmentation

### 5.1 Process intervals file
```
# $1 = panel.bed
# $2 = GCA_000001405.15_GRCh38_no_alt_analysis_set_100.bw

PURECN=/data/bcbio/anaconda/envs/r36/lib/R/library/PureCN/extdata

export PATH=/data/bcbio/anaconda/envs/r36/bin:$PATH

which Rscript

Rscript $PURECN/IntervalFile.R \
--infile $1 \
--fasta /data/genomes/Hsapiens/hg38/seq/hg38.fa \
--outfile intervals.txt \
--offtarget \
--genome hg38 \
--export baits_optimized_hg38.bed \
--mappability $2
```

### 5.2 Create normal DB
```
# $1 = snv_pon.vcf.gz from Mutect2 PON
# #2 = coverage.files.txt - list of target-coverage.hdf5
PURECN=/data/bcbio/anaconda/envs/r36/lib/R/library/PureCN/extdata

export PATH=/data/bcbio/anaconda/envs/r36/bin:$PATH

which Rscript

#Rscript $PURECN/NormalDB.R --help

Rscript $PURECN/NormalDB.R \
--outdir . \
--normal_panel $1 \
--assay exome_idt_v1 \
--genome hg38 \
--force \
--coveragefiles $2
```

### 5.3 Run pureCN
```
PURECN=/data/bcbio/anaconda/envs/r36/lib/R/library/PureCN/extdata

export PATH=/data/bcbio/anaconda/envs/r36/bin:$PATH

#which Rscript

# $1 = sample_id
# $2 = mapping_bias.rds
Rscript $PURECN/PureCN.R \
--sampleid $1 \
--out $1.out \
--tumor ${1}.counts.hdf5 \
--logratiofile $1.denoisedCR.tsv \
--segfile $1.modelFinal.seg \
--mappingbiasfile $2 \
--vcf $1.vcf.gz \
--statsfile $1.vcf.gz.stats \
--genome hg38 \
--funsegmentation Hclust \
--force \
--postoptimize \
--seed 123
```
