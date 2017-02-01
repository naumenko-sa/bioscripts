# bioscripts - scripts I am using on a daily basis

# Fastq,bam,cram

* bam2fq.sh converts bam to fastq, uses samtools and bedtools
* bam.reads_number.sh reports the number of paired reads in a bam file, uses samtools
* bam.remove_region.sh removes reads from a bam file located at regions specified by a bed file.
Badly filtered rRNA-depleted RNA-seq samples may have huge coverage of low complexity regions. 
It is better to filter those with prinseq, however sometimes it is necessary to remove a particular region.
* cram2fq.sh converts cram to fastq, uses cramtools wrapper from bcbio

# RNA-seq

* rnaseq.star.sh - 2pass on the fly STAR alignment for a single sample
* [rnaseq.feature_counts.sh](../master/rnaseq.feature_counts.sh) [file.bam] calculates features (reads) for RPKM calculation in R, outputs length of the genes. TPMs are generally better
(http://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/), because you can compare expression values across samples, and they are calculated by default in bcbio rnaseq pipeline, 
but [http://www.gtexportal.org](GTEX) values unlike [http://www.proteinatlas.org/)(Protein atlas) values are in RPKMs, and still many people think in terms of RPKMs. GTEX has much more samples than HPA.

### Differential expression

* dexpression.katie.R - edgeR DE, batch effect correction, pheatmap, GO, pathways

### Splicing analysis

#### [qorts](http://hartleys.github.io/QoRTs/index.html)
* rnaseq.qorts.makeflatgff.sh
* rnaseq.qorts.qc.sh calculate counts from bam file, discovers junctions
* rnaseq.qorts.merge_novel_splices.sh merges junctions from all samples
* rnaseq.splicing.junction_seq.R runs junction seq analysis in R
* rnaseq.splicing.junction_seq.sh runs R script in the queue

# Variant analysis

* vcf.validate.sh - validate variant calls with Genome in a bottle callset using RTG vcfeval tool
* [VT: biallelic sites decomposition](https://github.com/atks/vt)
* [RTG: accurate vcf comparison](https://github.com/RealTimeGenomics/rtg-tools)

### Gemini staff - for CHEO, MH, and Muscular projects

* gemini.decompose.sh decomposes and normalizes variants with vt
* gemini.vep.sh annotates vcf file with VEP
* gemini.vep2gemini.sh loads VEP annotated vcf to the GEMINI database
* gemini.gemini2txt.sh dumps gemini database into txt file with decompressed genotypes
* gemini.gemini2report.R creates nice report for import to excel from the gemini.txt dump
* gemini.from_rnaseq.sh creates gemini database and rare harmful variants report from bcbio's rna-seq pipeline output

# By project

## CHEO is a study of families with rare genetic conditions - [Care For Rare](http://care4rare.ca/) at [Children's Hospital of Eastern Ontario](http://www.cheori.org/)
For variant calling I use [bcbio ensemble approach](https://bcbio-nextgen.readthedocs.io/en/latest/contents/configuration.html#ensemble-variant-calling)
on per-family basis.  In brief, 2 out of 4 (gatk-haplotype, samtools, freebayes, and platypus) algorithms should be voting for a variant to be called.
This allows to achieve increased sensitivity required for research, compared to conservative strategy of the genetic testing laboratory.

1. [cheo.check_if_done.sh](../master/cheo.check_if_done.sh) [bcbio_job.output] check which bcbio jobs are done, useful when running 100x families
2. [cheo.postprocess.sh](../master/cheo.postprocess.sh) [family] cleans up after bcbio and prepares necessary tables for excel report generator
3. [gemini.gemini2report.R](../master/gemini.gemini2report.R) generates reports for excel import 
4. [cheo.c4r_database.sh](../master/cheo.c4r_database.sh) prepares variants from a family to be merged in a database seen_in_c4r
5. [cheo.c4r_database_merge.pl](../master/cheo.c4r_database_merge.pl) merges variant evidence from many samples

---

## Muscle is a study of genetic causes of muscular diseases