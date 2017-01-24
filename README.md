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

## CHEO

* cheo.check_if_done.sh [bcbio_job.output]