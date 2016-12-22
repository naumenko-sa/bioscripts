# bioscripts - scripts I am using on a daily basis

# RNA-seq

## rnaseq.star.sh - 2pass on the fly STAR alignment for a single sample

### Differential expression

* dexpression.katie.R - edgeR DE, batch effect correction, pheatmap, GO, pathways

### Splicing analysis

#### [qorts](http://hartleys.github.io/QoRTs/index.html)
* rnaseq.qorts.makeflatgff.sh
* rnaseq.qorts.qc.sh - calculate counts from bam file

# Variant analysis

* vcf.validate.sh - validate variant calls with Genome in a bottle callset using RTG vcfeval tool
* [VT: biallelic sites decomposition](https://github.com/atks/vt)
* [RTG: accurate vcf comparison](https://github.com/RealTimeGenomics/rtg-tools)

### Gemini staff - for CHEO, MH, and Muscular projects

* gemini.decompose.sh - decompose and normalize variants with vt
* gemini.vep.sh - annotate vcf file with VEP
* gemini.vep2gemini.sh - load VEP annotated vcf to the GEMINI database
* gemini.gemini2txt.sh - dump a gemini database into txt file with decompressed genotypes
* gemini.gemini2report.R  - create a nice report for import to excel from the gemini.txt dump

