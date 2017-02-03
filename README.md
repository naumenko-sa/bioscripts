# bioscripts - scripts I am using on a daily basis

[Fastq,Bam,Cram](https://github.com/naumenko-sa/bioscripts#fastqbamcram) [Genome assembly](https://github.com/naumenko-sa/bioscripts#genome-assembly)
[RNA-seq](https://github.com/naumenko-sa/bioscripts#rna-seq) [Variant analysis](https://github.com/naumenko-sa/bioscripts#variant-analysis)

[Ordered by project](https://github.com/naumenko-sa/bioscripts#by-project)

# Fastq,bam,cram

* [bam2fq.sh](../master/bam2fq.sh) converts bam to fastq, uses samtools and bedtools
* [bam.reads_number.sh](../master/bam.reads_number.sh) reports the number of paired reads in a bam file, uses samtools
* [bam.remove_region.sh](../master/bam.remove_region.sh) removes reads from a bam file located at regions specified by a bed file.
Badly filtered rRNA-depleted RNA-seq samples may have huge coverage of low complexity regions. 
It is better to filter those with [prinseq](http://http://prinseq.sourceforge.net/), however sometimes it is necessary to remove a particular region.
* [cram2fq.sh](../master/cram2fq.sh) converts cram to fastq, uses cramtools wrapper from bcbio

# Genome assembly
The wisdom here is to avoid large genomes, polyploid genomes, and creating your own genome assembly. 
See my [lecture](http://makarich.fbb.msu.ru/snaumenko/ngs_lecture/naumenko.genome_assembly-n.pdf) (in Russian).
For large genomes it is better to have multiple libraries, with substantial amoung of mate pairs with 5k,10k,20k insert size. 
For a serious work a special computing node is necessary (1-2T RAM). Surprisingly, such a node is not that expensive: just buy
a cheap 4CPU SuperMicro server capable to carry up to 4-8T RAM, buy RAM, and insert it into server. Avoid vendors and sales persons.
Look for engeneers to cooperate.


* [genome_assembly.spades.pbs](../master/genome_assembly.spades.pbs) runs [spades](http://bioinf.spbau.ru/spades) assembler. Spades is the best assembler
for genomes up to 100G.

# Phylogenetics

# RNA-seq

* [rnaseq.star.sh](../master/rnaseq.star.sh) - 2pass on the fly STAR alignment for a single sample. Two passes are recommended to enhance alignment and calculation of counts for novel splice junctions
(1st pass discovers junctions, 2nd pass makes an alignment).
* [rnaseq.feature_counts.sh](../master/rnaseq.feature_counts.sh) [file.bam] calculates features (reads) for RPKM calculation in R, outputs length of the genes. TPMs are generally better
(http://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/), because you can compare expression values across samples, and they are calculated by default in bcbio rnaseq pipeline, 
but [GTEX](http://www.gtexportal.org) values unlike [Protein Atlas](http://www.proteinatlas.org/) values are in RPKMs, and still many people think in terms of RPKMs. GTEX has much more samples than HPA.
* [rnaseq.load_rpkm_counts.R](../master/rnaseq.load_rpkm_counts.R) calculates RPKM counts from feature_counts output.

### De noto transcriptome assembly

### Differential expression

* dexpression.katie.R - edgeR DE, batch effect correction, pheatmap, GO, pathways

### Micro-RNA
I did a couple of analysis in *A.thaliana* and *S.tuberosum* (potato). First I run [smallRNA-seq pipeline from bcbio](https://bcbio-nextgen.readthedocs.io/en/latest/contents/pipelines.html#smallrna-seq).
It is much better tuned for human miRNA, where it does discovery of novel miRNA (because primary tools are developed mostly for human).
For non-human species it does a good job of mapping to mirBase and has an extensive [report](https://github.com/lpantano/seqcluster/blob/master/seqcluster/templates/report.rmd) in R.  
For novel miRNA discovery I tried [shortstack](http://sites.psu.edu/axtell/software/shortstack/).

### Splicing analysis

[qorts](http://hartleys.github.io/QoRTs/index.html) and [junctionseq](https://bioconductor.org/packages/release/bioc/html/JunctionSeq.html) do the job.
Finally they produce html report for all differentially expressed genes, plotting exon usage, junction counts, novel juctions, as a nice plots.
The problem is that RNA is quite noisy by its biological nature, so expect to have many events and many false positives.
Additionaly you have tracks for IGV. They are useful, because default exon usage plots are not for exons from the canonical isoform,
but for a set of chunks from flattened annotation of many isoforms, and people are asking questions, why here are 125 exons not 108. It is hard to see
the correspondence to exons of the canonical isoforms for long genes from those plots.

Additionaly I use sailfish relative isoform expression levels which [bcbio rna-seq pipeline](https://bcbio-nextgen.readthedocs.io/en/latest/contents/pipelines.html#rna-seq) 
outputs by default.

#### [qorts](http://hartleys.github.io/QoRTs/index.html)
* rnaseq.qorts.makeflatgff.sh
* [rnaseq.qorts.qc.sh](../master/rnaseq.qorts.qc.sh) [file.bam] calculates counts from bam file
* rnaseq.qorts.merge_novel_splices.sh merges junctions from all samples
* [rnaseq.qorts.get_novel_exons.sh](../master/rnaseq.qorts.get_novel_exons.sh) [sample] [max_novel_exon_lengths] 
discovers novel exons using qorts files as input. A novel exon is reported when two novel junctions are found at the distance less then max_novel_exon_lengths.    
* rnaseq.splicing.junction_seq.R runs junction seq analysis in R
* rnaseq.splicing.junction_seq.sh runs R script in the queue

# Variant analysis, vcf files

* vcf.validate.sh - validate variant calls with Genome in a bottle callset using RTG vcfeval tool
* [VT: biallelic sites decomposition](https://github.com/atks/vt)
* [RTG: accurate vcf comparison](https://github.com/RealTimeGenomics/rtg-tools)

### Gemini staff - for CHEO, MH, and Muscular projects

[Gemini](https://gemini.readthedocs.io/en/latest/) is a database and framework for variant analysis. [Bcbio](http://bcbio-nextgen.readthedocs.io/en/latest/)
variant calling pipeline outputs variants in gemini format.

* [gemini.decompose.sh](../master/gemini.decompose.sh) decomposes and normalizes variants with vt
* [gemini.vep.sh](../master/gemini.vep.sh) annotates vcf file with VEP
* [gemini.vep2gemini.sh](../master/gemini.vep2gemini.sh) loads VEP annotated vcf to the GEMINI database
* [gemini.gemini2txt.sh](../master/gemini.gemini2txt.sh) dumps gemini database into txt file with decompressed genotypes
* [gemini.gemini2report.R](../master/gemini.gemini2report.R) creates nice report for import to excel from the gemini.txt dump
* [gemini.from_rnaseq.sh](../master/gemini.from_rnaseq.sh) creates gemini database and rare harmful variants report from bcbio's rna-seq pipeline output

# By project

## 6. Muscle [2016-] is a study of genetic causes of muscular diseases

---

## 5. MH [2016-] studies [Malignant hyperthermia](https://en.wikipedia.org/wiki/Malignant_hyperthermia)
1. [project_mh.R](../master/project_mh.R)
2. [rnaseq.qorts.get_novel_exons.sh](../master/rnaseq.qorts.get_novel_exons.sh)
3. [project_mh.RYR1.isoforms.txt](../master/project_mh.RYR1.isoforms.txt)

## 4. CHEO [2016-] is a study of families with rare genetic conditions - [Care For Rare](http://care4rare.ca/) at [Children's Hospital of Eastern Ontario](http://www.cheori.org/)
For variant calling I use [bcbio ensemble approach](https://bcbio-nextgen.readthedocs.io/en/latest/contents/configuration.html#ensemble-variant-calling)
on per-family basis.  In brief, 2 out of 4 (gatk-haplotype, samtools, freebayes, and platypus) algorithms should be voting for a variant to be called.
This allows to achieve increased sensitivity required for research, compared to conservative strategy of the genetic testing laboratory.

1. [cheo.check_if_done.sh](../master/cheo.check_if_done.sh) [bcbio_job.output] check which bcbio jobs are done, useful when running 100x families
2. [cheo.postprocess.sh](../master/cheo.postprocess.sh) [family] cleans up after bcbio and prepares necessary tables for excel report generator
3. [gemini.gemini2report.R](../master/gemini.gemini2report.R) generates reports for excel import 
4. [cheo.c4r_database.sh](../master/cheo.c4r_database.sh) prepares variants from a family to be merged in a database seen_in_c4r
5. [cheo.c4r_database_merge.pl](../master/cheo.c4r_database_merge.pl) merges variant evidence from many samples

---
## 3. Gammaruses [2012-]
---

## 2. Spectrum [2010-2012]
---
## 1. Reversals [2009-2012]