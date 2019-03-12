|By subject |[By project](https://github.com/naumenko-sa/bioscripts#by-project)|
|---|---|
|[Alignment](https://github.com/naumenko-sa/bioscripts#alignment)||
|[Articles and resources](https://github.com/naumenko-sa/bioscripts#articles-and-resources)||
|[Annotation](https://github.com/naumenko-sa/bioscripts#annotation)||
|[bam, cram, bed, fastq](https://github.com/naumenko-sa/bioscripts#bam-cram-bed-fastq)||
|[Coverage](https://github.com/naumenko-sa/bioscripts#coverage)||
|[Genome assembly](https://github.com/naumenko-sa/bioscripts#genome-assembly)||
|[Phylogenetics](https://github.com/naumenko-sa/bioscripts#phylogenetics)||
|[MISC](https://github.com/naumenko-sa/bioscripts#misc)||
|[RNA-seq](https://github.com/naumenko-sa/bioscripts#rna-seq)||
|[Variant analysis](https://github.com/naumenko-sa/bioscripts#variant-analysis-vcf-files)||
|[Visualization](https://github.com/naumenko-sa/bioscripts#visualization)||

# Alignment
* [Jim Kent's utils UCSC](http://hgdownload.soe.ucsc.edu/admin/exe/)
* [Synteny maps](http://www.agcol.arizona.edu/software/symap/)
* alignment.avr_pw_dist.pl [alignment.fasta] [1=print distances] prints average pairwise distance for multiple alignment, and all pairwise distances if $2==1 
* alignment.check_3x.sh [alignment.fasta] prints how many sequences in an alignment are not 3x in length
* alignment.check_conservative_site.sh
* alignment.consensus.2seq.pl [alignment.fasta] prints a consensus of two DNA sequences
* alignment.consensus.pl
* alignment.count_gaps.sh counts % of positions with gaps in an alignment
* [alignment.detect_stops.pl](../master/alignment.detect_stops.pl) [alignment.fasta] [--print_pos] prints internal stops and their positions
* alignment.fa2fasta.pl converts fixed line lengths fasta to long line format.
* alignment.fasta2slice.pl
* alignment.filter_not3x.sh
* alignment.fix_ends.sh
* alignment.muscle.sh.
* alignment.nucleotide_diversity.pl
* alignment.occupancy.pl
* [alignment.protein2dna.sh](../master/alignment.protein2dna.sh) [name.aln.fasta] reverse translates amino acid alignment into DNA alignment.
* alignment.remove_ambiguity.pl
* alignment.remove_cons_n_gaps.pl
* alignment.remove_gaps.pl
* alignment.remove_seqs_having_stops.pl
* alignment.remove_sp4occupancy.pl
* [alignment.remove_stops_n_gaps2.pl](../master/alignment.remove_stops_n_gaps2.pl) removes stops, gaps, and sites around indels.
* [alignment.remove_stops_n_gaps.pl](../master/alignment.remove_stops_n_gaps.pl) removes triplets with stops or gaps.
* alignment.reverse_complement.pl
* alignments.concatenate1.pl
* alignments.concatenate.pl
* alignment.sort_by_id.pl
* alignment.substitution_profile.pl
* alignment.tr_trv.pl
* [alignment.UCSC.dm6droSim1.sh](../alignment.UCSC.dm6droSim1.sh) builds a pairwise alignment of D.melanogaster and D.simulans.

# Annotation
* [genes.R](../master/genes.R) - various ENSEMBL annotations from biomaRt

# bam, cram, bed, fastq
* [jvarkit](https://github.com/lindenb/jvarkit/wiki): bamstats04,bamstats05.
* [bedtools](http://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html) - genome arithmetics.
* [fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic)
* [PRINSEQ](http://prinseq.sourceforge.net/)
* [bam2fq.sh](../master/bam2fq.sh) converts bam to fastq, uses samtools and bedtools
* For coverage first try cre.coverage.sh in cre.
    * [bam.coverage.sh](../master/bam.coverage.sh) outputs reads mapped to the interval specified by a bed file.
    * [bam.coverage.bamstats](../master/bam.coverage.bamstats.sh) prints average coverage for a gene, uses bamstats05.
    * [bam.gene_coverage.sh](../master/bam.gene_coverage.sh) prints mean coverage per region (exon,gene).
* [bam.reads_number.sh](../master/bam.reads_number.sh) reports the number of paired reads in a bam file, uses samtools
* [bam.remove_region.sh](../master/bam.remove_region.sh) removes reads from a bam file located at regions specified by a bed file.
Badly filtered rRNA-depleted RNA-seq samples may have huge coverage of low complexity regions. 
It is better to filter those with [prinseq](http://http://prinseq.sourceforge.net/), however sometimes it is necessary to remove a particular region.
* [bam.sort.sh](../master/bam.sort.sh) - sorts a bam file before calculating coverage.
* [basespace-cli](https://help.basespace.illumina.com/articles/descriptive/basespace-cli/). New Illumina sequencers upload data into the basespace cloud.
bs utility allows to access this data from HPC. To copy bcl files: `bs cp //./Runs/<project_name>/Data .`
* [bcl2fastq.sh](../master/bcl2fastq.sh) converts raw bcl Illumina files to fastq.
* [cram2fq.sh](../master/cram2fq.sh) converts cram to fastq, uses cramtools wrapper from bcbio
* `samtools quickcheck -vvvv [file.bam]` checks the integrity of a bam file.

# Coverage
* bam.coverage.bamstats05.sh - very quick coverage calculation for a bam and a bed file with bamstats05
* bam.coverage.sh - slow base-wise coverage calculation with bedtools + median coverage statistics

# Genome assembly
The wisdom here is to avoid large genomes, polyploid genomes, and creating your own genome assembler.
See my [lecture](http://makarich.fbb.msu.ru/snaumenko/ngs_lecture/naumenko.genome_assembly-n.pdf) (in Russian).
For large genomes it is better to have multiple libraries, with the substantial amount of mate pairs with 5k,10k,20k insert size. 
For a serious work a special computing node is necessary (1-2T RAM). Surprisingly, such a node is not that expensive: just buy
a cheap 4CPU SuperMicro server capable to carry up to 4-8T RAM, buy RAM, and insert it into server. Avoid vendors and sales persons.
Look for engineers to cooperate. For smaller genomes I prefer spades, for larger ones velvet + platanus. Remember to clean up reads
(check for contamination, quality trimming).

* [genome_assembly.spades.pbs](../master/genome_assembly.spades.pbs) runs [spades](http://bioinf.spbau.ru/spades) assembler. Spades is the best assembler
for genomes up to 100G.
* [genome_assembly.n50.sh](../master/genome_assembly.n50.sh) [contigs.fasta] calculates N50 metrics.

# Phylogenetics
Nothing is comparable to the feeling when you just have plotted a phylogenetic tree. It is very rewarding however the tree might be misleading.
Read [The Phylogenetic handbook](https://books.google.ca/books/about/The_Phylogenetic_Handbook.html?id=DeD_lQ-kBPQC&redir_esc=y).
In brief, it is necessary to build a good alignment(!), concatenate many genes, fit the model with modeltest (GTR+Г+I is usuallly the winner),
and run RAXML and MrBayes to compare two trees. I visualize trees with [Dendroscope](http://dendroscope.org/). It becomes better as years pass.
* [tree.raxml_boot.pbs](../master/tree.raxml_boot.pbs) infers a phylogenetic tree for cdna alignment with 100 bootstrap replicates.

# MISC

* [hla.athlates.sh](../master/hla.athlates.sh) runs Athlates program for HLA allele typing.

# RNA-seq
* [bcbio rna-seq pipeline](https://bcbio-nextgen.readthedocs.io/en/latest/contents/pipelines.html#rna-seq) does STAR alignment, variant calling,
expression measurements, isoform reconstruction, and quality metrics. Run sample-wise to speed up.
* [rnaseq.star.sh](../master/rnaseq.star.sh) - 2pass on the fly STAR alignment for a single sample. Two passes are recommended to enhance alignment and calculation of counts for novel splice junctions
(1st pass discovers junctions, 2nd pass makes an alignment).
* [rnaseq.feature_counts.sh](../master/rnaseq.feature_counts.sh) [file.bam] calculates features (reads) for RPKM calculation in R, outputs length of the genes. TPMs are generally better
(http://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/), because you can compare expression values across samples, and they are calculated by default in bcbio rnaseq pipeline, 
but [GTEX](http://www.gtexportal.org) values unlike [Protein Atlas](http://www.proteinatlas.org/) values are in RPKMs, and still many people think in terms of RPKMs. GTEX has much more samples than HPA.
* [rnaseq.load_rpkm_counts.R](../master/rnaseq.load_rpkm_counts.R) calculates RPKM counts from feature_counts output.

### De novo transcriptome assembly

### Differential expression (DE)
* [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html).
* [rnaseq.dexpression.R](../master/rnaseq.dexpression.R) - common utils for DE.
* [rnaseq.dexpression.katie.R](../master/rnaseq.dexpression.katie.R) - edgeR DE, batch effect correction, pheatmap, GO, pathways.
* [Gene Set Enrichment Analysis (GSEA)](http://software.broadinstitute.org/gsea/index.jsp). Run with all protein-coding genes, expression in cpm.

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
but for a set of chunks from the flattened annotation of many isoforms, and people are asking questions, why here are 125 exons not 108. It is hard to see
the correspondence to exons of the canonical isoforms for long genes from those plots. Loading tracks in IGV helps: you see junctions, counts, and the alignment,
people begin to understand you.

Additionaly I use sailfish relative isoform expression levels which bcbio rna-seq pipeline outputs by default.

#### [qorts](http://hartleys.github.io/QoRTs/index.html)
* rnaseq.qorts.makeflatgff.sh
* [rnaseq.qorts.qc.sh](../master/rnaseq.qorts.qc.sh) [file.bam] calculates counts from bam file
* rnaseq.qorts.merge_novel_splices.sh merges junctions from all samples
* [rnaseq.qorts.get_novel_exons.sh](../master/rnaseq.qorts.get_novel_exons.sh) [sample] [max_novel_exon_lengths] 
discovers novel exons using qorts files as input. A novel exon is reported when two novel junctions are found at the distance less then max_novel_exon_lengths.    
* rnaseq.splicing.junction_seq.R runs junction seq analysis in R
* rnaseq.splicing.junction_seq.sh runs R script in the queue

# Variant analysis, vcf files
For variant calling I use [bcbio ensemble approach](https://bcbio-nextgen.readthedocs.io/en/latest/contents/configuration.html#ensemble-variant-calling)
on per-family basis.  In brief, 2 out of 4 (gatk-haplotype, samtools, freebayes, and platypus) algorithms should be voting for a variant to be called.
This allows to achieve increased sensitivity required for research, compared to conservative strategy of the genetic testing laboratory.

* [VT toolkit](http://genome.sph.umich.edu/wiki/Vt),[VT toolkit source](https://github.com/atks/vt):biallelic sites decomposition, info fields removal
* [bcftools](http://www.htslib.org/doc/bcftools.html)
* [RTG: accurate vcf comparison](https://github.com/RealTimeGenomics/rtg-tools), vcfstats.
* [vcf.rtg.validate.sh](../master/vcf.rtg.validate.sh) - how to run rtg vcfeval validation.
* [bcbio.pbs](../master/bcbio.pbs) [project] submits bcbio project to the queue.
* [bcbio.cleanup.sh](../master/bcbio.cleanup.sh) [family] cleans up after bcbio and prepares necessary tables for excel report generator.
* [bcbio.prepare_families.sh](../master/bcbio.prepare_families.sh) creates symlinks, folders, config files to run variant calling for multiple families.
* [bcbio.prepare_samples.sh](../master/bcbio.prepare_samples.sh) creates symlinks, folder structure, config files to run multiple projects for bam file generation.
* [vcf.split_multisample.sh](../master/vcf.split_multisample.sh) splits multisample vcfs.
* [vcf.kinship.R](../master/vcf.kinship.R) calculates kinship for a family using SNPRelate and plots a pedigree.
Sometimes it is useful when studying cohorts of families. When samples are mislabeled (it happens), the whole study makes no sense,
at least until relatives are placed together to call and interpret variants.

### Gemini staff - for CHEO, MH, and Muscular projects

[Gemini](https://gemini.readthedocs.io/en/latest/) is a database and framework for variant analysis. [Bcbio](http://bcbio-nextgen.readthedocs.io/en/latest/)
variant calling pipeline outputs variants in gemini format.

* [gemini.decompose.sh](../master/gemini.decompose.sh) decomposes and normalizes variants with vt.
* [gemini.vep.sh](../master/gemini.vep.sh) annotates vcf file with VEP.
* [gemini.vep2gemini.sh](../master/gemini.vep2gemini.sh) loads VEP annotated vcf to the GEMINI database.
* [gemini.gemini2txt.sh](https://github.com/naumenko-sa/cre/blob/master/gemini.gemini2txt.sh) dumps gemini database into txt file with decompressed genotypes.
* [cre.R](https://github.com/naumenko-sa/cre/blob/master/cre.R) creates report for import to excel from the gemini.txt dump.
* [gemini.from_rnaseq.sh](../master/gemini.from_rnaseq.sh) creates gemini database and rare harmful variants report from bcbio's rna-seq pipeline output.
* [gemini.refseq.sh](https://github.com/naumenko-sa/cre/blob/master/gemini.refseq.sh) from vcf to refseq table for excel report generator.
* [gemini.vep.parse.pl](../master/gemini.vep.parse.pl) parses VEP annotation.
* [gemini.vep.refseq.sh](../master/gemini.vep.refseq.sh) annotates variants with RefSeq transcripts (NMs) using VEP.
* [gemini.variant_impacts.sh](../master/gemini.variant_impacts.sh) extracts variant impacts for all transcripts.

# Articles and resources
* [https://precision.fda.gov](https://precision.fda.gov/) is an effort to standardize computational pipelines.
* [https://www.clinicalgenome.org/](https://www.clinicalgenome.org/) is created by authors of ACMG guidelines.
* [Viral Genomics & Bioinformatics - Glasgow](http://bioinformatics.cvr.ac.uk/)
* [Gene set enrichment analysis](http://software.broadinstitute.org/gsea/index.jsp)
* [Pedigree chart designer](http://www.cegat.de/en/for-physicians/pedigree-chart-designer/) draws pedigree diagrams.
* [List of journals for the journal/article club](../master/journals.md)

# By project
## 6. RNA-seq in Mendelian muscle diseases [2016-2018]
* Developed into [crt](https://github.com/naumenko-sa/crt)
* Published in [AJHG](https://www.sciencedirect.com/science/article/pii/S0002929719300126)

## 5. [Malignant hyperthermia](https://en.wikipedia.org/wiki/Malignant_hyperthermia) [2016]
* Continued by [Alper Celik](https://github.com/celalp)
* [crt.malignant_hyperthermia.R](https://github.com/naumenko-sa/crt/blob/master/malignant_hyperthermia.R)
* [crt.qorts.get_novel_exons.sh](https://github.com/naumenko-sa/crt/blob/master/rnaseq.qorts.get_novel_exons.sh)

## 4. CHEO [2016-2019]
Developed into cre
1. [bcbio.prepare_families.sh](../master/bcbio.prepare_families.sh) creates symlinks, folders, config files to run variant calling for multiple families,
or use bcbio.prepare_samples.sh if some samples need to re-generate bam files.
2. [bcbio.array.pbs](../master/bcbio.array.pbs) runs multiple bcbio projects as a job array. I have tried a parallel execution of BCBIO, with IPython,
quite successfully, but finally with enigmatic faults,maybe because of the cluster issues.
3. [cheo.check_if_done.sh](../master/cheo.check_if_done.sh) [bcbio_job.output] check which bcbio jobs are done, useful when running 100x families
4. [bcbio.cleanup.sh](../master/bcbio.cleanup.sh) [family] cleans up after bcbio and prepares necessary tables for excel report generator
5. [gemini.gemini2report.R](../master/gemini.gemini2report.R) generates reports for excel import 
6. [cheo.c4r_database.sh](../master/cheo.c4r_database.sh) prepares variants from a family to be merged in a database seen_in_c4r
7. [cheo.c4r_database_merge.pl](../master/cheo.c4r_database_merge.pl) merges variant evidence from many samples
8. gemini.gemini2report.R again, to add information about C4R frequencies.
9. [cheo.filtered_vcf.sh](../master/cheo.filtered_vcf.sh) prepares filtered vcf file corresponding to excel report.

## 3. Gammarus [2012-2016]
* [gam.filter.sh](../master/gam.filter.sh) - filters reads
* gam.blastp.pbs, gam.blastp.sh - blastp for othomcl
* [gam.ortho.sh](../master/gam.ortho.sh) - orthologization of transcript with orthomcl
* [gam.orthomcl.make_alignments.pl](../master/gam.orthomcl.make_alignments.pl) - writes fasta alignments for N species (1:1 orthologs) after orthomcl

## 2. Spectrum [2010-2012]

## 1. Reversals [2009-2012]
