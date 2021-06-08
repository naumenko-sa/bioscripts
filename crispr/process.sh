
# convert docker to singularity
singularity pull docker://pinellolab/crispresso2
# test
singularity run -e crispresso2_latest.sif CRISPResso -h

singularity run -e crispresso2_latest.sif CRISPResso --version

# otherwise r1 and r2 are not even
bamtofastq filename=sample-ready.bam F=sample_1.fq F2=sample_2.fq O=unmatched_1.fq O2=unmatched_2.fq S=single.fq

#crispresso needs small baits
cat baits.mm10.sorted.bed | awk '{print $1"\t"$2"\t"$3-60"\n"$1"\t"$3-60"\t"$3}' > baits.mm10.sorted.split.bed

bedtools getfasta -fi /path/to/reference.fa -bed baits.mm10.sorted.split.bed > baits.mm10.sorted.split.fasta

cat baits.mm10.sorted.split.fasta | awk '{name=$1;getline;print name"\t"$0"\tNA"}' | sed s/">"// > baits.mm10.sorted.split.tsv

singularity run -e crispresso2_latest.sif CRISPRessoPooled -r1 sample_1.fq -r2 sample_2.fq -f baits.mm10.sorted.split.tsv -p 10