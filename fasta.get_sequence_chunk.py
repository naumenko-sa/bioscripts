#!/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/anaconda/bin/python

from Bio import SeqIO
for record in SeqIO.parse("mm10.fa","fasta"):
    if(record.id == 'chrX'):
        print('>'+record.id)
	print(record.seq[83729320:83927998])
