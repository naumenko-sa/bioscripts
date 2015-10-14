#!/usr/bin/perl
use strict;
use Bio::Index::Fastq;

#gets a read with quality values from fastq file

#usage get_read.pl [read_name] [reads.database.fastq.index]
#works with index created by bp_index.pl -dir . [file.fasta.idx - new file to create] [file-fasta-to index]

# Print out several sequences present in the index
# in Fastq format
my $Index_File_Name = shift;
my $inx = Bio::Index::Fastq->new('-filename' => $Index_File_Name);
my $out = Bio::SeqIO->new('-format' => 'Fastq','-fh' => \*STDOUT);
                        
foreach my $id (@ARGV) 
{
        my $seq = $inx->fetch($id); # Returns Bio::Seq::Quality object
	$out->write_seq($seq);
}