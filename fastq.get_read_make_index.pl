#!/usr/bin/perl
use strict;
use Bio::Index::Fastq;

#index fastq file
#usage get_read_index_fastq.pl [idx_name] [files.fastq]

my $Index_File_Name = shift;
my $inx = Bio::Index::Fastq->new('-filename' => $Index_File_Name,'-write_flag' => 1);
$inx->make_index(@ARGV);