#!/usr/bin/perl
use strict;

# remove IUPAC ambiguity codes with random letter
# http://www.boekhoff.info/?pid=data&dat=fasta-codes

my $id; my $seq;

open (INFILE, $ARGV[0]);

while (<INFILE>) {
    chomp $_;
    if (/>/) 
    {
	print $_."\n";
    }
    else 
    {
	$seq = $_;
	$seq =~ tr/KMRYSWBVHD/GAATCTCATG/;
	print $seq."\n";
    }
}
