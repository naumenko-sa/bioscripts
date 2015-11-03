#!/usr/bin/perl

#intput fna + qual files
#outpu: fastq

use strict;

open (DNA,$ARGV[0]);
open (QUAL,$ARGV[1]);

my $name;

while(<DNA>)
{
    my $seq;
    my $qual;
    if ($_ =~ /^>/)
    {
	$name = $_;
	my $t = <QUAL>;
    }
    else
    {
	$seq = $_;
	$qual = <QUAL>;
	chomp $qual;
    
	print $name;
	print $seq;
	print "+\n";
	my @ar = split(' ',$qual);
	my $qual_str;
	foreach my $m (@ar)
	{
    	    $qual_str .= chr($m+64);
	}
	print $qual_str."\n";
    }
}

close(DNA);
close(QUAL);