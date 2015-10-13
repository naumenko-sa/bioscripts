#!/usr/bin/perl
#blast outpus many hits in the case num_alignments = 1
%genes;

open (IN, $ARGV[0]);
open (OUT, ">".$ARGV[0].".fixed");

while (<IN>)
{
    chomp;
    
    @ar = split (' ',$_);
    #print $ar[0]."\n";

    if (!exists ($genes{$ar[0]}))
    {
	print OUT $_."\n";
	$genes{$ar[0]}=1;
    }
}

close(IN);
close(OUT);