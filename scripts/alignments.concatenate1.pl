#!/usr/bin/perl

#concatenate alignments

@files=<*.fasta>;
%alignment;

foreach $file (@files)
{
    open(IN,$file);
    while(<IN>)
    {
	chomp;
	$name=$_;
	$dna = <IN>;
	chomp $dna;
	if (exists($alignment{$name}))
	{
	    $alignment{$name}.=$dna;
	}
	else
	{
	    $alignment{$name}=$dna;
	}
    }
    close (IN);
}
foreach $key (keys %alignment)
{
    print $key."\n";
    print $alignment{$key}."\n";
}