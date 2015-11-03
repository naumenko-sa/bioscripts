#!/usr/bin/perl

open(IN,$ARGV[0]);

@aln;

while(<IN>)
{
    chomp;
    $dna = <IN>;
    chomp $dna;
    push (@aln,$dna);
}

close(IN);

for ($i=0;$i<length($aln[0]);$i++)
{
    my %spectrum;
    for ($j=0;$j<@aln;$j++)
    {
	$nuc = substr($aln[$j],$i,1);
	$spectrum{$nuc}++;
    }
    $k=$i+1;
    print $k."\t".keys(%spectrum)."\n";
}