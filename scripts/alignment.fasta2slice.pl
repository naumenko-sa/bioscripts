#!/usr/bin/perl

#get sites
#prints 0-based site number (aa number) and all alleles
@aln;

open(IN,$ARGV[0]);
while(<IN>)
{
    chomp;
    $name=$_;
    $dna = <IN>;
    chomp $dna;
    push (@aln,$dna);
}

for ($i=0;$i<length($aln[0]);$i++)
{
    $site='';
    for ($j=0;$j<@aln;$j++)
    {
	$site.=substr($aln[$j],$i,1);
    }
    print $i."\t".$site."\n";
}