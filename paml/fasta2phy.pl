#!/usr/bin/perl

#converts fasta to relaxed phylip format (without 10 base limit for name)

open(IN, $ARGV[0]);

%align;
@names;

while(<IN>)
{
    chomp;
    $name = $_;
    push(@names,$name);
    $dna = <IN>;
    chomp $dna;
    $align{$name}=$dna;
}


print @names."\t".length($dna)."\n";
foreach $name (@names)
{
    print substr($name,1,length($name)-1)."  ".$align{$name}."\n";
}
close(IN);
