#!/usr/bin/perl

#sort align by id

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


#print @names."\t".length($dna)."\n";
foreach $name (sort @names)
{
    print $name."\n";
    print $align{$name}."\n";
}
close(IN);