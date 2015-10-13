#!/usr/bin/perl

%freq;

while(<STDIN>)
{
    chomp;
    @ar = split(' ',$_);
    $freq{$ar[0]}+=$ar[1];
}

$total=0;

foreach ( sort { $freq{$b} <=> $freq{$a}} keys %freq )
{
    print $_."\t".$freq{$_}."\n";
}