#!/usr/bin/perl

%freq;

while(<STDIN>)
{
    chomp;
    $freq{$_}++;
}

$total=0;
foreach $key (keys %freq)
{
    $total += $freq{$key};    
}

foreach ( sort { $freq{$b} <=> $freq{$a}} keys %freq )
{
    print $_."\t".$freq{$_}."\n";
    #substr($freq{$_}/$total,0,5)."\n";
}