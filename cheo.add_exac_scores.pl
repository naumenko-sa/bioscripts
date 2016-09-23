#!/usr/bin/perl

#adds the exac scores to the table from gemini

open(EXAC,"/home/naumenko/work/cheo.variants/data/BED/exac.scores.sorted");
open(GEMINI,$ARGV[0]);

%scores1;
%scores2;

while(<EXAC>)
{
    chomp;
    @arr = split("\t", $_);
    $scores1{$arr[0]}=$arr[2];
    $scores2{$arr[0]}=$arr[3];
}


while (<GEMINI>)
{
    chomp;
    @arr = split ("\t",$_);
    $transcript=$arr[11];
    print $_."\t".$scores1{$transcript}."\t",$scores2{$transcript}."\n";
}

close(EXAC);
close(GEMINI);