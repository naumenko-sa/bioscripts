#!/usr/bin/perl

#adds omim string to gene to gemini table by ENSG id

open(OMIM,"/home/naumenko/work/cheo.variants/data/BED/omim.orphanet/omim/omim.forannotation1");
open(GEMINI,$ARGV[0]);

%scores1;

while(<OMIM>)
{
    chomp;
    @arr = split("\t", $_);
    $scores1{$arr[0]}=$arr[1];
    #print $arr[0]."\t".$arr[1]."\n";
}


while (<GEMINI>)
{
    chomp;
    @arr = split ("\t",$_);
    $ensembl=$arr[10];
    print $_."\t\"".$scores1{$ensembl}."\"\n";
}

close(OMIM);
close(GEMINI);