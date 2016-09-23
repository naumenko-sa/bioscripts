#!/usr/bin/perl

#adds annotation for snps with ids in format : chr:pos-alt-ref to gemini table by ENSG id

open(VCF,$ARGV[1]);
open(GEMINI,$ARGV[0]);

%scores1;

while(<VCF>)
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
    $arr[0] =~ s/chr//;
    $id=$arr[0]."-".$arr[2]."-".$arr[4]."-".$arr[5];
    print $_."\t".$scores1{$id}."\n";
}

close(OMIM);
close(GEMINI);