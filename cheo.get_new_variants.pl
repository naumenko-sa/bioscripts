#!/usr/bin/perl

#get new variants from gemini database output

open(LIST,$ARGV[1]);
open(GEMINI,$ARGV[0]);

%variants;

while(<LIST>)
{
    chomp;
    $variants{$_}=1;
}


while (<GEMINI>)
{
    chomp;
    @arr = split ("\t",$_);
    $arr[0] =~ s/chr//;
    $id=$arr[0]."-".$arr[2]."-".$arr[4]."-".$arr[5];
    print $id."\n";
    if(exists($variants{$id}))
    {
	print $_."\n";
    }
}

close(OMIM);
close(GEMINI);