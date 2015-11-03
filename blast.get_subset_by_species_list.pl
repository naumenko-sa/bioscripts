#!/usr/bin/perl

#get part of blast result file having particular species in the first and the second column
#useage blast.get_subset_by_species_list.pl [result.blast] [species.list]

%target_species;

open(IN,$ARGV[1]);

while(<IN>)
{
    chomp;
    $target_species{$_}=1;
}

close(IN);

open(IN,$ARGV[0]);
while(<IN>)
{
    chomp;
    @columns=split("\t",$_);
    @sp1 = split(/\|/,$columns[0]);
    @sp2 = split(/\|/,$columns[1]);
    if (exists($target_species{$sp1[0]}) && exists($target_species{$sp2[0]}))
    {
	print $_."\n";
    }
}
close(IN);