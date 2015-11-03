#!/usr/bin/perl

#get some records for fasta file where 
#species name is the first field >gam11|gene100
#useage fasta.someRecordsField.pl [file.fasta] [species.list]

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
    $gene_name=$_;
    @sp1 = split(/\|/,$_);
    $sp_name = substr($sp1[0],1);
    #print $sp_name."\n";
    if (exists($target_species{$sp_name}))
    {
	print $gene_name."\n";
	$dna=<IN>;
	print $dna;
    }
}
close(IN);