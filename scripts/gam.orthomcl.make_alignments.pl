#!/usr/bin/perl

use List::MoreUtils qw(uniq);

# use after orthomcl steps
# extract alignments from compliant fasta and groups.txt
# prints only 1:1 orthologs

#usage gam.orthomclToAlignements.pl groups.txt goodproteins.fasta number_of_species

%proteins;

$number_of_species = $ARGV[2];

open(PROT,$ARGV[1]);
while(<PROT>)
{
    chomp;
    $name = substr($_,1,length($_)-1);
    $pro = <PROT>;
    chomp $pro;
    $proteins{$name}=$pro;
}
close(PROT);

open(GROUPS,$ARGV[0]);

my $cog_count = 1;

while(<GROUPS>)
{
    chomp;
    @ar = split(' ',$_);
    if (@ar == $number_of_species+1)
    {
	@ar_entries;
	for ($i=1;$i<@ar;$i++)
	{
	    @ar_tmp = split(/\|/,$ar[$i]);
	    push (@ar_entries,$ar_tmp[0]);
	}
	@unique_entries = uniq(@ar_entries); 
	
	if (@unique_entries == $number_of_species)
	{
	    #$file = $ar[0];
	    $file = "cog".$cog_count.".fasta";
	    $cog_count++;
	    #$file =~ s/:/.fasta/;
	    #print $file."\n";
	    open(OUT,">$file");
    	    for($i=1;$i<@ar;$i++)
	    {
		print OUT ">".$ar[$i]."\n";
		print OUT $proteins{$ar[$i]}."\n";
	    }
	    close(OUT);
	}
    }
}

close(GROUPS);