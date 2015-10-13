#!/usr/bin/perl

#use after orthomcl steps
#extract alignements from compliant fasta and groups.txt

#usage orthomclToAlignements.pl groups.txt goodproteins.fasta

%proteins;

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

while(<GROUPS>)
{
    chomp;
    @ar = split(' ',$_);
    $file = $ar[0];
    $file =~ s/:/.fasta/;
    #print $file."\n";
    open(OUT,">$file");
    for($i=1;$i<@ar;$i++)
    {
	print OUT ">".$ar[$i]."\n";
	print OUT $proteins{$ar[$i]}."\n";
    }
    close(OUT);
}

close(GROUPS);