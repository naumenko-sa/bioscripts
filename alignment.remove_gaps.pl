#!/usr/bin/perl
#deletes nucleotide sites with gaps in any sequence
use strict;


my @sp;
my @nuc;
my $original_length;
open (IN,$ARGV[0]) || die();
while (<IN>)
{
    chomp;
    push (@sp,$_);
    #print $_."\n";
    my $s = <IN>;
    #print $s;
    chomp($s);
    push (@nuc,$s);
}
close(IN);

$original_length = length($nuc[0]);
my @to_del;
my @to_del_tripl;
#search for gaps
foreach my $str (@nuc)
{
    for (my $i=0;$i < length($str);$i++)
    {
    	my $c = substr($str,$i,1);
        if ($c =~ m/-/)
	{
    	    push (@to_del,$i);
        }
    }
}

my %unique = ();
foreach my $item (@to_del)
{
    $unique{$item} ++;
}
my @pos_del = keys %unique;


#cutting
foreach my $dna (@nuc)
{
        foreach my $pos (@pos_del)
        {
	    my $t = substr($dna,$pos,1,'z');
        }
        my $t = $dna =~ s/z//g;
}
                                
open(OUT,">$ARGV[0].trimmed");
for(my $i=0;$i<@sp;$i++)
{
        print OUT $sp[$i]."\n";
        print OUT $nuc[$i]."\n";    
}
close (OUT);
