#!/usr/bin/perl
#deletes all triplets with stops and gaps (deletions)
#removes 5 codons around deletions >=2 nucleotides

use strict;

my @sp;
my @nuc;
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

my %to_del;
#search for internal stops
foreach my $str (@nuc)
{
    for (my $i=0;$i < length($str)/3;$i++)
    {
	my $codon = substr($str,$i*3,3);
	
	if ($codon =~ m/TAG|TAA|TGA/i)
	{
	    $to_del{3*$i}=1;    	    
    	}
    	elsif ($codon =~ m/-/)
    	{
    	    $to_del{3*$i}=1;
    	    
    	    my $c = ($codon =~ s/-/-/g);
    	    
    	    #print $codon."\t".$c."\n";
    	    if ($c > 1)
    	    {
    		for (my $z=$i-5;$z<=$i+5;$z++)
    		{
    		    if ($z>=0)
    		    {
    			$to_del{3*$z}=1;
    		    }
    		}
    	    }
    	}
    }
}
        
#cutting
my @res_align;
my $res;
my $prev=0;
foreach my $pos (sort {$a <=> $b} keys %to_del)
{
    for (my $i=0;$i<@sp;$i++)
    {
	$res_align[$i].=substr($nuc[$i],$prev,$pos-$prev);
    }
    #skip codons
    $prev=$pos+3;
}
for (my $i=0;$i<@sp;$i++)
{
    $res_align[$i].=substr($nuc[$i],$prev,length($nuc[0])-$prev);
}
                        
for(my $i=0;$i<@sp;$i++)
{
    print $sp[$i]."\n";
    print $res_align[$i]."\n";    
}