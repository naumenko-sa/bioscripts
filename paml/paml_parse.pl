#!/usr/bin/perl

# parses paml dnaml output for ancestral states

#tree with node labels for Rod Page's TreeView
#((((((((1_hg18, 2_mm9) 17 , 3_canFam2) 16 , 4_loxAfr2) 15 , 5_monDom4) 14 , 6_ornAna1) 13 , 7_galGal3) 12 , 8_xenTro2) 11 , 9_danRer5) 10 ;
#Nodes 10 to 17 are ancestral

use strict;

my $paml_filename = $ARGV[0];
open (IN,$paml_filename);

my $line;
#does not influence output by fact
my @print_order = ("hg18","17","16","15","14","13","12","11","10","mm9","canFam2","loxAfr2","monDom4","ornAna1","galGal3","xenTro2","danRer5");
my @align;
do
{
    $line = <IN>;
    #print $line."\n";
}
until ($line =~ /site(.*)Freq(.*)Data/);

while ($line = <IN>)
{
    chomp $line;
    if ($line =~ /Summary/)
    {
	last;
    }
    my $results;
    if ($line =~ /[atgc-]/i)
    {
	my @tokens = split(" ",$line);
	my $n = scalar @tokens;
	$results = substr($tokens[2],0,9);
	my $anc;
	for (my $i=8;$i>=1;$i--)
	{
	    $anc.=substr($tokens[2+$i],0,1);
	}
	my $res = substr($results,0,1).$anc.substr($results,1,8);
	for (my $i=0;$i < scalar @print_order;$i++)
	{
	    $align[$i].=substr($res,$i,1);
	}
    }
}
close (IN);

for (my $j=0; $j < scalar @print_order;$j++)
{
    print ">".$print_order[$j]."\n";
    print $align[$j]."\n";
}
