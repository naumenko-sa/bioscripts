#!/usr/bin/perl
#removes sequences with internal stops save order

use strict;

my @names;

my %align;
open (IN,$ARGV[0]) || die();
while (<IN>)
{
    chomp;
    my $name = $_;
    my $s = <IN>;
    chomp($s);
    $align{$name}=$s;
    push (@names,$name);
}
close(IN);

my %to_del;
#search for internal stops
my $istops=0;
foreach my $name (@names)
{
    my $has_stop=0;
    for (my $i=0;$i < length($align{$name})/3 - 1;$i++) #except of last triplet
    {
	my $codon = substr($align{$name},$i*3,3);
	if ($codon =~ m/TAG|TAA|TGA/i)
	{
	    #print "Stop: ".$name."\t".3*$i."\n";
	    $has_stop++;
	}
    }
    if ($has_stop == 0)
    {
	    print $name."\n";
	    print $align{$name}."\n";
    }
}
#print "Total: ".$istops."\n";