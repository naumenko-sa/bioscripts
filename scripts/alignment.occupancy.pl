#!/usr/bin/perl
#calculates the number of columns without gaps in alignment
use strict;

my @sp;
my %al;

open (IN,$ARGV[0]) || die();
while (<IN>)
{
    chomp;
    my $name=$_;
    push (@sp,$name);
    #print $_."\n";
    my $s = <IN>;
    #print $s;
    chomp($s);
    $al{$name}=$s;
}
close(IN);

my %to_del;
#search for  gaps
my $without_gaps=0;
COLUMN: for (my $i=0;$i < length($al{$sp[0]});$i++)
{
    foreach my $s (@sp)
    {
	my $codon = substr($al{$s},$i,1);
	if ($codon =~ m/-/i)
	{
	    next COLUMN;
	}
    }
    $without_gaps++;
}
my $ratio=$without_gaps/length($al{$sp[0]});
printf("%.3f\n",$ratio);