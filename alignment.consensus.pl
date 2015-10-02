#!/usr/bin/perl
######################################################################
#     alignment consensus of many sequences
#     gaps and N's are not counted if something else present          
#     first sequence has a priority when allele frequencies are equal 
#######################################################################
use strict;

my @sp;
my %al;
my $consensus='';

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

my %site;
for (my $i=0;$i < length($al{$sp[0]});$i++)
{
    my $first_a = substr($al{$sp[0]},$i,1);
    foreach my $s (@sp)
    {
	my $a = substr($al{$s},$i,1);
	$site{$a}++;
    }
    delete $site{'N'};
    delete $site{'-'};
    my $c = largest_hash_value(\%site);
    if ($c eq '')
    {
	$c = substr($al{$sp[0]},$i,1);
    }
    $consensus .= $c;
    %site = ();
}
print ">consensus\n";
print $consensus."\n";

sub largest_hash_value (\%)
{
    my $hash = shift;
    keys %$hash;       # reset the each iterator

    my ($large_key, $large_val) = each %$hash;

    while (my ($key, $val) = each %$hash)
    {
	if ($val > $large_val) 
        {
            $large_val = $val;
            $large_key = $key;
        }
    }
    return $large_key;
}
