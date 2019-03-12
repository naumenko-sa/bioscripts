#!/usr/bin/perl

#prints alignment consensus of 2 sequences 
#gaps and N's are not counted if something else present
#first sequence has a priority when alleles are different

# $1 - alignment.fasta

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

for (my $i=0;$i < length($al{$sp[0]});$i++)
{
    my $a = substr($al{$sp[0]},$i,1);
    my $b = substr($al{$sp[1]},$i,1);
    
    my $c = $a;
    
    if (($a eq 'N' or $a eq '-') and ($b ne 'N' and $b ne '-'))
    {
    	$c = $b;
    }
    
    $consensus .= $c;
}
print ">consensus\n";
print $consensus."\n";