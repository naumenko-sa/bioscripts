#!/usr/bin/perl

#split genome reference for maf2fasta utility

open(IN,$ARGV[0]);

my $basename=$ARGV[0];
my $i=1;
while(<IN>)
{
    my $s1 = $_;
    my $s2 = <IN>; 
    open(OUT,">".$basename.".".$i);
    print OUT $s1;
    print OUT $s2;
    close(OUT);
    $i++;
}

close(IN);