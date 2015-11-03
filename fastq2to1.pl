#!/usr/bin/perl
# to merge 2 fastq files for left and right read into single file
#usage: fastq2to1.pl left.fastq right.fastq single.fastq

use strict;

open(LEFT,$ARGV[0]);
open(RIGHT,$ARGV[1]);
open(OUT,">".$ARGV[2]);

while(<LEFT>)
{
    chomp;
    my $l1 = $_;
    my $l2 = <LEFT>;
    chomp $l2;
    my $l3 = <LEFT>;
    chomp $l3;
    my $l4 = <LEFT>;
    chomp $l4;
    
    my $r1 = <RIGHT>;
    chomp $r1;
    my $r2 = <RIGHT>;
    chomp $r2;
    my $r3 = <RIGHT>;
    chomp $r3;
    my $r4 = <RIGHT>;
    chomp $r4;
    
    print OUT $l1."\n";
    print OUT $l2."\n";
    print OUT $l3."\n";
    print OUT $l4."\n";
    
    print OUT $r1."\n";
    print OUT $r2."\n";
    print OUT $r3."\n";
    print OUT $r4."\n";
}

close(LEFT);
close(RIGHT);
close(OUT);
