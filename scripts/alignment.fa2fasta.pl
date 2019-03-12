#!/usr/bin/perl
use strict;

# get rid of unnecessary EOL's in an alignment

my $id; my $seq;

open (INFILE, $ARGV[0]);

$id = <INFILE>; 
chomp $id;
$seq = "";

while (<INFILE>) {
    chomp $_;
    if (/>/) 
    {
	print $id."\n";
	print $seq."\n";
	$id = $_;
	$seq = "";
    }
    else 
    {
	$seq = $seq.$_;
    }
}

print $id."\n";
print $seq."\n";
