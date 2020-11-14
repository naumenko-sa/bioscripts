#!/usr/bin/perl

#trim left n nucleotides from fastq file

if (@ARGV != 2)
{
    print "Trim left N nucleotides from fastq file\n";
    print "Usage: trim_left.pl file.fastq N\n";
    exit(1);
}

open(IN,$ARGV[0]);

$i=0;
while(<IN>)
{
    if ($i % 2 == 1)
    {
	print substr($_,$ARGV[1],length($_)-$ARGV[1]);
    }
    else
    {
	print $_;
    }
    $i++;
}

close(IN);