#!/usr/bin/perl

#trim right n nucleotides from fastq file

if (@ARGV != 2)
{
    print "Trim right N nucleotides from fastq file on the right\n";
    print "Usage: trim_right.pl file.fastq N\n";
    exit(1);
}

open(IN,$ARGV[0]);

$i=0;
while(<IN>)
{
    if ($i % 2 == 1)
    {
	print substr($_,0,length($_)-$ARGV[1])."\n";
    }
    else
    {
	print $_;
    }
    $i++;
}

close(IN);