#!/usr/bin/perl

if (@ARGV != 2)
{
    print "Converts blastn output to coverage\n";
    print "Usage: blastn_to_cov.pl file.blastn reference.length\n";
    exit(1);
}

open(IN,$ARGV[0]);

my $len = $ARGV[1];

my @cov;

for ($i=1;$i<=$len;$i++)
{
    push(@cov,0);
}

while(<IN>)
{
    chomp;
    
    my $left=0;
    my $right=0;
    
    @ar = split(' ',$_);
    if ($ar[6] >= $ar[7])
    {
	$left = $ar[7];
	$right = $ar[6];
    }
    else
    {
	$left=$ar[6];
	$right=$ar[7];
    }
    #print $left."\t".$right."\n";
    for ($i=$left;$i<=$right;$i++)
    {
	$cov[$i]++;
    }
}

for ($i=1;$i<=$len;$i++)
{
    print $i."\t".$cov[$i]."\n";
}

close(IN);