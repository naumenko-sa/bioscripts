#!/usr/bin/perl

#average pairwise distance between sequence in an alignment

use POSIX;

open(ALN,$ARGV[0]);

@aln;
@names;
while (<ALN>)
{
    chomp;
    my $name = $_;
    chomp $name;
    push(@names,substr($name,1));
    my $dna = <ALN>;
    chomp ($dna);
    push(@aln,$dna);
}
close(ALN);

$dist=0;
$count=0;

my $pr=0;
if ($ARGV[1] == 1)
{
    $pr=1;
    print "\t";
    for ($i=0;$i<@aln;$i++)
    {
	print $names[$i]."\t";
    }
    print "|AVERAGE|\n";
}

for ($i=0;$i<@aln;$i++)
{
    my $av_dist=0;
    if ($pr)
    { print $names[$i]."\t";}
    for ($j=0;$j<@aln;$j++)
    {
	$pdist = hd($aln[$i],$aln[$j]) / length($aln[0]);
	$av_dist+=$pdist;
	if ($pr)
	{
	    #print $i."<=>".$j.":";
	    printf("%.3f",$pdist);
	    print "\t";
	}
	$dist += $pdist;
	$count++;
    }
    if($pr)
    {
	printf("|%.3f|",$av_dist/@aln);
	print "\n";
    }
}

printf ("%.3f",$dist/$count);
print "\n";

sub hd 
{
    my ($k, $l) = @_;
    my $diff = $k ^ $l;
    my $num_mismatch = $diff =~ tr/\0//c;
}