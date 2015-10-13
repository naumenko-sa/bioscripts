#!/usr/bin/perl

#reports everything about two seqs: substitutions, types, 1,2,3 pos

open(ALN,$ARGV[0]);

@dna;

$d = <ALN>;
$d = <ALN>;
chomp $d;
push(@dna,$d);

$d = <ALN>;
$d = <ALN>;
chomp $d;
push (@dna,$d);

%by_type;

@pos = (0,0,0);

$len=0;

for (my $i=0;$i<length($dna[0]);$i++)
{
    $a = substr($dna[0],$i,1);
    $b = substr($dna[1],$i,1);
    
    if ($a ne "-" and $b ne "-")
    {
	$len++;
	if ($a ne $b)
	{
	    @type = ($a,$b);
	    @type = sort @type;
	    $by_type{$type[0].$type[1]}++;
	    $pos[ $i % 3 ]++;
	}
    }
}

close(ALN);

print "Len:\t".$len."\n";

for (my $i=0; $i < 3; $i++)
{
    $k=$i+1;
    print "Pos ".$k."\t".$pos[$i]."\n";
}

foreach my $key (sort keys %by_type)
{
    print $key."\t".$by_type{$key}."\n";
}