#!/usr/bin/perl
use POSIX;

# ag ct

open(ALN,$ARGV[0]);

@aln;

while (<ALN>)
{
    chomp;
    my $name = $_;
    my $dna = <ALN>;
    chomp ($dna);
    push(@aln,$dna);
}
close(ALN);

my @tt = get_tr_trv($aln[0],$aln[1]);
$count=0;

print "Tr:".$tt[0]."\t"."Trv:".$tt[1]."\t"."tr/trv:";

if ($tt[1] != 0)
{
    printf ("%.3f",$tt[0]/$tt[1]);
}
else
{
    print "Inf";
}
print "\n";

sub hd 
{
    my ($k, $l) = @_;
    my $diff = $k ^ $l;
    my $num_mismatch = $diff =~ tr/\0//c;
}

sub get_tr_trv()
{
    my $tr=0;
    my $trv=0;
    my ($a,$b) = @_;
    for (my $i=0;$i < length($a);$i++)
    {
	my $x = substr($a,$i,1);
	if ($x =~ /A/i or $x =~ /G/i)
	{
	    $t = "PU";
	}
	else
	{
	    $t = "PY";
	}
	
	my $y = substr($b,$i,1);
	if ($y =~ /A/i or $y =~ /G/i)
	{
	    $z = "PU";
	}
	else
	{
	    $z = "PY";
	}
	if ($x ne $y and $x ne "-" and $y ne "-")
	{
	    if ($t eq $z)
	    {
		$tr++;
	    }
	    else
	    {
		$trv++;
	    }
	}
    }
    return ($tr,$trv);
}