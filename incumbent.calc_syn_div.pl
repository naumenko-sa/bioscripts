#!/usr/bin/perl
#calcs table 1 by yegor method
#depends on fasta2fourfold.pl
#usage calc_syn_div.pl file.fourfold [overall|parallel] [AC|AG|AT|CG|CT|GT]
use Switch;

open(IN,$ARGV[0]);

$HUMAN=0;
$MOUSE=1;
$RAT=2;
$DOG=3;

$total_sites=0;
$div_sites=0;

$statistic=$ARGV[1];

@nuc=split(//,$ARGV[2]);

while(<IN>)
{
    chomp;
    @ar = split(' ',$_);
    $site = $ar[3];
    if( $ar[2] =~ '([ATGC])\g1{3}' and $ar[4]  =~ '^([ATGC])\g1{3}' )
    {
	$m = substr($site,$MOUSE,1);
	$r = substr($site,$RAT,1);
	$h = substr($site,$HUMAN,1);
	$d = substr($site,$DOG,1);
	
	switch ($statistic)
	{
	    case "overall"
	    {
		if (is_proper_pair($m,$r))
		{
		    if ($m ne $r) #mouse-rat divergence
		    {
			$div_sites++;
		    }
		    $total_sites++;
		}
	    }
	    case "parallel"
	    {
		if (($h ne $d) and is_proper_pair($h,$d))
		{
		    if (is_extended_proper_pair($m,$r))
		    {
			$total_sites++;
			if ($m ne $r)
			{
			    $div_sites++;
			}
		    }
		}
	    }
	}
    }
}

print $div_sites."\t".$total_sites."\t".$div_sites/$total_sites."\n";

close(IN);

sub is_extended_proper_pair()
{
    @nuc=split(//,$ARGV[2]);
        $result = 0;
            $a=$_[0];
                $b=$_[1];
                    if ( is_proper_pair($a,$b) or (($a eq $b) and ($a eq $nuc[0] or $a eq $nuc[1])))
                        {
                    	$result = 1;
                    	    }
                    	        return $result;
                    	        }

sub is_proper_pair()
{
    @nuc=split(//,$ARGV[2]);
    $result = 0;
    $a=$_[0];
    $b=$_[1];
    if ((($a.$b) eq ($nuc[0].$nuc[0])) or
	(($a.$b) eq ($nuc[1].$nuc[1])) or 
	(($a.$b) eq ($nuc[0].$nuc[1])) or
	(($a.$b) eq ($nuc[1].$nuc[0])))
    {
	$result = 1;
    }
    return $result;
}

