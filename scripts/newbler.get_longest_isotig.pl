#!/usr/bin/perl

#gets longest isotig from newbler
#headers in fasta file should be 
#>contig00021_isogroup00001_691
#so delete length= and gene= from newbler output

open(IN,$ARGV[0]);

%lengths;
%seqs;

while (<IN>)
{
    $name = $_;
    chomp $name;
    $dna = <IN>;
    chomp $dna;
    @ar = split("_",$name);
    $isoname=$ar[1];
    $isolen=$ar[2];
    #print $isoname."\t".$isolen."\n";
    if (!exists $lengths{$isoname})
    {
	$lengths{$isoname} = $isolen;
	$seqs{$isoname}= $dna;
    }
    else
    {
	if ($isolen > $lengths{$isoname})
	{
	    $lengths{$isoname} = $isolen;
	    $seqs{$isoname} = $dna;
	}
    }
}

foreach $key (sort keys %lengths)
{
    print ">".$key."_".$lengths{$key}."\n";
    print $seqs{$key}."\n";
}

close(IN);