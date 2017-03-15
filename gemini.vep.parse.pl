#!/usr/bin/perl

# parses vcf with VEP annotation - CSQ field
# works for VEP - refseq !
# for ensembl fields are different:
# f_exon = 9
# f_hgvsc = 16
# f_hgvsp = 17

$f_gene=3;
$f_exon=6;
$f_hgvsc=13;
$f_hgvsp=14;

open(IN,$ARGV[0]);

print "superindex\tInfo_refseq_no_gene\n";

while (<IN>)
{
    chomp;
    $impacts_string='';
    if ($_ !~ /CHROM/)
    {
	@ar = split("\t",$_);
	$id = "chr".$ar[0].":".$ar[1]."-".$ar[2]."-".$ar[3];
	@ar2 = split(",",@ar[4]);
	#print $id."\n";
	foreach $impact (@ar2)
	{
	    #print $impact."\n";
	    @ar3= split(/\|/,$impact);
	    #if there is hvgsc
	    if ($ar3[$f_hgvsc] ne '')
	    {
		if ($ar3[$f_hgvsp] eq '')
		{
		    $ar3[$f_hgvsp] = 'NA';
		}
		if ($ar3[$f_exon] eq '')
		{
		    $exon='NA';
		}
		else
		{
		    $exon = $ar3[$f_exon];
		    $exon =~ s/"\-"/"NA"/;
		}
		$impacts_string .= "exon".$exon.":".$ar3[$f_hgvsc].":".$ar3[$f_hgvsp].",";
	    }
	}
	print $id."\t".substr($impacts_string,0,-1)."\n";
    }
}

close(IN);