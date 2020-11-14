#!/usr/bin/perl

#checks where is adapter (r1, r2, both reads)
#usage: mp_check_adapter.pl r1.fq f2.fq adapter-seq [print|columns] [what-to-print: 00,10,01,11]
#print prints entire reads to r1.fq.filt r2.fq.filt
#columns prints reads as seq1	seq2 to stdout (for debug)

if (@ARGV <3 )
{
    print "Usage: mp_check_adapter.pl r1.fq r2.fq adapter-seq [print|columns] [what-to-print: 00,10,01,11]\n";
}

open(R1,$ARGV[0]);
open(R2,$ARGV[1]);

$pat = $ARGV[2];

#stat | print | columns
$mode="stat";

if (@ARGV == 5)
{
    if ($ARGV[3] eq "print")
    {
	$mode ="print";
    }
    elsif ($ARGV[3] eq "columns")
    {
	$mode = "columns";
    }
    else
    {
	print "wrong usage";
	exit(1);
    }
}

if ($mode eq "print")
{
    open(FR1,">".$ARGV[0].".filt");
    open(FR2,">".$ARGV[1].".filt");
}

%stat;

@r1 = ();
@r2 = ();

while(<R1>)
{
    chomp;
    push (@r1, $_);
    $sr1 = <R1>;
    chomp $sr1;
    push (@r1,$sr1);
    $sr2 = <R2>;
    chomp $sr2;
    push (@r2,$sr2);
    $sr2 = <R2>;
    chomp $sr2;
    push (@r2,$sr2);
    
    $m1 = ( $sr1 =~ /$pat/ ) ? 1:0;
    
    $m2 = ( $sr2 =~ /$pat/ ) ? 1:0; 
    
    $res = $m1.$m2;
    
    if ($mode eq "columns" && $res eq $ARGV[4])
    {
	print $sr1."\t".$sr2."\n";
    }
    $stat{$m1.$m2}++;
    
    $sr1 = <R1>;
    chomp $sr1;
    push (@r1,$sr1);
    $sr1 = <R1>;
    chomp($sr1);
    push (@r1,$sr1);
    
    $sr2 = <R2>; chomp $sr2; push (@r2,$sr2);
    $sr2 = <R2>; chomp $sr2; push (@r2,$sr2);
    
    if ($mode eq "print" && $res eq $ARGV[4])
    {
	for ($i=0;$i<4;$i++)
	{
	    print FR1 $r1[$i]."\n";
	    print FR2 $r2[$i]."\n";
	}
    }
    @r1 = ();
    @r2 = ();
}

close(R1);
close(R2);

if (@ARGV == 5 && $ARGV[3] eq "print")
{
    close(FR1);
    close(FR2);
}

if ($mode eq "stat")
{
    foreach $key (keys %stat)
    {
	print $key."\t".$stat{$key}."\n";
    }
}