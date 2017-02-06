#!/usr/bin/perl

#################################################################
#    default : show number of internal stop codons TAA, TAG,TGA #
#    alignment.detect_stops.pl file.fasta			#
#    --print_pos -also show positions of stops                  #
#################################################################

sub usage()
{
    print "print the number (and position) of internal stops in the alignment\n";
    print "alignment.detect_stops.pl file.fasta [--print_pos]\n";
}

use strict;

if (@ARGV==0)
{
    usage();
    exit(0);
}

my %align;
open (IN,$ARGV[0]) || die();
while (<IN>)
{
    chomp;
    my $name = $_;
    my $s = <IN>;
    chomp($s);
    $align{$name}=$s;
}
close(IN);

my %to_del;
#search for internal stops
my $istops=0;
foreach my $name (keys %align)
{
    for (my $i=0;$i < length($align{$name})/3 - 1;$i++) #except of the last triplet
    {
	my $codon = uc(substr($align{$name},$i*3,3));
	if (($codon eq "TAG") or ($codon eq "TAA") or ($codon eq "TGA"))
	{
	    if ($ARGV[1] eq '--print_pos')
	    {
		print "Stop: ".$name."\t".3*$i."\n";
	    }
	    $istops++;
	}
    }
}
print $istops."\n";