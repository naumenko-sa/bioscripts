#!/usr/bin/perl

#change fasta contig headers to tba format
#contig length is usually wrong

my $i=1;
open(IN,$ARGV[0]);
my $name=$ARGV[0];

$name =~ s/.fasta//;

open(OUT,">".$name);
while(<IN>)
{
    chomp;
    my $title = $_;
    my $dna = <IN>;
    chomp $dna;
    if ($_ =~ />/)
    {
	my @ar=split('_',$_); 
	print OUT '>'.$name.":".$i.":1:+:".length($dna)."\n";
	print OUT $dna."\n";
	$i++;
    }
}
close(IN);
close(OUT);