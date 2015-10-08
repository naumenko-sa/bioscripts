#!/usr/bin/perl
#deletes all columns with - and conservative columns in aa alignment

use strict;


my @sp;
my %al;

open (IN,$ARGV[0]) || die();
while (<IN>)
{
    chomp;
    my $name=$_;
    push (@sp,$name);
    #print $_."\n";
    my $s = <IN>;
    #print $s;
    chomp($s);
    $al{$name}=$s;
}
close(IN);

my %to_del;
#search for  gaps
for (my $i=0;$i < length($al{$sp[0]});$i++)
{
    my @aa;
    foreach my $s (@sp)
    {
	my $codon = substr($al{$s},$i,1);
	if ($codon =~ m/-/i)
	{
	    $to_del{$i}=1;    	    
    	}
    	push (@aa,$codon);
    }
    my @un = uniq(@aa);
    if (@un == 1 or @un > 2)
    {
	$to_del{$i}=1;
    }
    @aa=();
}
#cutting
my %res_align;
my $res;
my $prev=0;

for (my $i=0;$i<length($al{$sp[0]});$i++)
{
    if(!exists($to_del{$i}))
    {
	foreach my $s (@sp)
	{
	    $res_align{$s} .= substr($al{$s},$i,1);
	}
    }
}
                        
foreach my $s (@sp)
{
    print $s."\n";
    print $res_align{$s}."\n";    
}

sub uniq {
  my %seen;
      return grep { !$seen{$_}++ } @_;
}