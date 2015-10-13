#!/usr/bin/perl

#substites names gam0 ... gam10 to full names
#usage subst_names.pl species.list tree.newick

if (@ARGV < 2)
{
    print "Usage: subst_names.pl species.list tree.newick\n";
    exit(0);
}

my %names;
open(IN,$ARGV[0]);

while(<IN>)
{
    chomp;
    my @ar = split(" ",$_,2);
    $names{$ar[0]} = $ar[1];
}

close(IN);

my $tree='';
open (TREE,$ARGV[1]);
while (<TREE>)
{
    chomp;
    $tree .= $_;
}
#print $tree;
my @keys = sort {length $b <=> length $a}keys %names;

foreach $key (@keys)
{
    #print $key."=>".$names{$key}."\n";
    $tree =~ s/$key/$names{$key}/g;
}

print $tree;

close(TREE);