#!/usr/bin/perl
# oncatenate fasta files:
#>gam1    +   >gam1    =  >gam1
#aaa	         bbb		aaabbb

#chdir("/mnt/home/snaumenko1/project_gamarus/5tree_gam13_pulex/5.2.make_orf/prot_fasta"); 
#chdir("/mnt/home/snaumenko1/project_gamarus/6dnds");

my @files = <*>;
my %align;
foreach $file (@files)
{
    #print $file."\n";
    open(IN,$file);
    while(<IN>)
    {
	chomp;
	my @ar = split (/\|/);
	my $name = substr($ar[0],1,length($_)-1);
	#print $name."\n";
	my $dna = <IN>;
	chomp $dna;
	#print $name."\t".length($dna)."\n";
	$align{$name}.=$dna;
	#print $name."=>".length($align{$name})."\n";
    }
    close($file);
}

foreach $key (sort keys %align)
{
    print ">".$key."\n";
    print $align{$key}."\n";
}
