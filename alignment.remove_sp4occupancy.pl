#!/usr/bin/perl
#calculates a = the proportion of columns without gaps in alignment
#removes species with greater number of gaps until p>threshold (default=0.5)
#written for gammarus project 5.03.2015

my @sp;
my %al;
my %al_stat;

open (IN,$ARGV[0]) || die();
my $threshold = 0.5;

if (@ARGV == 2)
{
    $threshold = $ARGV[1];
}

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

while(calc_gaps_in_columns(\@sp,\%al) < $threshold)
{
    calc_gaps_per_specie(\@sp,\%al,\%al_stat);
    #print $al_stat{$sp[0]}."\n";
    my $bad_sp = largest_hash_value(\%al_stat);
    #print $bad_sp."\n";

    delete $al_stat{$bad_sp};
    delete $al{$bad_sp};
    my @del_indexes = reverse(grep { $sp[$_] eq $bad_sp } 0..$#sp);
    foreach $item (@del_indexes) {
       splice (@sp,$item,1);
    }
}

#printf "%.3f\n",calc_gaps_in_columns(\@sp,\%al);

foreach $sp (@sp)
{
    print $sp."\n";
    print $al{$sp}."\n";
}

close(IN);

sub calc_gaps_per_specie(\@,\%,\%)
{
    my ($p_sp,$p_al,$p_al_stat) = @_;
    foreach my $sp (@$p_sp)
    {
	my $count = ($$p_al{$sp} =~ s/-/-/g);
        $$p_al_stat{$sp}=$count;
    }
}

sub largest_hash_value (\%)
{
    my $hash = shift;
    keys %$hash;       # reset the each iterator
        
    my ($large_key, $large_val) = each %$hash;
            
    while (my ($key, $val) = each %$hash) 
    {
	if ($val > $large_val) 
	{
            $large_val = $val;
            $large_key = $key;
        }
    }
    return $large_key;
}

sub calc_gaps_in_columns()
{
    my ($p_sp,$p_al) =@_;
    my $without_gaps=0;
    COLUMN: for (my $i=0;$i < length($$p_al{$$p_sp[0]});$i++)
    {
	foreach my $s (@$p_sp)
	{
	    my $codon = substr($$p_al{$s},$i,1);
	    if ($codon =~ m/-/i)
	    {
		next COLUMN;
	    }
	}
	$without_gaps++;
    }
    my $ratio=$without_gaps/length($$p_al{$$p_sp[0]});
    return $ratio;
}