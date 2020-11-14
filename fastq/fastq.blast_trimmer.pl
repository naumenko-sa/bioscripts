#!/usr/bin/perl

#trims fastq by blast results
#in blastn result may be lesser reads then in fastq

if (@ARGV != 2)
{
    print "Trim reads using blast result\n";
    print "Usage: fastq_blast_trimmer.pl file.fq blastn.result\n";
    exit(1);
}
            

open (FQ,$ARGV[0]);
open (BLAST,$ARGV[1]);

%from;
%to;

while (<BLAST>)
{
    my @ar = split ("\t",$_);
    #print $ar[4]."\t".$ar[5]."\n";
    $read_name = $ar[0];
    $from{$read_name} = $ar[4];
    $to{$read_name} = $ar[5];
}    

while(<FQ>)
{
    $s1 = $_;
    $read_name=substr($_,1,length($_)-1);
    chomp $read_name;
    $s2 = <FQ>;
    $s3 = <FQ>;
    $s4 = <FQ>;

    if (exists($from{$read_name}))
    {    
	print $s1;
	print substr($s2,$from{$read_name}-1,$to{$read_name}-$from{$read_name}+1)."\n";
	print $s3;
	print substr($s4,$from{$read_name}-1,$to{$read_name}-$from{$read_name}+1)."\n";
    }
}

close(BLAST);
close(FQ);