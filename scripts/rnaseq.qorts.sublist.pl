#!/usr/bin/perl

#FILE1 - list key - value
#1-10895-10969	1
#1-14829-14929	3
#1-15038-15310	1
#1-16310-16606	2
#1-17055-17605	3
#1-17055-17914	2
#1-92240-165883	3
#1-129223-142946	1
#1-135169-138163	1
#1-136708-136903	2

#FILE2 - list of keys
#1-10895-10969
#1-14829-14929
#1-15038-15310

#RESULT = key - value list for keys found in FILE2

open(IN,$ARGV[0]);
open(IN1,$ARGV[1]);

%list;
@array;

while (<IN>)
{
    $str = $_;
    chomp $str;
    @ar = split("\t",$str);
    $list{$ar[0]}=$ar[1];
}

while (<IN1>)
{
    chomp $_;
    push(@array,$_);
}

foreach $key (@array)
{
    print $key."\t".$list{$key}."\n";
}

close(IN);
close(IN1);