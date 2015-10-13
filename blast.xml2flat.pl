#!/usr/bin/perl

#PBS -d .

use strict;
use Bio::SearchIO; 
my $in = new Bio::SearchIO(-format => 'blastxml', 
                           -file   => 'gam3_1.fasta_vs_nr.xml');
while( my $result = $in->next_result ) {
  ## $result is a Bio::Search::Result::ResultI compliant object
  while( my $hit = $result->next_hit ) {
    ## $hit is a Bio::Search::Hit::HitI compliant object
    while( my $hsp = $hit->next_hsp ) {
      ## $hsp is a Bio::Search::HSP::HSPI compliant object
#      if( $hsp->length('total') > 50 ) {
#        if ( $hsp->percent_identity >= 75 ) {
          print "Query=",   $result->query_description,
            " Hit=",        $hit->name,
            " evalue=",     $hsp->evalue,
           # " aln=", $hsp->get_aln,
            " start query=", $hsp->start('query'),
            " end query=", $hsp->end('query'),
            " start hit=", $hsp->start('hit'),
             " end hit=", $hsp->end('hit'),
            " score=", $hsp->score,  "\n";
        }
      }
    }  

