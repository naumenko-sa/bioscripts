#!/usr/bin/perl 
 
# calculate mean and CIs from bootstrapped distribution 
# produced by MK_v2.pl

# v. 2: in addition to bootstrap percentiles, output also all the stats in one line


use strict;
use Statistics::Descriptive;

# resampled values
my (@dnds, @pnps_freq, @pnps_all, @alpha_freq, @alpha_all);

my $our_part_almost_started;
my $our_part_started;


# real values
my ($dnds_r,  $pnps_freq_r, $pnps_all_r, $alpha_freq_r, $alpha_all_r);

# load distibution
while (<>) {
    chomp;

    # Get actual data
    my $delme;
    if (/dN\/dS/) {
	($delme, $dnds_r) = split /:\s+/;
    }
    if (/pN\/pS, all/) {
	($delme, $pnps_all_r) = split /:\s+/;
    }
    if (/pN\/pS, frequent/) {
	($delme, $pnps_freq_r) = split /:\s+/;
    }
    if (/Alpha, all/) {
	($delme, $alpha_all_r) = split /:\s+/;
    }
    if (/Alpha, freq/) {
	($delme, $alpha_freq_r) = split /:\s+/;
    }
    
    # Get bootstrapped data
    if (/Resampled distributions/) {$our_part_almost_started = 1};
    if (/dn\/ds/ && $our_part_almost_started) {$our_part_started = 1};
    
    next unless $our_part_started;

    my ($dnds, $pnps_freq, $pnps_all, $alpha_freq, $alpha_all) = split /\t/;
    
    push @dnds, $dnds;
    push @pnps_freq, $pnps_freq;
    push @pnps_all, $pnps_all;
    push @alpha_freq, $alpha_freq;
    push @alpha_all, $alpha_all;
}    


# define statistics objects
my $dnds_stat = Statistics::Descriptive::Full->new();
my $pnps_freq_stat = Statistics::Descriptive::Full->new();
my $pnps_all_stat = Statistics::Descriptive::Full->new();
my $alpha_freq_stat = Statistics::Descriptive::Full->new();
my $alpha_all_stat = Statistics::Descriptive::Full->new();


$alpha_freq_stat->add_data(@alpha_freq);

# get percentiles
my $alpha_freq_low_percent  = $alpha_freq_stat->percentile(2.5);
my $alpha_freq_high_percent = $alpha_freq_stat->percentile(97.5);


# output real data
print "\n";
print "$dnds_r\t";
print "$pnps_freq_r\t";
print "$pnps_all_r\t";
print "$alpha_freq_r\t";
print "$alpha_freq_low_percent\t";
print "$alpha_freq_high_percent\t";
print "$alpha_all_r\t";
