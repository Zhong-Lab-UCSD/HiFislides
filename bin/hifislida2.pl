#!/usr/bin/perl -w
use strict;


my $datestring = localtime();
print STDERR $datestring,"\n";

my $hifi_N_file;
$hifi_N_file = shift @ARGV;

my %hifi_N;
open IN,$hifi_N_file;
while(<IN>) {
	chomp;
	if(m/^MN00185/) {
		my @a = split /\t/;
		my ($rd,$spot,$NSpot,$N,$as) = @a;
		if($NSpot == 1 && $N == 1 && $spot=~m/_1$/) {
			$hifi_N{$rd} = $spot;
		}
		# NSpot spots aligned by this hifi read tied with the highest AS
		# this hifi_N_file was outputed by hifia.pl
	}
}
close IN;

# SAM per tile
my $sam = shift @ARGV;

my %read_to_kept;
my $read_to_kept = "NAN";
$read_to_kept = shift @ARGV if @ARGV > 0;
if($read_to_kept ne "NAN") {
	open IN,$read_to_kept;
	while(<IN>) {
		chomp;
		if(m/(MN00185:\d+:\S+:\d:\d+:\d+:\d+)/) {
			my $rd = $1;
			$read_to_kept{$rd} = 1;
		}
	}
	close IN;
} else {
	%read_to_kept = %hifi_N;
}
my %spotcount_per_tile;
my %readcount_per_tile;

my %read_per_tile;
my %mapped_spot_per_tile;

foreach my $samI ($sam) {
	open IN,$samI;
	while(<IN>) {
		chomp;
		my $go = 0;
		if(m/\t0\tVH00454:/) {
			$go = 1;
		}
		if(m/\t256\tVH00454:/) {
			$go = 1;
		}
		if($go == 1) {
			my @a = split /\t/;
			my $spot = $a[2];
			my $hifi = $a[0];
			if(exists $read_to_kept{$hifi}) {
				if(exists $hifi_N{$hifi}) {
					#if($NSpot == 1 && $N == 1 && $spot=~m/_1$/) {
					# $hifi_N{$rd} = $spot;
					if($spot eq $hifi_N{$hifi}) {						
						if($spot=~m/:1:(1\d\d\d):(\d+):(\d+)/) {
							my ($TILE,$x,$y) = ($1,$2,$3);
							$TILE = "T".$TILE;
							$mapped_spot_per_tile{$TILE}->{$spot} = 1;							
							$read_per_tile{$TILE}->{$hifi} = 1;
						}
					}
				}
			}
		}
	}
	print STDERR $sam,"\n";
	close IN;
}

my @tile = keys %read_per_tile;

foreach my $tile (@tile) {
	my $nhifi = scalar keys %{$read_per_tile{$tile}};
	$readcount_per_tile{$tile} = $nhifi;
	my $nspot = scalar keys %{$mapped_spot_per_tile{$tile}};
	$spotcount_per_tile{$tile} = $nspot;
}

@tile = sort {$readcount_per_tile{$b} <=> $readcount_per_tile{$a}} @tile;

foreach my $tile (@tile) {
	print ">TILE","\t",$tile,"\t",$readcount_per_tile{$tile},"\t",$spotcount_per_tile{$tile},"\n";
}

$datestring = localtime();
print STDERR $datestring,"\n";

