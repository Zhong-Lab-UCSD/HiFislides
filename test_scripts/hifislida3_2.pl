#!/usr/bin/perl -w
use strict;

my %hifi2as;
my %hifi2highest;
my %tile4roi;
my @tile4roi;


my $datestring = localtime();
print STDERR $datestring,"\n";
my ($aout,$a2out,$dupspots) = @ARGV;

open IN,$a2out;
while(<IN>) {
	chomp;
	if(m/(\d{4,5})/) {
		my $ti = $1;
		$tile4roi{$ti} = {};
	}
}
close IN;
my %dupspots;
open IN,$dupspots;
while(<IN>) {
	chomp;
	my @a = split /\t/;
	if(defined $a[1]) {
		$dupspots{$a[0]}->{$a[1]} = 1;	
	}
}
close IN;

my %spot_in_roi;
open IN,$aout;
while(<IN>) {
	chomp;
	# print $hifi,"\t",$spot[0],"\t",scalar @spot,"\t",$N,"\t",$as,"\n";
	my ($hifi,$spot,$NSpot,$N,$as) = split /\t/;
	if($N < 1000) {
		$spot=~s/_\d+//;
		if($spot=~m/:1:(\d{4,5}):/) {
			my $tile = $1;
			#my $Tile = "T".$tile;
			if(exists $tile4roi{$tile}) {
				$spot_in_roi{$hifi}->{$spot} = 1;
			}
		}
		my @dupspots = keys %{$dupspots{$spot}};
		foreach my $spot_1 (@dupspots) {
			if($spot_1=~m/:1:(\d{4,5}):/) {
				my $tile = $1;
				#my $Tile = "T".$tile;
				if(exists $tile4roi{$tile}) {
					$spot_in_roi{$hifi}->{$spot_1} = 1;
				}
			}
		}
	}
}
close IN;

print "HiFi_read_id","\t","tile_id","\t","col","\t","row","\t","N","\n";

foreach my $hifi (keys %spot_in_roi) {
	my @spot = keys %{$spot_in_roi{$hifi}};
	my $N = scalar @spot;
	foreach my $spot (@spot) {
		if($spot=~m/:1:(\d{4,5}):(\d+):(\d+)/) {
			# range of y is larger than x.
			# Y: rows, 119 grids.
			# X: cols, 85 grids.
			my ($tile,$y,$x) = ($1,$2,$3);
			print $hifi,"\t",$tile,"\t",$x,"\t",$y,"\t",$N,"\n";
		}
	}
}

print STDERR "Number of spatially resolved HiFi reads within ROI","\n";
print STDERR scalar keys %spot_in_roi,"\n";

$datestring = localtime();
print STDERR $datestring,"\n";
