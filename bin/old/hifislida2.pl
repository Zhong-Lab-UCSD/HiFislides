#!/usr/bin/perl -w
use strict;


my $datestring = localtime();
print STDERR $datestring,"\n";

my $hifi_N_file;
$hifi_N_file = shift @ARGV;

my %hifi_N;
open IN,$hifi_N_file;

my %read_per_tile;
my %mapped_spot_per_tile;

while(<IN>) {
	chomp;
	if(m/^\S+:\d+:\S+:\d+:\d+:\d+:\d+/) {
		my @a = split /\t/;
		my ($rd,$spot,$NSpot,$N,$as) = @a;
		if($NSpot == 1 && $N == 1 && $spot=~m/_1$/) {
			$hifi_N{$rd} = $spot;
			if($spot=~m/:1:(\d+):(\d+):(\d+)/) {
				my ($TILE,$x,$y) = ($1,$2,$3);
				# $TILE = "T".$TILE;
				$read_per_tile{$TILE}->{$rd} = 1;
				$mapped_spot_per_tile{$TILE}->{$spot} = 1;
			}
		}
		# NSpot spots aligned by this hifi read tied with the highest AS
		# this hifi_N_file was outputed by hifia.pl
	}
}
close IN;

# SAM per tile

my @tile = keys %read_per_tile;
my %readcount_per_tile;
my %spotcount_per_tile;

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

