#!/usr/bin/perl -w
use strict;

my $datestring = localtime();
print STDERR $datestring,"\n";

my @sam = @ARGV;

my %hifi2as;
my %hifi2highest;
my %hifi2spot0;
my %mapped_spot;

my %mappedbc_per_tile;

foreach my $sam (@sam) {
	open IN,$sam;
	while(<IN>) {
		chomp;
		my $go = 0;
		my $sum_of_flag;
		if(m/\t0\t\S+:/) {
			$go = 1;
			$sum_of_flag = 0;
		}
		if(m/\t256\t\S+:/) {
			$go = 1;
			# in cases that 256 indicate a supplementary/secondary alignment but with alignment score tie with the 1st one.
			$sum_of_flag = 256;
		}
		# 0 and 256 are +/+ aligned.
		if($go == 1) {
			if(m/AS:i:(\d+)/){
				my $as = $1;
				my @a = split /\t/;
				my $spot = $a[2];
				my $hifi = $a[0];

				$hifi2highest{$hifi} = 0 unless exists $hifi2highest{$hifi};

				#a[3]	POS	1-based leftmost POSition/coordinate of clipped sequence
				#a[4]	MAPQ	MAPping Quality (Phred-scaled)
				#a[5]	CIAGR	extended CIGAR string

				# a[3] is the 1-based offset into the forward reference strand where leftmost character of the alignment occurs
				# a[4] is the Mapping quality. It is said that 
				# a[5] is the CIGAR string representation of alignment
				
				$hifi2as{$hifi}->{$as}->{$spot} = 1;
				if($as > $hifi2highest{$hifi}) {
					$hifi2highest{$hifi} = $as;
				}
			}
		}
	}
	print STDERR $sam,"\n";
	close IN;
}

print STDERR "Total_Number_of_HiFi_reads_aligned_with_barcode","\t";
print STDERR scalar keys %hifi2as,"\n";

my $n1 = 0;
foreach my $hifi (keys %hifi2highest) {
	my $as = $hifi2highest{$hifi};
	my @spot = keys %{$hifi2as{$hifi}->{$as}};
	my $N = 0;

	foreach my $spot (@spot) {
		$spot=~m/_(\d+)$/;
		my $n = $1;
		# n > 1 when multiple L1R1 read share the same sequence.
		$N = $N + $n;
	}

	# barcode(s) with highest alignment score mapped by each HIFI read
	# when N = 1 and scalar @spot == 1, spot[0] is the only one barcode exactly mapped by this HIFI read with the highest alignment score.

	# print $hifi,"\t",$spot[0],"\t",scalar @spot,"\t",$N,"\t",$as,"\n";
	
	if($N <= 1000) {
		foreach my $spot (@spot) {
			print $hifi,"\t",$spot,"\t",scalar @spot,"\t",$N,"\t",$as,"\n";	
		}
	}
}

$datestring = localtime();
print STDERR $datestring,"\n";
