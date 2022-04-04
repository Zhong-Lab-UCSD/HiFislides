#!/usr/bin/perl -w
use strict;

my @fq = @ARGV;

my %hd;
my %rd;

foreach my $fq (@fq) {

	if($fq=~m/\.gz$/) {
		open(IN,"gunzip -c $fq |");
	} else {
		open IN,$fq;
	}

	while(<IN>) {
		chomp;
		if(m/^\S(\w+:\d+:\w+:\d+:\d+:\d+:\d+)/) {
			# machine,i,flowcell,lane,tile,x,y
			my $hd = $1;
			my $rd = <IN>;
			chomp $rd;
			$rd{$rd}->{$hd} = 1;
		}
	}
	close IN;
	print STDERR $fq,"\n";
}

foreach my $rd (keys %rd) {
	my @hd = keys %{$rd{$rd}};
	my @coord;
	foreach my $hd (@hd) {
		$hd=~m/:(\d+:\d+:\d+:\d+)$/;
		my $coord = $1;
		push @coord,'L'.$coord;
	}
	my $coord = join "_",@coord;
	print ">",$hd[0],"|",scalar @hd,"|",$coord,"\n";
	print $rd,"\n";	
}
