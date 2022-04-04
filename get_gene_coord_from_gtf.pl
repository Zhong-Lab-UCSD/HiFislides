#!/usr/bin/perl -w
use strict;

my %chrom;
foreach my $i (1..26) {
	$chrom{$i} = 1;
}
$chrom{"MT"} = 1;
$chrom{"X"} = 1;
$chrom{"Y"} = 1;

my $gtf = shift @ARGV; # the gene-centric gtf
my $geneid = shift @ARGV;

open IN,$gtf;
while(<IN>) {
	chomp;
	if(m/\tgene\t/) {
		my @a = split /\t/;
		if(m/($geneid\d+)/) {
			my $gene = $1;
			if(m/gene_biotype "(\S+)"/) {
				my $bioty = $1;
				my $chro;
				if(exists $chrom{$a[0]}) {
					$chro = $a[0];
				} else {
					$chro = 23;
				}
				my $chr = $a[0];
				print $chr,"\t",$a[3],"\t",$a[4],"\t",$gene,":",$chro,":",$a[3],":",$a[4],":",$a[6],":",$bioty,"\n";
			}
		}
	}
}
close IN;

