#!/usr/bin/perl -w
use strict;

my %ens;
open IN,"/mnt/extraids/OceanStor-0/linpei/genome/gene/gene2ensembl";
while(<IN>) {
	chomp;
	# entrezgene to ensg
	# both human and mouse.
	my $ensg = 'NA';
	if(m/(ENSG\d+)/) {
		$ensg = $1;
	}
	if(m/(ENSMUSG\d+)/) {
		$ensg = $1;
	}
	if($ensg=~m/\d+/) {
		my @a = split /\t/;
		my $entrezg = $a[1];
		$ens{$entrezg."_at"} = $ensg;
	}
}
close IN;

open IN,shift @ARGV;
# CORUM core data
<IN>;
while(<IN>) {
	chomp;
	my @a = split /\t/;

	my $cplx = $a[1];
	$cplx=~s/\s/_/g;

	my $geneid = $a[6];
	my @geneid = split /;/,$geneid;
	foreach my $gene (@geneid) {
		my $entrezg = $gene."_at";
		if(exists $ens{$entrezg}) {
			my $ensg = $ens{$entrezg};
			my $species = "xxx";
			$species = "mmu" if $ensg=~m/ENSMUSG/;
			$species = "hsa" if $ensg=~m/ENSG/;
			print $cplx,"\t",$species,"\t",$ensg,"\n";
		}
	}
}
close IN;
