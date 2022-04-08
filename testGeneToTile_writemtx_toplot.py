#!/usr/bin/python3.6

import sys
import numpy as np
import pandas as pd

RColor10 = ["burlywood4","aquamarine4","black","blue3","brown3","chocolate","red","darkorchid3","gold","orange"]

read2gene2coord = sys.argv[1]
tile = sys.argv[2]

mycolnames = ["Lane","Rowi","L2Read","Gene","Tile","L1Coord"]

mat = pd.read_csv(read2gene2coord,sep="\t",names=mycolnames,skiprows=1)

gene2spot = dict()
spot2gene = {}

genex = {}
genex["ENSMUSG00000099364"] = 1;
genex["ENSMUSG00000000000"] = 1;
genex["ENSMUSGNNNNNNNNNNN"] = 1;

for i, rowi in mat.iterrows():
	if rowi["Tile"] == tile:
		spot = rowi["L1Coord"];
		gene = rowi["Gene"];
		if gene in genex:
			next
		else:
			if spot in spot2gene:
				spot2gene[spot][gene] = 1;
			else:
				spot2gene[spot] = {}
				spot2gene[spot][gene] = 1;
			if gene in gene2spot:
				gene2spot[gene] = gene2spot[gene] + 1
			else:
				gene2spot[gene] = 1

nspot = sorted(gene2spot.values(),reverse=True)
gene2color = {}
for gene in gene2spot:
	if gene2spot[gene] > nspot[9]:
		if gene not in gene2color:
			if len(RColor10) > 0:
				gene2color[gene] = RColor10.pop()
for ns in nspot:
	for gene in gene2color:
		if gene2spot[gene] == ns:
			if gene in gene2color:
				print("# ",gene,gene2spot[gene],gene2color[gene])

gene2spot["NA"] = 0;

for i, rowi in mat.iterrows():
	if rowi["Tile"] == tile:
		spot = rowi["L1Coord"];
		gene = rowi["Gene"];
		if gene in genex:
			next
		else:
			L1coord = spot.split(':');
			if len(spot2gene[spot]) == 1:
				colo = "grey"
				if gene in gene2color:
					colo = gene2color[gene]
				print(gene,tile,L1coord[5],L1coord[6],colo)
			else:
				gene_1 = "NA";
				for gene in spot2gene[spot]:
					if gene2spot[gene] > gene2spot[gene_1]:
						gene_1 = gene;		 
				colo = "grey"
				if gene_1 in gene2color:
					colo = gene2color[gene_1]
				print(gene_1,tile,L1coord[5],L1coord[6],colo)
