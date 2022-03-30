###
rm(list=ls())

date()
cat("Start: \n")

library(biomaRt)
library(reactome.db)
library(org.Mm.eg.db)
library(pathview)

#
#
#
# this obj is a list
reactomeGene2Pathways = as.list(reactomeEXTID2PATHID)
#### 
idmap = org.Mm.egENSEMBL

ensgtopathways_1 = array()
ensgtopathways_2 = array()
ensgtopathways_3 = array()

k = 1
for(j in sample(names(reactomeGene2Pathways),5)) {
	if(length(idmap[[j]]) > 0) {
		mypath = reactomeGene2Pathways[[j]]
		for(p in mypath) {
			ensgtopathways_1[k] = idmap[[j]]
			ensgtopathways_3[k] = paste(j,"_at",sep="")
			ensgtopathways_2[k] = p
			k = k + 1
		}
	}
}
############################################
#
# [1] "/mnt/extraids/OceanStor-0/linpei/genome"
LRhuman = as.matrix(read.table("/mnt/extraids/OceanStor-0/linpei/genome/human_lr_pair.txt",header=T,sep="\t"))
LRmouse = as.matrix(read.table("/mnt/extraids/OceanStor-0/linpei/genome/mouse_lr_pair.txt",header=T,sep="\t"))

LRhuman_1 = LRhuman[which(!is.na(LRhuman[,"ligand_ensembl_gene_id"]) & !is.na(LRhuman[,"receptor_ensembl_gene_id"])),]
LRmouse_1 = LRmouse[which(!is.na(LRmouse[,"ligand_ensembl_gene_id"]) & !is.na(LRmouse[,"receptor_ensembl_gene_id"])),]

lihuman = LRhuman_1[,"ligand_ensembl_gene_id"]
rehuman = LRhuman_1[,"receptor_ensembl_gene_id"]
i="hsapiens"
j="mmusculus"

ensembl_i = useMart("ensembl",dataset = paste(i,"_gene_ensembl",sep=""));

ams = getBM(attributes = c("ensembl_gene_id",paste(j,"_homolog_ensembl_gene",sep=""),paste(j,"_homolog_orthology_type",sep="")),mart=ensembl_i);
                                                              ams11 = ams[which(ams[,paste(j,"_homolog_orthology_type",sep="")] == "ortholog_one2one"),]

lihuman2mus_1 = array()
lihuman2mus_2 = array()
k = 1
for(j in unique(lihuman)) {
	if(sum(ams11[,"ensembl_gene_id"] == j) > 0) {
		for(i in which(ams11[,"ensembl_gene_id"] == j)) {
			lihuman2mus_2[k] = ams11[i,"mmusculus_homolog_ensembl_gene"]
			lihuman2mus_1[k] = j
			k = k + 1
		}
	}
}

rehuman2mus_1 = array()
rehuman2mus_2 = array()
k = 1
for(j in unique(rehuman)) {
	if(sum(ams11[,"ensembl_gene_id"] == j) > 0) {
		for(i in which(ams11[,"ensembl_gene_id"] == j)) {
			rehuman2mus_2[k] = ams11[i,"mmusculus_homolog_ensembl_gene"]
			rehuman2mus_1[k] = j
			k = k + 1
		}
	}
}

LRhuman_1_mus_1 = array()
LRhuman_1_mus_2 = array()

for(i in 1:nrow(LRhuman_1)) {
	xli = LRhuman_1[i,"ligand_ensembl_gene_id"]
	xre = LRhuman_1[i,"receptor_ensembl_gene_id"]
	if(sum(lihuman2mus_1 == xli) == 1 & sum(rehuman2mus_1 == xre) == 1) {
		LRhuman_1_mus_1[i] = lihuman2mus_2[which(lihuman2mus_1 == xli)]			
		LRhuman_1_mus_2[i] = rehuman2mus_2[which(rehuman2mus_1 == xre)]
	}
}

############################################
data(paths.hsa)
paths.mmu = gsub("hsa","",names(paths.hsa))

idmap = org.Mm.egENSEMBL

agene = array()
aensg = array()
apath = array()
k = 1

for(pathid in sample(paths.mmu,5)) {
	msg = download.kegg(pathway.id = pathid, species = "mmu", kegg.dir = ".",file.type="xml")
	if(msg == "succeed") {
		node.data=node.info(paste("mmu",pathid,".xml",sep=""))	
		for(j in names(node.data$kegg.names)) {
			genepernode = node.data$kegg.names[[j]]
			genepernode = intersect(genepernode,mappedkeys(idmap))
			for(g in genepernode) {
				agene[k] = paste(g,"_at",sep="")
				aensg[k] = idmap[[g]]
				apath[k] = pathid
				k = k + 1
			}
		}
	}
}

gene2pathway = cbind(agene,aensg,apath)


