rm(list=ls())
date()
cat("Start: \n")
library(biomaRt)
# https://github.com/grimbough/biomaRt/issues/31
httr::set_config(httr::config(ssl_verifypeer = FALSE))

i="hsapiens"
j="mmusculus"
ensembl_i = useMart("ensembl",dataset = paste(i,"_gene_ensembl",sep=""));
ams = getBM(attributes = c("ensembl_gene_id",paste(j,"_homolog_ensembl_gene",sep=""),paste(j,"_homolog_orthology_type",sep="")),mart=ensembl_i);
                                                              ams11 = ams[which(ams[,paste(j,"_homolog_orthology_type",sep="")] == "ortholog_one2one"),]

LRhuman = as.matrix(read.table("/mnt/extraids/OceanStor-0/linpei/genome/human_lr_pair.txt",header=T,sep="\t"))
LRmouse = as.matrix(read.table("/mnt/extraids/OceanStor-0/linpei/genome/mouse_lr_pair.txt",header=T,sep="\t"))

LRhuman_1 = LRhuman[which(!is.na(LRhuman[,"ligand_ensembl_gene_id"]) & !is.na(LRhuman[,"receptor_ensembl_gene_id"])),]
LRmouse_1 = LRmouse[which(!is.na(LRmouse[,"ligand_ensembl_gene_id"]) & !is.na(LRmouse[,"receptor_ensembl_gene_id"])),]

lihuman = LRhuman_1[,"ligand_ensembl_gene_id"]
rehuman = LRhuman_1[,"receptor_ensembl_gene_id"]

lihumanuniq = unique(lihuman)
rehumanuniq = unique(rehuman)

lihuman2mus = matrix(data=NA,ncol=2,nrow=length(lihumanuniq))
for(k in 1:length(lihumanuniq)) {
	j = lihumanuniq[k]
	if(sum(ams11[,"ensembl_gene_id"] == j) == 1) {
		i = which(ams11[,"ensembl_gene_id"] == j)
		lihuman2mus[k,2] = ams11[i,"mmusculus_homolog_ensembl_gene"]
		lihuman2mus[k,1] = j
}}

rehuman2mus = matrix(data=NA,ncol=2,nrow=length(rehumanuniq))
for(k in 1:length(rehumanuniq)) {
	j = rehumanuniq[k]
	if(sum(ams11[,"ensembl_gene_id"] == j) == 1) {
		i = which(ams11[,"ensembl_gene_id"] == j)
		rehuman2mus[k,2] = ams11[i,"mmusculus_homolog_ensembl_gene"]
		rehuman2mus[k,1] = j
}}

lihuman2mus = lihuman2mus[which(!is.na(lihuman2mus[,1])),]
rehuman2mus = rehuman2mus[which(!is.na(rehuman2mus[,1])),]

LRhuman_1_mus_1 = matrix(data=NA,ncol=2,nrow=nrow(LRhuman_1))
for(i in 1:nrow(LRhuman_1)) {
	xli = LRhuman_1[i,"ligand_ensembl_gene_id"]
	xre = LRhuman_1[i,"receptor_ensembl_gene_id"]
	if(sum(lihuman2mus[,1] == xli) == 1 & sum(rehuman2mus[,1] == xre) == 1) {
		LRhuman_1_mus_1[i,1] = lihuman2mus[which(lihuman2mus[,1] == xli),2]			
		LRhuman_1_mus_1[i,2] = rehuman2mus[which(rehuman2mus[,1] == xre),2]
	}
}



data = read.table("MatrixL001_geneXtile.o",sep="\t")

tiles = rep('t',168)
i = 1
tilenumeric = rep(1,168)
for(j in c(1,2)) {
        for(k in c(1,2,3,4,5,6)) {
                for(u in c("01","02","03","04","05","06","07","08","09","10","11","12","13","14")) {
                        ti = paste("T",j,k,u,sep="")
                        tiles[i] = ti
                        tilenumeric[i] = as.numeric(paste(j,k,u,sep=""))
                        i = i + 1
                }
        }
}
colnames(data) = c("L2Read","Gene",tiles)
nn = 1

library(doParallel)
registerDoParallel(32)

genegroup = list()
genegroup[["Ligand"]] = lihuman2mus
genegroup[["Receptor"]] = rehuman2mus

resultpgroup = list()

for(U in c("Ligand","Receptor")) {
	genehi = genegroup[[U]][,2]
	
	resultp = foreach(j=tilenumeric,.combine=rbind) %dopar% {
	        sumo = rep(0,nrow(data))
	        if(nn == 1) {
	                sumo = data[,paste("T",j,sep="")]
	                ntile = 1
	        } else {
	                tile9 = matrix(data=0,nrow=3,ncol=3)
	                tile9[2,2] = j
	                tile9[2,1] = j - 100
	                tile9[2,3] = j + 100
	                tile9[1,] = tile9[2,] - 1
	                tile9[3,] = tile9[2,] + 1
	                ntile = sum(paste("T",as.vector(tile9),sep="") %in% tiles)
	                if(ntile >= 2) {
	                        tilex = intersect(paste("T",as.vector(tile9),sep=""),tiles)
	                        sumo = apply(data[,tilex],1,sum)
	                }
	        }
	        if(sum(sumo) == 0) {
	        } else {
	                ngenehij = length(unique(data[which(data[,"Gene"] %in% genehi & sumo == 1),"Gene"]))
	                a = sum(data[,"Gene"] %in% genehi & sumo == 1)
	                b = sum(data[,"Gene"] == "ENSMUSG00000000000" & sumo == 1)
	                c = sum(data[,"Gene"] %in% genehi & sumo == 0)
	                d = sum(data[,"Gene"] == "ENSMUSG00000000000" & sumo == 0)
	                m2 = matrix(data=c(a,b,c,d),2,2)
	                if(sum(m2 == 0) <= 1) {
	                        m2test = chisq.test(m2)
	                        obshi = 0
	                        if(m2test$observed[1,1] > m2test$expected[1,1]) {
	                                obshi = 1
	                        }
	                        msg = c(j,ntile,ngenehij,m2test$observed[1,1],m2test$expected[1,1],b,c,d,m2test$p.value,obshi)
	                        msg
	                }
	        }
	}
	resultp = resultp[order(resultp[,9],decreasing=F),]
	rownames(resultp) = paste("Tile_",as.vector(resultp[,1]),sep="")
	resultpgroup[[U]] = resultp
}

resultp_1 = resultpgroup[["Ligand"]]
resultp_2 = resultpgroup[["Receptor"]]

x = intersect(rownames(resultp_1),rownames(resultp_2))

pdf("pval_Li_Re.pdf");
plot(x=-1*log10(resultp_1[x,9]),y=-1*log10(resultp_2[x,9]),pch=20,xlab="Ligand",ylab="Receptor");
dev.off()

