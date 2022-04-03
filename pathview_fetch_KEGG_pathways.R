rm(list=ls())

library(pathview)
library(org.Mm.eg.db)

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
		# node.info is a parser function
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

# column 1,2,3: EntrezGene ID, Ensembl Gene ID, KEGG Pathway ID
gene2pathway = cbind(agene,aensg,apath)

