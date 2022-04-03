rm(list=ls())

library(pathview)
data(paths.hsa)
paths.mmu = gsub("hsa","",names(paths.hsa))

library(org.Mm.eg.db)

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

