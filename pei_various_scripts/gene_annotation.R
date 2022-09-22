# Reference:
# https://gist.github.com/slowkow/395448fe1094c8ea4c7d

library(biomaRt)

ensembl_ids = c("ENSG00000243485", "ENSG00000237613", "ENSG00000186092", "ENSG00000238009")

mart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

dat = getBM(values = ensembl_ids,filters = c("ensembl_gene_id"),attributes = c("ensembl_gene_id", "external_gene_name"),mart = mart)

# Reference:
# https://www.biostars.org/p/178726/

library(org.Hs.eg.db)
annot = select(org.Hs.eg.db,keys = keys(org.Hs.eg.db),columns = c('ENTREZID','SYMBOL','ENSEMBL'),keytype = 'ENTREZID')
  
