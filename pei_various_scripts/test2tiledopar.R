library(doParallel)
library(plyr)
# > detectCores()
# [1] 192

args = commandArgs()

filein = args[which(args == "--args") + 1]

# nn: num of neighboring tiles
nn = as.numeric(args[which(args == "--args") + 2])

ntopgene = as.numeric(args[which(args == "--args") + 3])

# nt: num of threads
nt = as.numeric(args[which(args == "--args") + 4]);registerDoParallel(nt)

fileou = paste(filein,"_",nn,"_",ntopgene,".txt",sep="")

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

data = read.table(filein,sep="\t");
colnames(data) = c("L2Read","Gene",tiles)

generdcnt = count(data[,"Gene"])
generdcnt = generdcnt[order(generdcnt[,'freq'],decreasing=T),]

ntopgene = ntopgene + 3
genehi = setdiff(generdcnt[1:ntopgene,"x"],c("ENSMUSG00000000000","ENSMUSGNNNNNNNNNNN","ENSMUSG00000099364"))

# genehi_500 = setdiff(generdcnt[1:503,"x"],c("ENSMUSG00000000000","ENSMUSGNNNNNNNNNNN","ENSMUSG00000099364"))
# genehi1000 = setdiff(generdcnt[1:1003,"x"],c("ENSMUSG00000000000","ENSMUSGNNNNNNNNNNN","ENSMUSG00000099364"))
# genehi2000 = setdiff(generdcnt[1:2003,"x"],c("ENSMUSG00000000000","ENSMUSGNNNNNNNNNNN","ENSMUSG00000099364"))

# ans = list()
# result_tile2grid = matrix(data=0,nrow=length(tilenumeric),ncol=9)
# result_tile2grid[,8] = 1
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
# result_tile2grid = result_tile2grid[order(result_tile2grid[,8],decreasing=F),]
# ans[[paste(x,"_N",nn,sep="")]] = result_tile2grid
colnames(resultp) = c("Tile","GroupSize","NGene","Obs","Exp","x21","x12","x22","pval","Is_It_Overrepresented")
write.table(resultp,file=fileou,sep="\t",row.names=F,quote=F)
