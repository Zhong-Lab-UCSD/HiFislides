library(grDevices)

args = commandArgs()
filein1 = args[which(args == "--args") + 1]
filein2 = args[which(args == "--args") + 2]
fileop = args[which(args == "--args") + 3]

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
topsurface = tilenumeric[which(tilenumeric < 2000)]
botsurface = tilenumeric[which(tilenumeric > 2000)]

twosurface = list()
twosurface[["top"]] = topsurface
twosurface[["bot"]] = botsurface

# pvall = c(ans[["hi500"]][,4],ans[["hi1000"]][,4],ans[["hi2000"]][,4])
mains = c("top-500 genes","top-1000","top-2000")
names(mains) = c("hi500","hi1000","hi2000")

# pvall = pvall[order(pvall,decreasing=F)]
# pvallqt = quantile(pvall,probs=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))

#cololist = list()
#cololist[["hi"]] = hcl.colors(5,palette = "Reds 3")[2:5]
#cololist[["lo"]] = hcl.colors(5,palette = "Blues 3")[2:5]

# a resource for color:
# https://developer.r-project.org/Blog/public/2019/04/01/hcl-based-color-palettes-in-grdevices/

cololist=list()
cololist[["hi"]] = c("#CC1C2F","#FF7078","#FFBEC1")
cololist[["lo"]] = c("#0072B4","#79ABE2","#C3DBFD")

colegend = c("#CC1C2F","#FF7078","#FFBEC1","#0072B4","#79ABE2","#C3DBFD","grey","black")

for(i in 1:8) {
	cat(i,colegend[i],"\n")
}

pdf(fileop)
for(gene in names(mains)) {
	d_surf = 0
	plot(axes=F,main=mains[gene],1:40,1:40,col="white",pch=20,xlab="6 Swaths per surface,2 surfaces per lane",ylab="")
	text(y=8,x=10,"Lane 1 \n top")
	text(y=8,x=18,"Lane 1 \n bottom")
	text(y=8,x=26,"Lane 2 \n top")
	text(y=8,x=34,"Lane 2 \n bottom")
	legend(x=5,y=40,pch=c(20,20,20,20,20,20,1,4),
		col=colegend,
		legend=c("< 1e-50","> 1e-50 & < 1e-10","> 1e-10 & < 0.05","< 1e-50","> 1e-50 & < 1e-10","> 1e-10 & < 0.05","> 0.05","missing tile"))
	for(filein in c(filein1,filein2)) {
		load(filein)
		for(surf in c("top","bot")) {
			surface = twosurface[[surf]]
			ansm = ans[[gene]]
			tilm = matrix(data=surface,ncol=6,nrow=14)
			for(j in 1:6) {
				for(i in 1:14) {
					cologrp = "hi";
					if(sum(ansm[,1] == tilm[i,j]) == 1) {
					} else {
						points(x=j+6.5+d_surf,y=i+10,pch=4)
						next;
					}	
					if(ansm[which(ansm[,1] == tilm[i,j]),2] > ansm[which(ansm[,1] == tilm[i,j]),3]) {
						cologrp = "hi"
					}
					if(ansm[which(ansm[,1] == tilm[i,j]),2] < ansm[which(ansm[,1] == tilm[i,j]),3]) {
						cologrp = "lo"
					}
					colopvbin = cololist[[cologrp]]
			
					pval = ansm[which(ansm[,1] == tilm[i,j]),4]
					mycolo = "grey";
					if(pval < 1e-50) {
						mycolo = colopvbin[1]
					} else {
						if(pval > 1e-50 & pval < 1e-10) {
							mycolo = colopvbin[2]
						}
						if(pval > 1e-10 & pval < 0.05) {
							mycolo = colopvbin[3]
						}
					}
					mydot = 20
					if(pval > 0.05) {
						mycolo = "grey"
						mydot = 1
					}
					points(x=j+6.5+d_surf,y=i+10,pch=mydot,col=mycolo)
				}
			}
			d_surf = d_surf + 8
		}
	}
}
dev.off()
