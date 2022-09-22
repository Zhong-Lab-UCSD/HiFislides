library(ggplot2)

 data = read.table("mappedbarcode_per_tile.o",sep=" ")
 tilestop = matrix(data=0,ncol=3,nrow=10)
 tilestop[,1] = 11110:11101
 tilestop[,2] = 12110:12101
 tilestop[,3] = 13110:13101
 mappedbarcodetop = matrix(data=0,ncol=3,nrow=10)
 for(i in 1:10) {
 for(j in 1:3) {
 mappedbarcodetop[i,j] = data[which(data[,1] == tilestop[i,j]),2]
 }
 }
 
 tilesbot = matrix(data=0,ncol=3,nrow=10)
 tilesbot[,1] = 21110:21101
 tilesbot[,2] = 22110:22101
 tilesbot[,3] = 23110:23101
 mappedbarcodebot = matrix(data=0,ncol=3,nrow=10)
 for(i in 1:10) {
 for(j in 1:3) {
 mappedbarcodebot[i,j] = data[which(data[,1] == tilesbot[i,j]),2]
 }
 }

mappedbc = cbind(mappedbarcodetop,mappedbarcodebot)

vx = array()
vy = array()
mapped_barcode = array()
k = 1
for(i in 1:nrow(mappedbc)) {
for(j in 1:ncol(mappedbc)) {
vx[k] = i
vy[k] = j
mapped_barcode[k] = mappedbc[i,j]
k = k + 1
}
}
df = data.frame(vx,vy,mapped_barcode)
pdf("fig.pdf")
ggplot(df,aes(vy,vx)) + geom_point(aes(size = 20,colour = mapped_barcode)) + scale_colour_gradient2() + 
theme(axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()) + xlab("") + ylab("")
dev.off()

