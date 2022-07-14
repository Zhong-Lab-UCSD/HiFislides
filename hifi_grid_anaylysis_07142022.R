###
##
#

#
##
###
data1 = read.table("Output_spot_to_gene_for_in.txt",sep="\t")
data2 = read.table("Output_spot_to_gene_for_ex.txt",sep="\t")
x = 1114:1101
tilem = cbind(x,x+100,x+200,x+300,x+400,x+500)
plot(1:10,1:10,xlim=c(0,80000*12),ylim=c(0,80000*12),pch=20,col="white")
colo = c("purple","red","green3")
for(i in 1:11) {
for(j in 1:6) {
	# the first one is 1101.
	Tile = tilem[i,j]
	dx = (j - 1) * 80000
	dy = (i - 1) * 60000
	data1i = data1[which(data1[,1] == Tile),]
	data2i = data2[which(data2[,1] == Tile),]
	colo = c("purple","red","green3")
	x0 = 1000
	y0 = 1000
	# d: controls the size of grid
	while(x0 < 80000) {
		x1 = x0 + d
		while(y0 < 70000) {
			y1 = y0 + d
			# cat(x0,x1,y0,y1,"\n")
			n1 = sum(data1i[,2] > x0 & data1i[,2] < x1 & data1i[,3] > y0 & data1i[,3] < y1)
			n2 = sum(data2i[,2] > x0 & data2i[,2] < x1 & data2i[,3] > y0 & data2i[,3] < y1)
			if(n1 == 0 & n2 == 0) {
				rect(xleft=x0+dx, ybottom=y0+dy, xright=x1+dx, ytop=y1+dy,col="white",lwd=0.1,border="white")
			} else {
				if(n1 > n2) {
					rect(xleft=x0+dx, ybottom=y0+dy, xright=x1+dx, ytop=y1+dy,col=colo[1],lwd=0.1,border="white")
				} else {
					if(n2 > n1) {
						rect(xleft=x0+dx, ybottom=y0+dy, xright=x1+dx, ytop=y1+dy,col=colo[2],lwd=0.1,border="white")
					} else {
						rect(xleft=x0+dx, ybottom=y0+dy, xright=x1+dx, ytop=y1+dy,col="grey",lwd=0.1,border="white")
					}
				}
			}
			y0 = y1	
		}
		x0 = x1
		y0 = 1000
	}
}
}


