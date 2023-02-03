data = read.table("Output_by_hifislida2.o",sep="\t")

x = 1101:1111


tilm = cbind(x,x+100,x+200,x+300,x+400,x+500)

# tilm is the matrix for arranging Tile IDs as their relative position on the surface
# we have 11 x 6 tiles.

tile = as.numeric(gsub("T","",data[,2]))

data2 = cbind(tile,data[,2:4])

mat = matrix(data=NA,nrow=1,ncol=6)

# w: the size of ROI. 
# w = 3 indicate ROI has 3 x 3 tiles

# we do not know the exact shape and location of ROI.
# so, we rank potential ROIs to see which is most significantly different from the remaining tiles.

for(w in c(2,3,4)) {
	d = w - 1
	for(i in 1:(6-d)) {
		for(j in 1:(11-d)) {

			w1 = which(data2[,1] %in% tilm[j:(j+d),i:(i+d)])
			w0 = which(!(data2[,1] %in% tilm[j:(j+d),i:(i+d)]))
			pv = ks.test(data2[w1,3],data2[w0,3])$p.value
			# print(tilm[j:(j+2),i:(i+2)])
		
			mat = rbind(mat,c(d,i,j,mean(data2[w1,3]),mean(data2[w0,3]),pv))
		}
	}
}


for(i in 1:nrow(mat)) {
	if(mat[i,4] > mat[i,5] & mat[i,6] < 0.05) {
		i = mat[i,2]
		j = mat[i,3]
		d = mat[i,1]
		# print(tilm[j:(j+d),i:(i+d)])

		roi = tilm[j:(j+d),i:(i+d)]
		write.table(file="Tiles_in_ROI.txt",col.names=F,row.names=F,sep="",quote=F,as.vector(unlist(roi)))
		break
	}
}

save(mat,file="ks_test_to_identify_ROI_mat.RData")
