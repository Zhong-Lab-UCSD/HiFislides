library(sf)
library(dplyr)
library(mapview)
library(tmap)

args = commandArgs()
filt = args[which(args == '--args')+1]
nbin = args[which(args == '--args')+2]
nbin = as.numeric(nbin)

data1 = read.table(gsub("=",filt,"/mnt/extraids/OceanStor-0/linpei/hifi/data_3/lib2/m2to1_using_all_=_1.o"),sep="\t")
data2 = read.table(gsub("=",filt,"/mnt/extraids/OceanStor-0/linpei/hifi/data_3/lib2/m2to1_using_all_=_2.o"),sep="\t")

tilebot_1 = 21101:21110
tilebot_2 = 22101:22110
tilebot_3 = 23101:23110

tilebot = rbind(tilebot_1,tilebot_2,tilebot_3)

data1_new_coord = matrix(data=0,ncol=2,nrow=nrow(data1))
data2_new_coord = matrix(data=0,ncol=2,nrow=nrow(data2))
colnames(data1_new_coord) = c('myX','myY')
colnames(data2_new_coord) = c('myX','myY')

rowk1 = 1
rowk2 = 1

dy = 0
for(myi in 1:3) {	
	dx = 0
	for(myj in 1:10) {

		ti = tilebot[myi,myj]
		ti = paste('T',ti,sep='')
		if(sum(data1[,1] == ti) > 0) {
			data1i = data1[which(data1[,1] == ti),]
			data2i = data2[which(data2[,1] == ti),]
			for(i in 1:nrow(data1i)) {
				# points(y=data_i[i,4]+dy,x=data_i[i,5]+dx,pch=20,cex=2,col=colo[i])
				data1_new_coord[rowk1,'myY'] = data1i[i,4] + dy
				data1_new_coord[rowk1,'myX'] = data1i[i,5] + dx
				rowk1 = rowk1 + 1
			}
			for(i in 1:nrow(data2i)) {
				# points(y=data_i[i,4]+dy,x=data_i[i,5]+dx,pch=20,cex=2,col=colo[i])
				data2_new_coord[rowk2,'myY'] = data2i[i,4] + dy
				data2_new_coord[rowk2,'myX'] = data2i[i,5] + dx
				rowk2 = rowk2 + 1
			}
		}
		# go to the next tile
		dx = dx + 20000
	}
	dy = dy + 25000
}	

matr_1 = data1_new_coord
colnames(matr_1) = c('myX','myY')
matr_1 = as.data.frame(matr_1)

matr_2 = data2_new_coord
colnames(matr_2) = c('myX','myY')
matr_2 = as.data.frame(matr_2)

test_points = st_as_sf(matr_1,coords=c('myX','myY'))
test_points_st_make_grid = st_make_grid(test_points,n=c(nbin,nbin),what="polygons",square = FALSE)
test_points_st_make_grid_st_sf = st_sf(test_points_st_make_grid)
test_points_st_make_grid_st_sf$nspot = lengths(st_intersects(test_points_st_make_grid_st_sf,test_points))

test_reads = st_as_sf(matr_2,coords=c('myX','myY'))
test_points_st_make_grid_st_sf$num_of_read = lengths(st_intersects(test_points_st_make_grid_st_sf,test_reads))
test_points_st_make_grid_st_sf$num_of_read_per_spot = test_points_st_make_grid_st_sf$num_of_read/test_points_st_make_grid_st_sf$nspot
test_points_st_make_grid_st_sf_filter = filter(test_points_st_make_grid_st_sf,num_of_read>0)

tmap_mode("plot")
map_honeycomb2 = tm_shape(test_points_st_make_grid_st_sf_filter) + 
tm_fill(
col = "num_of_read",
palette = "Reds",
style = "cont",
alpha=0.6)
# map_honeycomb2
tmap_save(map_honeycomb2,gsub("=",filt,"BotSurface_num_of_=_per_hexagonal_bin.png"))

# tmap_mode("plot")
# map_honeycomb = tm_shape(test_points_st_make_grid_st_sf_filter) + 
# tm_fill(
# col = "num_of_read",
# palette = "Reds",
# style = "cont",
# alpha=0.6)
#  tmap_save(map_honeycomb,gsub("=",filt,"BotSurface_num_of_=.png"))

#pdf("hexa_nspot.pdf")
#map_honeycomb3 = tm_shape(test_points_st_make_grid_st_sf) + 
#tm_fill(
#col = "nspot",
#palette = "Reds",
#style = "cont",
#alpha=0.6)
#map_honeycomb3
#dev.off()
