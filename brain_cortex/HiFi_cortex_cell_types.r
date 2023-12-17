# TO BE RUN FIRST: /mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/result/HiFi_global_functions.r
# TO BE RUN SECOND: /mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/result/cortex/HiFi_cortex_init.r

##### Setup directories
OUT_DIR = paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/result/cortex/", SAMPLE_NAME, "/", ROI_label, "/cell_types/")
DATA_DIR=paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/IGM/", SAMPLE_NAME)
ROI_tiles = read.table(paste0(DATA_DIR, "/", ROI_label, "_tiles.txt"))$V1
dir.create(OUT_DIR, recursive = T, showWarnings = F)

##### Input cell type markers
df_cell_type_markers = read.xlsx("/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/result/cortex/human_protein_atlas_markers.xlsx", sheet = "Sheet1")
colnames(df_cell_type_markers)[1] = "gene"

# Table showing the number of marker genes for each cell type
df = data.frame(V1=names(table(as.factor(df_cell_type_markers$cell_type))),
                V2=as.numeric(table(as.factor(df_cell_type_markers$cell_type))))
write.xlsx(df, paste0(paste(strsplit(OUT_DIR, "/")[[1]][1:(length(strsplit(OUT_DIR, "/")[[1]])-3)], collapse = "/"), "/cell_type_markers_table.xlsx"))


##### Assign regions to spots

# Data frame with spot to region associations (each row is a spot)
# spot expression level: sum of the expression level of marker genes expressed in that spot
spot2cell_roi = assign_feature_to_spots(input_data,
                                        df_markers=df_cell_type_markers,
                                        feature_label="cell_type", # cell_type or layer
                                        with_mixed_features=F)
rownames(spot2cell_roi) = spot2cell_roi$spot_id
saveRDS(spot2cell_roi, paste0(OUT_DIR, "/spot2cell_roi.rds"))

df1 = aggregate(spot2cell_roi$spot_id, by = list(spot2cell_roi$cell_type), FUN=length)
df1$perc_spots = df1$x / nrow(spot2cell_roi)
colnames(df1)[1:2] = c("cell_type", "n_spots")
df1 = rbind(df1, data.frame(cell_type="Total", n_spots=nrow(spot2cell_roi), perc_spots = 1))

df2 = aggregate(spot2cell_roi$spot_id, by = list(spot2cell_roi$tile_id, spot2cell_roi$cell_type), FUN=length)
df2$perc_spots = sapply(1:nrow(df2), function(x){df2[x,3]/nrow(spot2cell_roi[which(spot2cell_roi$tile_id == df2[x,1]),])})
colnames(df2)[1:3] = c("tile_id", "cell_type", "n_spots")

wb <- createWorkbook()
addWorksheet(wb, "ROI")
writeData(wb, "ROI", df1, startRow = 1, startCol = 1, rowNames = F, colNames = T)
addWorksheet(wb, "by_tile")
writeData(wb, "by_tile", df2, startRow = 1, startCol = 1, rowNames = F, colNames = T)
saveWorkbook(wb, file = paste0(OUT_DIR, "/spot2cell_table.xlsx"), overwrite = TRUE)


##### Subset input data to select only spot-gene pairs where spots are assigned to cell types
input_data.filter = input_data[which(input_data$spot_id %in% spot2cell_roi$spot_id),]
input_data.filter = merge(input_data.filter, spot2cell_roi[,c("spot_id", "cell_type")], by = "spot_id", all.x = T)
saveRDS(input_data.filter, paste0(OUT_DIR, "/input_data.filter.rds"))

# Stats
print(paste0("Total number of genes expressed in spots assigned to cell types: ", length(unique(input_data.filter$gene_name))))
print(paste0("Median number of genes expressed per spot in spots assigned to cell types: ", median(table(as.factor(input_data.filter$spot_id)))))

### Plot data over the 2D place colored by cell type
df = spot2cell_roi
df$cell_type = factor(df$cell_type)

# Colormaps
colfunc_red <- colorRampPalette(c("#660000", "#ffb3b3")) # red
colfunc_blue <- colorRampPalette(c("#000080", "#8080ff")) # blue
colfunc_gray <- colorRampPalette(c("#333333", "#cccccc")) # gray

custom_colors = c("blue", "green", "red", "black", "orange", "purple")

# Axis labels (850 and 1180 are the estimated sizes in um of a tile)
n_ticks_x = 4
step_x = max(df$col_i-min(df$col_i)) / (n_ticks_x - 1)
breaks_x = c(0, step_x, step_x*(n_ticks_x-2), max(df$col_i-min(df$col_i)))
labels_x = round(breaks_x * 850 * n_tiles_col_ROI / max(df$col_i-min(df$col_i)))

n_ticks_y = 5
step_y = max(df$row_i-min(df$row_i)) / (n_ticks_y - 1)
breaks_y = c(0, step_y, step_y*(n_ticks_y-3), step_y*(n_ticks_y-2), max(df$row_i-min(df$row_i)))
labels_y = round(breaks_y * 1180 * n_tiles_row_ROI / max(df$row_i-min(df$row_i)))

png(paste0(OUT_DIR, "/spot_cell_type.png"), width = 4 * n_tiles_col_ROI, height = 4.5 * n_tiles_row_ROI, res = 400, units = "in")
ggplot(data=df) +
  geom_point(aes(x=col_i-min(col_i), y=max(row_i)-row_i, color=cell_type), size = 0.5) +
  scale_color_manual(values = custom_colors) +
  scale_x_continuous(breaks=breaks_x,
                     labels=paste0(as.character(labels_x), " um")) +
  scale_y_continuous(breaks=breaks_y,
                     labels=paste0(as.character(labels_y), " um")) +
  xlab("") +
  ylab("") +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  theme_bw()
dev.off()

# To add transparency
p <- ggplot(data=df) +
  geom_point(aes(x=col_i-min(col_i), y=max(row_i)-row_i, color=cell_type), size = 0.5) +
  scale_color_manual(values = custom_colors) +
  scale_x_continuous(breaks=breaks_x,
                     labels=paste0(as.character(labels_x), " um")) +
  scale_y_continuous(breaks=breaks_y,
                     labels=paste0(as.character(labels_y), " um")) +
  xlab("") +
  ylab("") +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='transparent'),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.position = "none")
ggsave(p, width = 3.7 * n_tiles_col_ROI, height = 5 * n_tiles_row_ROI, dpi = 400, units = "in",
       filename=paste0(OUT_DIR, "/spot_cell_type_transparent.png"),
       bg="transparent")

# Each cell type separately
for (i in 1:length(unique(as.character(df$cell_type)))){
  p <- ggplot(data=df[which(df$cell_type == unique(as.character(df$cell_type))[i]),]) +
    geom_point(aes(x=col_i-min(col_i), y=max(row_i)-row_i, color=cell_type), size = 0.5) +
    scale_color_manual(values = custom_colors[i]) +
    scale_x_continuous(breaks=breaks_x,
                       labels=paste0(as.character(labels_x), " um")) +
    scale_y_continuous(breaks=breaks_y,
                       labels=paste0(as.character(labels_y), " um")) +
    xlab("") +
    ylab("") +
    guides(colour = guide_legend(override.aes = list(size=2))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill='transparent'),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          plot.background = element_rect(fill='transparent', color=NA),
          legend.position = "none")
  ggsave(p, width = 3.7 * n_tiles_col_ROI, height = 5 * n_tiles_row_ROI, dpi = 400, units = "in",
         filename=paste0(OUT_DIR, "/spot_cell_type_transparent_", unique(as.character(df$cell_type))[i], ".png"),
         bg="transparent")
}


##### Assign cell types to every spot using nearest neighbor approach
input_data.cell_type_NN = input_data[, c(1,2,3,4,8,9)]
input_data.cell_type_NN = input_data.cell_type_NN[!duplicated(input_data.cell_type_NN$spot_id),]
rownames(input_data.cell_type_NN) = input_data.cell_type_NN$spot_id

out_list = pbmcapply::pbmclapply(rownames(input_data.cell_type_NN), 
                                 assign_feature_nearest_neighbor, 
                                 input_data.feature_NN=input_data.cell_type_NN, 
                                 spot2feature_roi=spot2cell_roi, 
                                 feature_label="cell_type",
                                 mc.cores = 32)
input_data.cell_type_NN$cell_type_NN = unlist(out_list)

# Plot data over the 2D place colored by cell type
df = input_data.cell_type_NN
df$cell_type_NN = factor(df$cell_type_NN)

png(paste0(OUT_DIR, "/spot_cell_type_ALL.png"), width = 5.5 * n_tiles_col_ROI, height = 5 * n_tiles_row_ROI, res = 400, units = "in")
ggplot(data=df) +
  geom_point(aes(x=col_i-min(col_i), y=max(row_i)-row_i, color=cell_type_NN), size = 0.1) +
  scale_color_manual(values = custom_colors) +
  scale_x_continuous(breaks=breaks_x,
                     labels=paste0(as.character(labels_x), " um")) +
  scale_y_continuous(breaks=breaks_y,
                     labels=paste0(as.character(labels_y), " um")) +
  xlab("") +
  ylab("") +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  theme_bw()
dev.off()



#################### Spaceflow
create_spaceflow_input(input_data,
                       input_data.filter,
                       df_markers=df_cell_type_markers,
                       feature_label="cell_type",
                       spaceflow_option=1)

create_spaceflow_input(input_data,
                       input_data.filter,
                       df_markers=df_cell_type_markers,
                       feature_label="cell_type",
                       spaceflow_option=2)

# Option 3
df = input_data.filter
df$cell_type = factor(df$cell_type)

xmin_rect = 10000
xmax_rect = 40000
ymin_rect = 30000
ymax_rect = 60000

ggplot(data=df) +
  geom_point(aes(x=col_i-min(col_i), y=max(row_i)-row_i, color=cell_type), size = 0.5) +
  scale_color_manual(values = custom_colors) +
  geom_rect(xmin=xmin_rect, xmax=xmax_rect, ymin=ymin_rect, ymax=ymax_rect, color="black", fill=NA) +
  scale_x_continuous(breaks=breaks_x,
                     labels=paste0(as.character(labels_x), " um")) +
  scale_y_continuous(breaks=breaks_y,
                     labels=paste0(as.character(labels_y), " um")) +
  xlab("") +
  ylab("") +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  theme_bw()

create_spaceflow_input(input_data,
                       input_data.filter,
                       df_markers=df_cell_type_markers,
                       feature_label="cell_type",
                       spaceflow_option=3,
                       xmin_rect,
                       xmax_rect,
                       ymin_rect,
                       ymax_rect)


##### Run SpaceFlow on Python: /mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/result/spaceflow.py


##### Analyze SpaceFlow output

# Read output and create dataframe with spatial domains
df_spatial_domains = read_spaceflow_output(input_data,
                                           input_data.filter,
                                           df_markers=df_cell_type_markers,
                                           feature_label="cell_type",
                                           spaceflow_option=1,
                                           OUT_DIR,
                                           xmin_rect,
                                           xmax_rect,
                                           ymin_rect,
                                           ymax_rect)

# Assign cell types to spatial domains
assign_cell_types_spatial_domains(df_spatial_domains,
                                  spaceflow_option=1,
                                  OUT_DIR,
                                  OR_thres=1,
                                  p_value_thres=0.01)

plot_cluster_feature_dotplot(OUT_DIR,
                             spaceflow_option=1,
                             feature_label="cell_type",
                             p_value_thres=0.01,
                             p_value_cutoff=NA)


# Plot spatial domains
df_spatial_domains = plot_spatial_domains(df_spatial_domains,
                                          OUT_DIR,
                                          out_filename="spot_spatial_domain_map.png",
                                          spaceflow_option=1,
                                          feature_label="cell_type",
                                          only_cluster_assigned_to_feature=T,
                                          spatial_domains_updated_names=c(),
                                          cluster_colors=c(),
                                          xmin_rect,
                                          xmax_rect,
                                          ymin_rect,
                                          ymax_rect)


# Assign spatial domains to all the spots using nearest neighbor
input_data.spatial = assign_spatial_domain_nearest_neighbor(input_data,
                                                            df_spatial_domains,
                                                            spatial_domains_updated_names=c())

input_data.spatial = plot_spatial_domains(df_spatial_domains=input_data.spatial,
                                          OUT_DIR,
                                          out_filename="spot_spatial_domain_map_ALL.png",
                                          spaceflow_option,
                                          feature_label="cell_type",
                                          spatial_domains_updated_names=c(),
                                          cluster_colors=c(),
                                          xmin_rect,
                                          xmax_rect,
                                          ymin_rect,
                                          ymax_rect)



##### Spatial cell-cell interaction enrichment and ligand-receptor analysis
for (i in c(0.1, 0.05, 0.01)){
  giotto_analysis(input_data.filter.rds = paste0(OUT_DIR, "/input_data.filter.rds"),
                  tiles = my_giotto_tiles,
                  pct_thres = i,
                  OUT_DIR)
}

for (i in c(0.1, 0.05, 0.01)){
  giotto_plot(input_dir=paste0(OUT_DIR, "/giotto_", i*100))
}

# Stats
for (i in c(0.1, 0.05, 0.01)){
  temp = data.frame(fread(paste0(OUT_DIR, "/giotto_", i*100, "/expression_matrix.txt")))
  print(paste0("Number of genes - ", i, ": ", nrow(temp)))
  print(paste0("Number of spots - ", i, ": ", ncol(temp)-1))
}


















