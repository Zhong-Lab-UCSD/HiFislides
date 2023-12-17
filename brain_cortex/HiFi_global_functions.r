library(data.table)
library(ggplot2)
library(openxlsx)
library(readxl)
library(reshape2)
library(org.Hs.eg.db)
library(Giotto)
library(circlize)
library(ComplexHeatmap)
library(Matrix)
library(cowplot)
library(ggplotify)
library(igraph)
library(clusterProfiler)

######################## Global functions for HiFi-slides

# Annotation hg38
gencode.v41.annotation.gene = data.frame(fread("/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/hg38_annotation/gencode.v41.annotation.gene.gtf", nThread = 16))
gencode.v41.annotation.gene$gene_id = sapply(gencode.v41.annotation.gene$V9, function(x){gsub('\\"','', strsplit(strsplit(x, "; ")[[1]][1], " ")[[1]][2])})
gencode.v41.annotation.gene$gene_type = sapply(gencode.v41.annotation.gene$V9, function(x){gsub('\\"','', strsplit(strsplit(x, "; ")[[1]][2], " ")[[1]][2])})
gencode.v41.annotation.gene$gene_name = sapply(gencode.v41.annotation.gene$V9, function(x){gsub('\\"','', strsplit(strsplit(x, "; ")[[1]][3], " ")[[1]][2])})
gencode.v41.annotation.gene = gencode.v41.annotation.gene[,-9]

# Annotation mm10
gencode.vM25.annotation.gene = data.frame(fread("/mnt/extraids/SDSC_NFS/rcalandrelli/Genomes/Mus_musculus/gencode/gencode.vM25.annotation.gene.gtf", nThread = 16))
gencode.vM25.annotation.gene$gene_id = sapply(gencode.vM25.annotation.gene$V9, function(x){gsub('\\"','', strsplit(strsplit(x, "; ")[[1]][1], " ")[[1]][2])})
gencode.vM25.annotation.gene$gene_type = sapply(gencode.vM25.annotation.gene$V9, function(x){gsub('\\"','', strsplit(strsplit(x, "; ")[[1]][2], " ")[[1]][2])})
gencode.vM25.annotation.gene$gene_name = sapply(gencode.vM25.annotation.gene$V9, function(x){gsub('\\"','', strsplit(strsplit(x, "; ")[[1]][3], " ")[[1]][2])})
gencode.vM25.annotation.gene = gencode.vM25.annotation.gene[,-9]

'%!in%' <- function(x,y)!('%in%'(x,y))

#### Create the matrix with tiles ####
# Offsets were chosen in order to combine spots from different tiles into one figure.
# There is no physical distance between each two neighboring tiles.
create_tile_matrix <- function(flowcell_type,
                               first_tile_position,
                               n_tiles_row,
                               n_tiles_col,
                               surface,
                               n_bin=NA){
  if (flowcell_type == "NextSeq"){
    
    if (first_tile_position == "top-left"){
      x = 1101:(1101 + n_tiles_row-1)
      tile_matrix = x
      for (i in 1:(n_tiles_col-1)){
        tile_matrix = cbind(tile_matrix, x+i*100)
      }
    }
    else if (first_tile_position == "bottom-right"){
      x = (1000 + n_tiles_col*100 + n_tiles_row):(1000 + n_tiles_col*100 + 1)
      tile_matrix = x
      for (i in 1:(n_tiles_col-1)){
        tile_matrix = cbind(tile_matrix, x-i*100)
      }
    }
    
  } else if (flowcell_type == "MiniSeq"){
    if (surface == "1"){
      x = 11101:11104
      tile_matrix = cbind(x, x+1000, x+2000)
    } else {
      x = 21101:21104
      tile_matrix = cbind(x, x+1000, x+2000)
    }
  }
  
  # Bin tiles
  if (!is.na(n_bin)){
    template_matrix = matrix(0,nrow=n_bin,ncol=n_bin)
    template_matrix_expanded = transform(expand.grid(i = seq(nrow(template_matrix)), j = seq(ncol(template_matrix))), v = c(template_matrix))
    template_matrix_expanded$v = paste0(template_matrix_expanded$i, "_", template_matrix_expanded$j)
    
    temp_list = list()
    for (i in as.numeric(tile_matrix)){
      temp_list[[as.character(i)]] = matrix(paste0(i, "_", template_matrix_expanded$v), n_bin, n_bin)
    }
    
    if (flowcell_type == "NextSeq"){
      if (first_tile_position == "top-left"){
        tile_matrix = do.call("rbind", temp_list[as.character(x)])
        for (i in 1:(n_tiles_col-1)){
          tile_matrix = cbind(tile_matrix, do.call("rbind", temp_list[as.character(x+i*100)]))
        }
      }
      else if (first_tile_position == "bottom-right"){
        tile_matrix = do.call("rbind", temp_list[as.character(x)])
        for (i in 1:(n_tiles_col-1)){
          tile_matrix = cbind(tile_matrix, do.call("rbind", temp_list[as.character(x-i*100)]))
        }
      }
    }
    else if (flowcell_type == "MiniSeq"){
        tile_matrix = do.call("rbind", temp_list[as.character(x)])
        for (i in c(1000,2000)){
          tile_matrix = cbind(tile_matrix, do.call("rbind", temp_list[as.character(x+i)]))
        }
    }
    
  }
  
  return(tile_matrix)
}


#### Plot spot and gene heatmaps of the flowcell #### 
plot_flowcell_heatmaps <- function(INTEGRATE_DIR=paste0(DATA_DIR, "/L2R1_L2R2_integrate/"),
                                   SAMPLE_NAME,
                                   custom_tiles=NA, # "top" or "bottom" to only plot those tiles
                                   normalization_file=NA, # to normalize tiles (tab separated with tile ID in column 1 and normalization values in column 2)
                                   scaling_factor=1,
                                   n_tiles_row,
                                   n_tiles_col,
                                   n_min_read=NA,
                                   n_max_read=NA,
                                   n_min_spot=NA,
                                   n_max_spot=NA,
                                   n_min_gene=NA,
                                   n_max_gene=NA)
  
  {
  
  ### Plot normalizing factor heatmap if present
  if(!is.na(normalization_file)){
    temp = read.table(normalization_file)
    rownames(temp) = temp$V1
    
    temp_matrix = matrix(nrow=n_tiles_row, ncol=n_tiles_col)
    for (i in 1:ncol(tile_matrix)){
      temp_matrix[,i] = temp[as.character(tile_matrix[,i]),2]
    }
    
    col_fun = colorRamp2(breaks = seq(min(temp_matrix[!is.na(temp_matrix)]), 
                                      max(temp_matrix[!is.na(temp_matrix)]),
                                      length=3), 
                         colors=c("blue","white","red"))
    
    
    p1<-Heatmap(temp_matrix, cluster_rows = F, cluster_columns = F, na_col = "black", col = col_fun, column_title = strsplit(tail(strsplit(normalization_file, "/")[[1]], 1), "\\.")[[1]][1],
                heatmap_legend_param = list(title = ""))
    png(paste0(paste(strsplit(OUT_DIR, "/")[[1]][1:(length(strsplit(OUT_DIR, "/")[[1]])-1)], collapse = "/"), "/",
               strsplit(tail(strsplit(normalization_file, "/")[[1]], 1), "\\.")[[1]][1], "_heatmap.png"), width = 2.5, height = 5, res = 200, units = "in")
    print(p1)
    dev.off()
    
  }

  ### Read pair number
  temp = read.table(paste0(INTEGRATE_DIR, "/tile_read_number_table.txt"))
  rownames(temp) = temp$V1
 
  if (!is.na(custom_tiles)){
    temp = temp[which(temp$V1 %in% custom_tiles),]
  }
  
  if(!is.na(normalization_file)){
    norm_data = read.table(normalization_file)
    rownames(norm_data) = norm_data$V1
    temp$V2 = temp$V2 / merge(temp, norm_data, by = "V1")$V2.y * scaling_factor
  }
  
  temp_matrix = matrix(nrow=n_tiles_row, ncol=n_tiles_col)
  for (i in 1:ncol(tile_matrix)){
    temp_matrix[,i] = temp[as.character(tile_matrix[,i]),2]
  }
  # temp_matrix[which(is.na(temp_matrix), arr.ind = T)] = 0
  
  col_fun = colorRamp2(breaks = seq(ifelse(is.na(n_min_read), min(temp_matrix[!is.na(temp_matrix)]), n_min_read), 
                                    ifelse(is.na(n_max_read), max(temp_matrix[!is.na(temp_matrix)]), n_max_read),
                                    length=3), 
                       colors=c("blue","white","red"))
  
  
  p1<-Heatmap(temp_matrix, cluster_rows = F, cluster_columns = F, na_col = "black", col = col_fun, column_title = paste0("Number of read pairs"),
              heatmap_legend_param = list(title = ""))
  png(paste0(paste(strsplit(OUT_DIR, "/")[[1]][1:(length(strsplit(OUT_DIR, "/")[[1]])-1)], collapse = "/"), "/read_number_heatmap.png"), width = 2.5, height = 5, res = 200, units = "in")
  print(p1)
  dev.off()
  
  
  ### Spot number
  temp = read.table(paste0(INTEGRATE_DIR, "/tile_spot_number_table.txt"))
  rownames(temp) = temp$V1
  
  if (!is.na(custom_tiles)){
    temp = temp[which(temp$V1 %in% custom_tiles),]
  }
  
  if(!is.na(normalization_file)){
    norm_data = read.table(normalization_file)
    rownames(norm_data) = norm_data$V1
    temp$V2 = temp$V2 / merge(temp, norm_data, by = "V1")$V2.y * scaling_factor
  }
  
  temp_matrix = matrix(nrow=n_tiles_row, ncol=n_tiles_col)
  for (i in 1:ncol(tile_matrix)){
    temp_matrix[,i] = temp[as.character(tile_matrix[,i]),2]
  }
  # temp_matrix[which(is.na(temp_matrix), arr.ind = T)] = 0
  
  col_fun = colorRamp2(breaks = seq(ifelse(is.na(n_min_spot), min(temp_matrix[!is.na(temp_matrix)]), n_min_spot), 
                                    ifelse(is.na(n_max_spot), max(temp_matrix[!is.na(temp_matrix)]), n_max_spot),
                                    length=3), 
                       colors=c("blue","white","red"))
  
  
  p1<-Heatmap(temp_matrix, cluster_rows = F, cluster_columns = F, na_col = "black", col = col_fun, column_title = paste0("Number of spots"),
              heatmap_legend_param = list(title = ""))
  png(paste0(paste(strsplit(OUT_DIR, "/")[[1]][1:(length(strsplit(OUT_DIR, "/")[[1]])-1)], collapse = "/"), "/spot_number_heatmap.png"), width = 2.5, height = 5, res = 200, units = "in")
  print(p1)
  dev.off()
  
  
  ### Expressed gene number
  temp = read.table(paste0(INTEGRATE_DIR, "/tile_gene_number_table.txt"))
  rownames(temp) = temp$V1
  
  if (!is.na(custom_tiles)){
    temp = temp[which(temp$V1 %in% custom_tiles),]
  }
  
  if(!is.na(normalization_file)){
    norm_data = read.table(normalization_file)
    rownames(norm_data) = norm_data$V1
    temp$V2 = temp$V2 / merge(temp, norm_data, by = "V1")$V2.y * scaling_factor
  }
  
  temp_matrix = matrix(nrow=n_tiles_row, ncol=n_tiles_col)
  for (i in 1:ncol(tile_matrix)){
    temp_matrix[,i] = temp[as.character(tile_matrix[,i]),2]
  }
  # temp_matrix[which(is.na(temp_matrix), arr.ind = T)] = 0
  
  temp_matrix[which(temp_matrix == min(temp_matrix), arr.ind = T)] = NA
  col_fun = colorRamp2(breaks = seq(ifelse(is.na(n_min_gene), min(temp_matrix[!is.na(temp_matrix)]), n_min_gene), 
                                    ifelse(is.na(n_max_gene), max(temp_matrix[!is.na(temp_matrix)]), n_max_gene),
                                    length=3), 
                       colors=c("blue","white","red"))
  
  p2<-Heatmap(temp_matrix, cluster_rows = F, cluster_columns = F, na_col = "black", col = col_fun, column_title = paste0("Number of expressed genes"),
              heatmap_legend_param = list(title = ""))
  png(paste0(paste(strsplit(OUT_DIR, "/")[[1]][1:(length(strsplit(OUT_DIR, "/")[[1]])-1)], collapse = "/"), "/gene_number_heatmap.png"), width = 2.5, height = 5, res = 200, units = "in")
  print(p2)
  dev.off()
  
  return()
}


#### Reading input data #### 
read_HiFi_input_data <- function(INPUT_DIR,
                                 SAMPLE_NAME,
                                 ROI_label,
                                 tile_matrix){
  
  input_data = data.frame(fread(paste0(INPUT_DIR, "/", SAMPLE_NAME, ".L2R2_genome_spatial.final.", ROI_label, ".txt"), nThread=24))
  colnames(input_data) = c("spot_id", "tile_id", "col", "row", "gene_name", "gene_biotype", "expression")
  
  # input_data_list = list()
  # for (i in names(ROI_tiles_list)){
  #   input_data = data.frame(fread(paste0(DATA_DIR, "/L2R1_L2R2_integrate/", SAMPLE_NAME, "/HiFi_L2R2_genome_spatial.final.", i, ".txt")))
  #   colnames(input_data) = c("spot_id", "tile_id", "col", "row", "gene_name", "gene_biotype", "expression")
  #   input_data_list[[i]] = input_data
  # }
  # 
  # if (length(input_data_list) > 1){
  #   input_data = do.call("rbind", input_data_list)
  # } else {
  #   input_data = input_data_list[[1]]
  #   rm(input_data_list)
  # }
  
  if (flowcell_type == "NextSeq"){
    row_tile_offset = 74500
    col_tile_offset = 56500
  } else if (flowcell_type == "MiniSeq"){
    row_tile_offset = 30000
    col_tile_offset = 22000
  }
  
  # 2D row value
  calculate_row_i <- function(i){
    tile_row = which(tile_matrix == input_data[i,"tile_id"], arr.ind=T)[1]
    d_row = row_tile_offset * (tile_row - 1)
    row_i = as.numeric(input_data[i,"row"]) + d_row
    return(row_i)
  }
  input_data$row_i = unlist(pbmcapply::pbmclapply(1:nrow(input_data), calculate_row_i, mc.cores = 32))
  
  # 2D col value
  calculate_col_i <- function(i){
    tile_col = which(tile_matrix == input_data[i,"tile_id"], arr.ind=T)[2]
    d_col = col_tile_offset * (tile_col - 1)
    col_i = as.numeric(input_data[i,"col"]) + d_col
    return(col_i)
  }
  input_data$col_i = unlist(pbmcapply::pbmclapply(1:nrow(input_data), calculate_col_i, mc.cores = 32))
  
  return(input_data)
}

#### Add binned tiles to input data #### 
assign_binned_tile <- function(input_data,
                               tile_matrix,
                               n_bin){
  
  n_tiles = length(unique(input_data$tile_id))
  
  assign_bin_to_spot <- function(my_tile){
    temp = input_data[which(input_data$tile_id == my_tile),]
    temp$tile_id_bin = ""
    
    split_ranges_row = levels(ggplot2::cut_number(min(temp$row):max(temp$row), n_bin))
    split_ranges_col = levels(ggplot2::cut_number(min(temp$col):max(temp$col), n_bin))
    
    for (i in 1:length(split_ranges_row)){
      min_row = as.numeric(gsub("\\[|\\(", "", strsplit(split_ranges_row[i], ",")[[1]][1]))
      max_row = as.numeric(gsub("\\]|\\)", "", strsplit(split_ranges_row[i], ",")[[1]][2]))
      for (j in 1:length(split_ranges_col)){
        min_col = as.numeric(gsub("\\[|\\(", "", strsplit(split_ranges_col[j], ",")[[1]][1]))
        max_col = as.numeric(gsub("\\]|\\)", "", strsplit(split_ranges_col[j], ",")[[1]][2]))
        temp[which(temp$row >= min_row & temp$row <= max_row & temp$col >= min_col & temp$col <= max_col), "tile_id_bin"] = paste0(my_tile, "_", i, "_", j)
      }
    }
    
    return(temp)
    
  }
  
  out_list = pbmcapply::pbmclapply(unique(input_data$tile_id), assign_bin_to_spot, mc.cores = n_tiles)
  out_df = do.call("rbind", out_list)
  return(out_df)
  
}



#### ROI statistics #### 
calculate_ROI_stats <- function(INPUT_DIR,
                                SAMPLE_NAME,
                                OUT_DIR,
                                input_data){
  
  df1 = data.frame(spot_gene_pairs = nrow(input_data),
                   n_spots = length(unique(input_data$spot_id)),
                   n_genes = length(unique(input_data$gene_name)),
                   median_genes_per_spot = median(table(as.factor(input_data$spot_id))))
  
  
  temp = read.table(paste0(INPUT_DIR, "/tile_spot_number_table.txt"))
  temp = temp[which(temp$V1 %in% unique(input_data$tile_id)),]
  rownames(temp) = temp$V1
  df2 = data.frame(avg_spot_tile = mean(temp$V2),
                   avg_spot_10um2 = mean(temp$V3),
                   tile_name_max_spot = temp[which(temp$V2 == max(temp$V2)),"V1"],
                   tile_max_spot = max(temp$V2),
                   tile_max_spot_10um2 = max(temp$V3))
  
  temp = read.table(paste0(INPUT_DIR, "/tile_gene_number_table.txt"))
  temp = temp[which(temp$V1 %in% unique(input_data$tile_id)),]
  rownames(temp) = temp$V1
  df3 = data.frame(avg_gene_tile = mean(temp$V2),
                   avg_gene_10um2 = mean(temp$V3),
                   tile_name_max_gene = temp[which(temp$V2 == max(temp$V2)),"V1"],
                   tile_max_gene = max(temp$V2),
                   tile_max_gene_10um2 = max(temp$V3))
  
  df = cbind(df1, df2, df3)
  df = data.frame(t(df))
  df$V1 = rownames(df)
  df$V2 = as.numeric(df[,1])
  df = df[,-1]
  df$V2 = round(df$V2, 2)
  write.xlsx(df, paste0(OUT_DIR, "/ROI_stats.xlsx"), rowNames = F, colNames = F)
  
  return()
}


#### Plot spots with gene expression or tile features related to genes ####
plot_gene_feature_spots <- function(input_data,
                                    my_genes,
                                    feature,
                                    tile_matrix,
                                    aggregate_field="tile_id",
                                    n_bin_tile=NA,
                                    normalization_file,
                                    scaling_factor=1,
                                    expr_threshold=NA, # minimum gene expression level to filter out spots
                                    quantile_plotting_thres=NA,
                                    out_label){
  
  input_data_gene = input_data[which(input_data$gene_name %in% my_genes),]
  
  if (!is.na(expr_threshold)){
    input_data_gene = input_data_gene[which(input_data_gene$expression > expr_threshold),]
  }
  
  input_data_gene = aggregate(input_data_gene$expression, by=list(tile_id=input_data_gene[,aggregate_field], spot_id=input_data_gene$spot_id, row_i=input_data_gene$row_i, col_i=input_data_gene$col_i), FUN=sum)
  colnames(input_data_gene)[5] = "expression"
  
  # Normalize data
  if (!is.na(normalization_file)){
    temp = read.table(normalization_file)
    rownames(temp) = temp$V1
    
    tile_matrix_norm = create_tile_matrix(flowcell_type,
                                          first_tile_position,
                                          n_tiles_row,
                                          n_tiles_col,
                                          surface)
    
    norm_matrix = matrix(nrow=nrow(tile_matrix_norm), ncol=ncol(tile_matrix_norm))
    for (i in 1:ncol(tile_matrix_norm)){
      norm_matrix[,i] = temp[as.character(tile_matrix_norm[,i]),2]
    }
    
    if (!is.na(n_bin_tile)){
      n1 = norm_matrix[rep(seq_len(nrow(norm_matrix)), each = n_bin_tile),]
      norm_matrix = n1[,rep(seq_len(ncol(n1)), each = n_bin_tile)]
    }
    
  }
  
  ### Plot spots with gene expression
  if (feature == "spots"){
    df = input_data_gene
    
    n_ticks_x = 4
    step_x = max(df$col_i-min(df$col_i)) / (n_ticks_x - 1)
    breaks_x = c(0, step_x, step_x*(n_ticks_x-2), max(df$col_i-min(df$col_i)))
    labels_x = round(breaks_x * width_tile * n_tiles_col_ROI / max(df$col_i-min(df$col_i)))
    
    n_ticks_y = 5
    step_y = max(df$row_i-min(df$row_i)) / (n_ticks_y - 1)
    breaks_y = c(0, step_y, step_y*(n_ticks_y-3), step_y*(n_ticks_y-2), max(df$row_i-min(df$row_i)))
    labels_y = round(breaks_y * height_tile * n_tiles_row_ROI / max(df$row_i-min(df$row_i)))
    
    if (!is.na(expr_threshold)){
      df = df[which(df$expression > expr_threshold),]
    }
    
    p1 <- ggplot(data=df) +
      geom_point(aes(x=col_i-min(col_i), y=max(row_i)-row_i, color=log(expression))) +
      scale_colour_continuous(low = "blue", high = "red") +
      scale_x_continuous(breaks=breaks_x,
                         labels=paste0(as.character(labels_x), " um")) +
      scale_y_continuous(breaks=breaks_y,
                         labels=paste0(as.character(labels_y), " um")) +
      xlab("") +
      ylab("") +
      theme_bw() +
      ggtitle(out_label)
    
    return(list(p1, input_data_gene))
    
  }
  
  ### Plot number of spots per tile expressing my_genes
  if (feature == "number_spots_per_tile"){
    temp1 = input_data_gene[!duplicated(input_data_gene$spot_id),]
    temp1 = aggregate(temp1$expression, by=list(temp1$tile_id), FUN=length)
    tile_matrix_index = transform(expand.grid(i = seq(nrow(tile_matrix)), j = seq(ncol(tile_matrix))), v = c(tile_matrix))
    temp_matrix = matrix(nrow=nrow(tile_matrix), ncol=ncol(tile_matrix))
    for (i in temp1$Group.1){
      temp_matrix[tile_matrix_index[which(tile_matrix_index$v == i), "i"],tile_matrix_index[which(tile_matrix_index$v == i), "j"]] = temp1[which(temp1$Group.1 == i), "x"]
    }
    
    # if (!is.na(normalization_file)){
    #   if (aggregate_field == "tile_id"){
    #     temp_matrix = temp_matrix / norm_matrix
    #   } else {
    #     # In this case each binned tile is normalized by the tile value
    #     tile_matrix_norm_index = transform(expand.grid(i = seq(nrow(tile_matrix_norm)), j = seq(ncol(tile_matrix_norm))), v = c(tile_matrix_norm))
    #     temp_matrix_index = which(!is.na(temp_matrix), arr.ind = T)
    #     for (z in 1:nrow(temp_matrix_index)){
    #       # Find original tile
    #       temp_tile = strsplit(tile_matrix_index[which(tile_matrix_index$i == as.numeric(temp_matrix_index[z,"row"]) & 
    #                                          tile_matrix_index$j == as.numeric(temp_matrix_index[z,"col"])), "v"], "_")[[1]][1]
    #       temp_norm_value = norm_matrix[tile_matrix_norm_index[which(tile_matrix_norm_index$v == temp_tile),"i"], tile_matrix_norm_index[which(tile_matrix_norm_index$v == temp_tile),"j"]]
    #       temp_matrix[temp_matrix_index[z,"row"], temp_matrix_index[z,"col"]] = temp_matrix[temp_matrix_index[z,"row"], temp_matrix_index[z,"col"]] / temp_norm_value
    #     }
    #   }
    # }
    
    out_matrix = temp_matrix # to output the number of spots expressing the genes
    out_matrix[which(out_matrix==0), arr.ind = T] = NA
    
    if (!is.na(normalization_file)){
      temp_matrix = temp_matrix / norm_matrix * scaling_factor
    }
    
    temp_matrix[which(temp_matrix==0), arr.ind = T] = NA
    
    if (is.na(quantile_plotting_thres)){
      col_fun = colorRamp2(breaks = seq(min(temp_matrix[!is.na(temp_matrix)]),
                                        max(temp_matrix[!is.na(temp_matrix)]),
                                        length=3),
                           colors=c("blue","white","red"))
    } else {
      col_fun = colorRamp2(breaks = seq(quantile(temp_matrix[!is.na(temp_matrix)], quantile_plotting_thres),
                                        quantile(temp_matrix[!is.na(temp_matrix)], 1-quantile_plotting_thres),
                                        length=3),
                           colors=c("blue","white","red"))
    }

    p1 <- Heatmap(temp_matrix, cluster_rows = F, cluster_columns = F, na_col = "black", col = col_fun, column_title = out_label,
                  heatmap_legend_param = list(title = ""))
    return(list(p1, out_matrix))
  }
  
  
  ### Plot total gene expression of spots per tile
  if (feature == "gene_expression_per_tile"){
    temp1 = aggregate(input_data_gene$expression, by=list(input_data_gene$tile_id), FUN=sum)
    tile_matrix_index = transform(expand.grid(i = seq(nrow(tile_matrix)), j = seq(ncol(tile_matrix))), v = c(tile_matrix))
    temp_matrix = matrix(nrow=n_tiles_row, ncol=n_tiles_col)
    for (i in temp1$Group.1){
      temp_matrix[tile_matrix_index[which(tile_matrix_index$v == i), "i"],tile_matrix_index[which(tile_matrix_index$v == i), "j"]] = temp1[which(temp1$Group.1 == i), "x"]
    }
    
    if (!is.na(normalization_file)){
      if (aggregate_field == "tile_id"){
        temp_matrix = temp_matrix / norm_matrix
      } else {
        # In this case each binned tile is normalized by the tile value
        tile_matrix_norm_index = transform(expand.grid(i = seq(nrow(tile_matrix_norm)), j = seq(ncol(tile_matrix_norm))), v = c(tile_matrix_norm))
        temp_matrix_index = which(!is.na(temp_matrix), arr.ind = T)
        for (z in 1:nrow(temp_matrix_index)){
          # Find original tile
          temp_tile = strsplit(tile_matrix_index[which(tile_matrix_index$i == as.numeric(temp_matrix_index[z,"row"]) & 
                                                         tile_matrix_index$j == as.numeric(temp_matrix_index[z,"col"])), "v"], "_")[[1]][1]
          temp_norm_value = norm_matrix[tile_matrix_norm_index[which(tile_matrix_norm_index$v == temp_tile),"i"], tile_matrix_norm_index[which(tile_matrix_norm_index$v == temp_tile),"j"]]
          temp_matrix[temp_matrix_index[z,"row"], temp_matrix_index[z,"col"]] = temp_matrix[temp_matrix_index[z,"row"], temp_matrix_index[z,"col"]] / temp_norm_value
        }
      }
    }
    
    col_fun = colorRamp2(breaks = seq(min(temp_matrix[!is.na(temp_matrix)]),
                                      max(temp_matrix[!is.na(temp_matrix)]),
                                      length=3),
                         colors=c("blue","white","red"))
    
    p1<-Heatmap(temp_matrix, cluster_rows = F, cluster_columns = F, na_col = "black", col = col_fun, column_title = out_label,
                heatmap_legend_param = list(title = ""))
    return(list(p1, out_matrix))
  }
  
}



#### Gradient of a feature across tiles along a direction ####
gradient_gene_feature_tiles <- function(my_gene,
                                        input_data,
                                        ROI_tiles,
                                        direction, # vert, horiz, diag1 (from top-right to bottom-left), diag2 (from top-left to bottom-right)
                                        feature, # spot_number, gene_expr
                                        input_tile_matrix,
                                        aggregate_field="tile_id",
                                        n_bin_tile=NA,
                                        normalization_file=NA,
                                        scaling_factor=1,
                                        expr_threshold=NA){
  
  
  # Subset data to select spots expressing the gene of interest
  input_data_temp = input_data[which(input_data$gene_name == my_gene),]
  
  if (!is.na(expr_threshold)){
    input_data_temp = input_data_temp[which(input_data_temp$expression > expr_threshold),]
  }
  
  if (feature == "gene_expr"){
    # Expression level by tile
    temp1 = aggregate(input_data_temp$expression, by=list(input_data_temp[,aggregate_field]), FUN=sum)
  }
  if (feature == "spot_number"){
    # Number of spots per tile
    temp1 = aggregate(input_data_temp$expression, by=list(input_data_temp[,aggregate_field]), FUN=length)
  }
  
  tile_matrix_index = transform(expand.grid(i = seq(nrow(input_tile_matrix)), j = seq(ncol(input_tile_matrix))), v = c(input_tile_matrix))
  temp_matrix = matrix(0, nrow=nrow(input_tile_matrix), ncol=ncol(input_tile_matrix))
  rownames(temp_matrix) = 1:nrow(temp_matrix)
  colnames(temp_matrix) = 1:ncol(temp_matrix)
  for (i in temp1$Group.1){
    temp_matrix[tile_matrix_index[which(tile_matrix_index$v == i), "i"],tile_matrix_index[which(tile_matrix_index$v == i), "j"]] = temp1[which(temp1$Group.1 == i), "x"]
  }
  
  # Selection of the ROI
  if (aggregate_field == "tile_id"){
    temp = tile_matrix_index$v %in% ROI_tiles
  } else {
    temp = sapply(tile_matrix_index$v, function(x){strsplit(x,"_")[[1]][1] %in% ROI_tiles})
  }
  
  tile_matrix_index_ROI = tile_matrix_index[which(temp==1),] # indexes of entries under ROI
  temp_matrix_ROI = temp_matrix[min(tile_matrix_index[temp,"i"]):max(tile_matrix_index[temp,"i"]),
                                min(tile_matrix_index[temp,"j"]):max(tile_matrix_index[temp,"j"])]
  # temp_matrix_ROI = temp_matrix[min(which(!is.na(temp_matrix), arr.ind = T)[,1]):max(which(!is.na(temp_matrix), arr.ind = T)[,1]),
  #                               min(which(!is.na(temp_matrix), arr.ind = T)[,2]):max(which(!is.na(temp_matrix), arr.ind = T)[,2])]
  
  # Put NA where there are tile entries not under ROI (these NAs are used to remove entries from the normalization matrix, if used)
  index.temp_matrix_ROI = transform(expand.grid(i = rownames(temp_matrix_ROI), j = colnames(temp_matrix_ROI)))
  index.temp_matrix_ROI$ROI = paste0(index.temp_matrix_ROI$i,index.temp_matrix_ROI$j) %in% paste0(tile_matrix_index_ROI$i,tile_matrix_index_ROI$j)
  index.temp_matrix_noROI = index.temp_matrix_ROI[which(index.temp_matrix_ROI$ROI == FALSE),]
  temp_matrix_ROI[as.matrix(index.temp_matrix_noROI[,1:2])] = NA
  
  # Normalize data
  if (!is.na(normalization_file)){
    temp = read.table(normalization_file)
    rownames(temp) = temp$V1
    
    tile_matrix_norm = create_tile_matrix(flowcell_type,
                                          first_tile_position,
                                          n_tiles_row,
                                          n_tiles_col,
                                          surface)
    
    norm_matrix = matrix(nrow=nrow(tile_matrix_norm), ncol=ncol(tile_matrix_norm))
    for (i in 1:ncol(tile_matrix_norm)){
      norm_matrix[,i] = temp[as.character(tile_matrix_norm[,i]),2]
    }
    
    if (!is.na(n_bin_tile)){
      n1 = norm_matrix[rep(seq_len(nrow(norm_matrix)), each = n_bin_tile),]
      norm_matrix = n1[,rep(seq_len(ncol(n1)), each = n_bin_tile)]
    }
    
    if (aggregate_field == "tile_id"){
      temp = tile_matrix_index$v %in% ROI_tiles
    } else {
      temp = sapply(tile_matrix_index$v, function(x){strsplit(x,"_")[[1]][1] %in% ROI_tiles})
    }
    norm_matrix = norm_matrix[min(tile_matrix_index[temp,"i"]):max(tile_matrix_index[temp,"i"]),
                              min(tile_matrix_index[temp,"j"]):max(tile_matrix_index[temp,"j"])]
    # norm_matrix = norm_matrix[min(which(!is.na(temp_matrix), arr.ind = T)[,1]):max(which(!is.na(temp_matrix), arr.ind = T)[,1]),
    #                           min(which(!is.na(temp_matrix), arr.ind = T)[,2]):max(which(!is.na(temp_matrix), arr.ind = T)[,2])]
    norm_matrix[which(is.na(temp_matrix_ROI), arr.ind = T)] = NA
  }
  
  
  if (direction == "vert"){
    
    temp_mat = temp_matrix_ROI
    temp_mat[which(is.na(temp_mat), arr.ind = T)] = 0
    
    if (!is.na(normalization_file)){
      temp_norm_matrix = norm_matrix
      temp_norm_matrix[which(is.na(temp_norm_matrix), arr.ind = T)] = 0
      out_vect = rowSums(temp_mat) / rowSums(temp_norm_matrix)
    } else {
      out_vect = rowSums(temp_mat)
    }
    
  }
  
  else if (direction == "horiz"){
    
    temp_mat = temp_matrix_ROI
    temp_mat[which(is.na(temp_mat), arr.ind = T)] = 0
    
    if (!is.na(normalization_file)){
      temp_norm_matrix = norm_matrix
      temp_norm_matrix[which(is.na(temp_norm_matrix), arr.ind = T)] = 0
      out_vect = colSums(temp_mat) / colSums(temp_norm_matrix)
    } else {
      out_vect = colSums(temp_mat)
    }
    
  }
  
  else if (direction == "diag1"){
    
    out_vect = c()
    for (k in (ncol(temp_matrix_ROI)-1) : -(nrow(temp_matrix_ROI)-1)){
      
      # This matrix is used to check what are the values in the segment under analysis, if they are all NA (region outside the ROI) then skip analysis
      a = triu(matrix(1, nrow(temp_matrix_ROI), ncol(temp_matrix_ROI)), k=k) - triu(matrix(1, nrow(temp_matrix_ROI), ncol(temp_matrix_ROI)), k=k+1)
      segment_values = temp_matrix_ROI[which(a == 1, arr.ind = T)]
      
      if (sum(is.na(segment_values)) != length(segment_values)){
        
        temp_mat = triu(temp_matrix_ROI, k=k) - triu(temp_matrix_ROI, k=k+1)
        # temp_mat[which(temp_mat == 0, arr.ind = T)] = NA
        
        if (!is.na(normalization_file)){
          temp_norm = triu(norm_matrix, k=k) - triu(norm_matrix, k=k+1)
          out_vect = c(out_vect, sum(temp_mat[which(!is.na(temp_mat))])/sum(temp_norm[which(!is.na(temp_mat))]))
          # out_vect = c(out_vect, sum(temp_mat)/sum(temp_norm))
        } else {
          out_vect = c(out_vect, sum(temp_mat[which(!is.na(temp_mat))]))
          # out_vect = c(out_vect, sum(temp_mat))
        }
      }
    }
    
  }
  
  else if (direction == "diag2"){
    
    temp_matrix_ROI = temp_matrix_ROI[,ncol(temp_matrix_ROI):1] # to flip the matrix horizontally to go from top-left to bottom-right
    norm_matrix = norm_matrix[,ncol(norm_matrix):1] # to flip the matrix horizontally to go from top-left to bottom-right
    
    out_vect = c()
    for (k in (ncol(temp_matrix_ROI)-1) : -(nrow(temp_matrix_ROI)-1)){
      
      # This matrix is used to check what are the values in the segment under analysis, if they are all NA (region outside the ROI) then skip analysis
      a = triu(matrix(1, nrow(temp_matrix_ROI), ncol(temp_matrix_ROI)), k=k) - triu(matrix(1, nrow(temp_matrix_ROI), ncol(temp_matrix_ROI)), k=k+1)
      segment_values = temp_matrix_ROI[which(a == 1, arr.ind = T)]
      
      if (sum(is.na(segment_values)) != length(segment_values)){
        
        temp_mat = triu(temp_matrix_ROI, k=k) - triu(temp_matrix_ROI, k=k+1)
        # temp_mat[which(temp_mat == 0, arr.ind = T)] = NA
        
        if (!is.na(normalization_file)){
          temp_norm = triu(norm_matrix, k=k) - triu(norm_matrix, k=k+1)
          out_vect = c(out_vect, sum(temp_mat[which(!is.na(temp_mat))])/sum(temp_norm[which(!is.na(temp_mat))]))
          # out_vect = c(out_vect, sum(temp_mat)/sum(temp_norm))
        } else {
          out_vect = c(out_vect, sum(temp_mat[which(!is.na(temp_mat))]))
          # out_vect = c(out_vect, sum(temp_mat))
        }
      }
    }
    
  }
  
  # Linear regression
  df = data.frame(x = 0:(length(out_vect)-1),
                  y = out_vect)
  
  out.lm = lm(y ~ x, data=df)
  summary.lm = summary(out.lm)
  summary.lm.coeff = summary.lm$coefficients
  
  out_list = list(my_gene, out_vect, summary.lm.coeff)
  
  return(out_list)
  
  # ggplot(df,aes(x, y)) +
  #   geom_point() +
  #   geom_smooth(method='lm', formula= y~x)
  
}



#### Assign features to spots #### 
# If a spot expresses one or more markers of a cell type, it is assigned to that cell type. 
# If a spot expresses marker genes associated with several cell types, it is assigned to the cell type corresponding to the marker with the highest expression level. 
# If a spot expresses marker genes associated with several cell types and two or more marker genes have the same highest expression level, then it is not assigned.
assign_feature_to_spots <- function(input_data,
                                    df_markers,
                                    feature_label, # cell_type or layer
                                    with_mixed_features=F){
  
  spot_feature_list = list() # list with all the spots associated to cell types (spot expression level: sum of the expression level of marker genes expressed in that spot)
  
  for (i in names(table(as.factor(df_markers[,feature_label])))){
    temp = input_data[which(input_data$gene_name %in% df_markers[which(df_markers[,feature_label] == i), "gene"]),]
    if (nrow(temp) > 0){
      temp[,feature_label] = i
      temp1 = aggregate(temp$expression, 
                        list(temp$spot_id, temp$tile_id, temp$col, temp$row, temp$col_i, temp$row_i, temp[,feature_label]),
                        sum)
      colnames(temp1) = c("spot_id", "tile_id", "col", "row", "col_i", "row_i", feature_label, "expression")
      spot_feature_list[[i]] = temp1
    }
  }
  
  df_spot_feature = do.call("rbind", spot_feature_list)
  
  ### Treat eventual duplicated spots
  '%!in%' <- function(x,y)!('%in%'(x,y))
  
  spot_id_duplicated = names(table(as.factor(df_spot_feature$spot_id))[(table(as.factor(df_spot_feature$spot_id)) > 1)])
  
  if (length(spot_id_duplicated) == 0){
    spot2feature_roi = df_spot_feature
  
  } else {
    
    find_rownames_to_keep <- function(x){
      dup_spots = df_spot_feature[which(df_spot_feature$spot_id == x),]
      
      # With mixed cells
      if (with_mixed_features == T){
        out_df = data.frame(V1="Mix",
                            V2=mean(dup_spots$expression),
                            V3=rownames(dup_spots[1,]))
        
        return(out_df)
        
      } else {
        dup_spots_agg = aggregate(dup_spots$expression, 
                                  list(dup_spots[,feature_label]),
                                  sum)
        colnames(dup_spots_agg) = c(feature_label, "expression")
        
        # there is only one feature associated with the max expression value
        if (length(unique(dup_spots_agg[which(dup_spots_agg$expression == max(dup_spots_agg$expression)), feature_label])) == 1){
          out_df = data.frame(V1=dup_spots_agg[which(dup_spots_agg$expression == max(dup_spots_agg$expression)), feature_label],
                              V2=max(dup_spots_agg$expression),
                              V3=rownames(dup_spots[1,]))
          return(out_df)
        }

        return()
      }
      
    }
    
    out_list = pbmcapply::pbmclapply(spot_id_duplicated, find_rownames_to_keep, mc.cores=32)
    out_df = do.call("rbind", out_list)
    
    spot2feature_1 = df_spot_feature[which(df_spot_feature$spot_id %!in% spot_id_duplicated),]
    spot2feature_2 = df_spot_feature[which(df_spot_feature$spot_id %in% spot_id_duplicated),]
    spot2feature_2 = spot2feature_2[out_df$V3,]
    spot2feature_2[,feature_label] = out_df$V1
    spot2feature_2[,"expression"] = out_df$V2

    spot2feature_roi = rbind(spot2feature_1, spot2feature_2)
    
    # rownames_to_keep = c()
    # for (i in spot_id_duplicated){
    #   dup_spots = df_spot_feature[which(df_spot_feature$spot_id == i),]
    #   
    #   # With mixed cells
    #   if (with_mixed_features == T){
    #     df_spot_feature[rownames(dup_spots[1,]), feature_label] = "Mix"
    #     df_spot_feature[rownames(dup_spots[1,]), "expression"] = mean(dup_spots$expression)
    #     rownames_to_keep = c(rownames_to_keep, rownames(dup_spots[1,]))
    #   } else {
    #     dup_spots_agg = aggregate(dup_spots$expression, 
    #                               list(dup_spots[,feature_label]),
    #                               sum)
    #     colnames(dup_spots_agg) = c(feature_label, "expression")
    #     
    #     # there is only one feature associated with the max expression value
    #     if (length(unique(dup_spots_agg[which(dup_spots_agg$expression == max(dup_spots_agg$expression)), feature_label])) == 1){
    #       df_spot_feature[rownames(dup_spots[1,]), feature_label] = dup_spots_agg[which(dup_spots_agg$expression == max(dup_spots_agg$expression)), feature_label]
    #       df_spot_feature[rownames(dup_spots[1,]), "expression"] = max(dup_spots_agg$expression)
    #       rownames_to_keep = c(rownames_to_keep, rownames(dup_spots[1,]))
    #     }
    #     
    #     # # there is only one cell type associated with the max expression value
    #     # if (length(unique(dup_spots[which(dup_spots$expression == max(dup_spots$expression)), feature_label])) == 1){
    #     #   df_spot_feature[rownames(dup_spots[1,]), feature_label] = dup_spots[which(dup_spots$expression == max(dup_spots$expression)), feature_label][1]
    #     #   df_spot_feature[rownames(dup_spots[1,]), "expression"] = max(dup_spots$expression)
    #     #   rownames_to_keep = c(rownames_to_keep, rownames(dup_spots[1,]))
    #     # }
    #   }
    #   
    # }
    # 
    # spot2feature_1 = df_spot_feature[which(df_spot_feature$spot_id %!in% spot_id_duplicated),]
    # spot2feature_2 = df_spot_feature[which(df_spot_feature$spot_id %in% spot_id_duplicated),]
    # spot2feature_2 = spot2feature_2[rownames_to_keep,]
    # 
    # spot2feature_ = rbind(spot2feature_1, spot2feature_2)
  }
  
  return(spot2feature_roi)
  
}


#### Assign feature to a spot using nearest neighbor #### 
assign_feature_nearest_neighbor <- function(spot_id,
                                            input_data.feature_NN,
                                            spot2feature_roi,
                                            feature_label){
  if (spot_id %in% spot2feature_roi$spot_id){
    out = as.character(spot2feature_roi[spot_id, feature_label])
  } else {
    row_i = input_data.feature_NN[spot_id, "row_i"]
    col_i = input_data.feature_NN[spot_id, "col_i"]
    
    dist_vect = (row_i - spot2feature_roi$row_i)^2 + (col_i - spot2feature_roi$col_i)^2
    out = as.character(spot2feature_roi[rownames(spot2feature_roi)[which(dist_vect == min(dist_vect))], feature_label])[1]
  }
  return(out)
}


#### Test enrichment of spots with feature across tiles (TBD) ####
feature_enrichment_tiles <- function(spot2feature_roi,
                                     feature_label){
  
  # Chi-square between feature and tiles
  out_list = list()
  for (i in sort(unique(spot2feature_roi$tile_id))){
    for (j in sort(as.character(unique(spot2feature_roi$layer)))){
      a = nrow(spot2feature_roi[which(spot2feature_roi$tile_id == i & spot2feature_roi[,feature_label] == j),])
      b = nrow(spot2feature_roi[which(spot2feature_roi$tile_id == i & spot2feature_roi[,feature_label] != j),])
      c = nrow(spot2feature_roi[which(spot2feature_roi$tile_id != i & spot2feature_roi[,feature_label] == j),])
      d = nrow(spot2feature_roi[which(spot2feature_roi$tile_id != i & spot2feature_roi[,feature_label] != j),])
      
      temp = as.table(matrix(data=c(a,b,c,d), nrow=2, ncol=2, byrow = T))
      out = chisq.test(temp)
      df = data.frame(feature = j,
                      tile_id = i,
                      odds_ratio = (a+1)*(d+1)/((b+1)*(c+1)),
                      p_value = out$p.value)
      out_list[[paste0(i,"_",j)]] = df
    }
  }
  
  df = do.call("rbind", out_list)
  
  # domain.cell_type_matrix_OR = reshape2::dcast(df, cell_type ~ spatial_domain, value.var = "odds_ratio", fill = 0)
  # rownames(domain.cell_type_matrix_OR) = domain.cell_type_matrix_OR$cell_type
  # domain.cell_type_matrix_OR = domain.cell_type_matrix_OR[,-1]
  # 
  # domain.cell_type_matrix_pval = reshape2::dcast(df, cell_type ~ spatial_domain, value.var = "p_value", fill = NA)
  # rownames(domain.cell_type_matrix_pval) = domain.cell_type_matrix_pval$cell_type
  # domain.cell_type_matrix_pval = domain.cell_type_matrix_pval[,-1]
  # 
  # temp = round(domain.cell_type_matrix_OR,2) * ((domain.cell_type_matrix_OR > OR_thres) * (domain.cell_type_matrix_pval < p_value_thres))
  # temp[which(temp == 0 | is.na(temp), arr.ind = T)] = ""
  # 
  # wb <- createWorkbook()
  # addWorksheet(wb, "OR")
  # writeData(wb, "OR", domain.cell_type_matrix_OR, startRow = 1, startCol = 1, rowNames = T, colNames = T)
  # addWorksheet(wb, "p_value")
  # writeData(wb, "p_value", domain.cell_type_matrix_pval, startRow = 1, startCol = 1, rowNames = T, colNames = T)
  # addWorksheet(wb, "enrichment")
  # writeData(wb, "enrichment", temp, startRow = 1, startCol = 1, rowNames = T, colNames = T)
  # 
  # saveWorkbook(wb, file = paste0(OUT_DIR, "/spaceflow_", spaceflow_option, "/chisq_domain_cell_type.xlsx"), overwrite = TRUE)
  

}


#### Test statistical association between feature and spot location (TBD) ####
feature_enrichment_spatial_location <- function(spot2feature_roi,
                                                feature_label,
                                                coordinate_field){
  
  spot2feature_roi[,feature_label] = factor(spot2feature_roi[,feature_label])
  
  temp = spot2feature_roi[which(spot2feature_roi$expression > 1),]
  temp$y = temp$row_i
  temp$y = temp$col_i
  temp$y = sqrt(temp$row_i^2 + temp$col_i^2)
  
  
  ggplot(data=temp) +
    geom_violin(aes(x=layer,y=y)) +
    #stat_boxplot(aes(x=layer,y=y), geom ='errorbar', width=0.2) +
    geom_boxplot(aes(x=layer,y=y), width=0.2) +
    xlab("")
  
  library(dplyr)
  group_by(temp, layer) %>%
    summarise(
      count = n(),
      mean = mean(y, na.rm = TRUE),
      sd = sd(y, na.rm = TRUE)
    )
  
  # Compute the analysis of variance
  res.aov <- aov(y ~ layer, data = temp)
  # Summary of the analysis
  summary(res.aov)
  
}





#### Create spot x gene matrix ####  
make_spot_gene_matrix <- function(input_dataframe,
                                  n_HVG=0 # set to 0 to use all the genes, otherwise to compute and use a certain number of high variable genes
)
{
  
  if (n_HVG > 0){
    ### Compute gene variability across spots
    temp1 = aggregate(input_dataframe$expression, list(input_dataframe$gene_name), mean)
    temp2 = aggregate(input_dataframe$expression, list(input_dataframe$gene_name), var)
    
    input_data.filter_stats = merge(temp1, temp2, by = "Group.1")
    colnames(input_data.filter_stats) = c("gene_name", "mean", "var")
    input_data.filter_stats = input_data.filter_stats[order(input_data.filter_stats$var/input_data.filter_stats$mean, decreasing = T),]
    input_data.filter_stats$HVG = 0
    input_data.filter_stats[1:n_HVG,"HVG"] = 1
    input_data.filter_stats$HVG = factor(input_data.filter_stats$HVG)
    
    p <- ggplot(data = input_data.filter_stats) +
      geom_point(aes(x=mean, y=log(var/mean + 1), color=HVG)) +
      scale_color_manual(values = c("black","red"))
    
    # png(paste0(OUT_DIR, "/HVG.png"), width = 5, height = 4, res = 200, units = "in")
    # print(p)
    # dev.off()
    
    # Select the top high variable genes 
    hvgs = input_data.filter_stats[1:n_HVG, "gene_name"]
    input_dataframe = input_dataframe[!duplicated(input_dataframe[,c("spot_id","gene_name")]),] # to remove potential duplicates given by the same gene name but different biotypes
    out_matrix = reshape2::dcast(input_dataframe[which(input_dataframe$gene_name %in% hvgs),], spot_id ~ gene_name, value.var = "expression", fill = 0)
  } else {
    # If using all the genes
    input_dataframe = input_dataframe[!duplicated(input_dataframe[,c("spot_id","gene_name")]),] # to remove potential duplicates given by the same gene name but different biotypes
    out_matrix = reshape2::dcast(input_dataframe, spot_id ~ gene_name, value.var = "expression", fill = 0)
  }
  
  ### Run this for both
  # rownames(out_matrix) = out_matrix[,1]
  # out_matrix = out_matrix[,-1]
  
  # To fix the issue when there is only one gene
  temp_rownames = out_matrix[,1]
  temp_colnames = colnames(out_matrix)[-1]
  out_matrix = data.frame(out_matrix[,-1])
  rownames(out_matrix) = temp_rownames
  colnames(out_matrix) = temp_colnames
  
  return(out_matrix)
}

##### Test function to multicore dcast
# df = data.frame(spot_id = rep(c(1:5),4),
#                 gene_name = paste0("a",1:20),
#                 expression = c(100:119))
# reshape2::dcast(df, spot_id ~ gene_name, value.var = "expression", fill = 0)
# 
# 
# n_sets = 50
# 
# df = input_data.filter_1
# df = df[!duplicated(df[,c("spot_id","gene_name")]),]
# 
# my_spots = sort(unique(df$spot_id))
# my_genes = sort(unique(df$gene_name))
# 
# # Split spots into sets
# my_spots_split = split(my_spots, ceiling(seq_along(my_spots) / round(length(my_spots)/n_sets)))
# my_genes_split = split(my_genes, ceiling(seq_along(my_genes) / round(length(my_genes)/n_sets)))
# 
# dcast_indexes = expand.grid(1:length(my_spots_split), 1:length(my_genes_split))
# 
# my_dcast <- function(x){ # x is a row of dcast index
#   df_temp = df[which(df$spot_id %in% my_spots_split[dcast_indexes[x,1]][[1]] &
#                        df$gene_name %in% my_genes_split[dcast_indexes[x,2]][[1]]),]
#   
#   temp_mat = matrix(0, nrow=length(my_spots_split[dcast_indexes[x,1]][[1]]), ncol = length(my_genes_split[dcast_indexes[x,2]][[1]]))
#   rownames(temp_mat) = my_spots_split[dcast_indexes[x,1]][[1]]
#   colnames(temp_mat) = my_genes_split[dcast_indexes[x,2]][[1]]
#   
#   if (nrow(df_temp) > 0){
#     out_dcast = reshape2::dcast(df_temp, spot_id ~ gene_name, value.var = "expression", fill = 0)
#     temp_mat[as.character(out_dcast$spot_id), setdiff(colnames(out_dcast), "spot_id")] = as.matrix(out_dcast[,-1])
#     # for (n in 1:nrow(out_dcast)){
#     #   for (i in setdiff(colnames(out_dcast), "spot_id")){
#     #     temp_mat[as.character(out_dcast[n,"spot_id"]), i] = out_dcast[n,i]
#     #   }
#     # }
#   }
#   return(temp_mat)
# }
# 
# lapply(1:nrow(dcast_indexes), my_dcast)
# 
# out_list = pbmcapply::pbmclapply(1:nrow(dcast_indexes), 
#                                  my_dcast,
#                                  mc.cores = 32,
#                                  mc.cleanup = T)
# 
# out_matrix = matrix(0, nrow=length(my_spots), ncol=length(my_genes))
# rownames(out_matrix) = my_spots
# colnames(out_matrix) = my_genes
# 
# for (i in 1:length(out_list)){
#   out_matrix[rownames(out_list[[i]]), colnames(out_list[[i]])] = out_list[[i]]
# }



#### Create Seurat input files (TBD) ####
# create_seurat_input_files <- function(input_data,
#                                       annotation,
#                                       out_dir,
#                                       SAMPLE_NAME,
#                                       ROI_label){
#   
#   dir.create(paste0(INPUT_DIR, "/", SAMPLE_NAME, "/", ROI_label, "/feature_bc_matrix"), showWarnings = F)
#   
#   # features
#   input_data$gene_id = annotation[input_data$gene_name, "gene_id"]
#   input_data$feature = paste0(input_data$gene_id, "_", input_data$gene_name)
#   
#   # Save seurat files for the sample
#   barcodes = sort(unique(out$barcode))
#   write.table(barcodes, paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/MUSIC/data/feature_bc_matrix/", i, "/barcodes.tsv"), row.names = F, col.names = F, quote = F)
#   
#   temp = out[,c("gene_id","gene_name")]
#   temp = temp[!duplicated(temp),]
#   temp = temp[order(temp$gene_name),]
#   features = paste0(temp$gene_id, "_", temp$gene_name)
#   write.table(data.frame(V1 = temp$gene_id,
#                          V2 = temp$gene_name,
#                          V3 = "Gene Expression"),
#               paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/MUSIC/data/feature_bc_matrix/", i, "/features.tsv"), row.names = F, col.names = F, quote = F, sep = "\t")
#   
#   expr_matrix = data.frame(V1 = as.numeric(sapply(out$feature, function(x){which(features == x)})),
#                            V2 = as.numeric(sapply(out$barcode, function(x){which(barcodes == x)})),
#                            V3 = out$expression)
#   expr_matrix_mtx = Matrix::sparseMatrix(i = expr_matrix$V1, j = expr_matrix$V2, x = expr_matrix$V3)
#   Matrix::writeMM(obj = expr_matrix_mtx, paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/MUSIC/data/feature_bc_matrix/", i, "/matrix.mtx"))
#   
#   out_list[[i]] = out
#   
# }

#### Create SpaceFlow input data #### 
# spaceflow_option 1: centroids and marker genes are used
# spaceflow_option 2: centroids and all the expressed genes in centroids are used
# spaceflow_option 3: all the spots and all the expressed genes within the area delimited by rect (to reduce computational cost)
# input_data is used only in option 3, same for the rect coordinates
create_spaceflow_input <- function(input_data,
                                   input_data.filter,
                                   df_markers,
                                   feature_label,
                                   spaceflow_option,
                                   min_spot, # Include genes detected in at least this many spots.
                                   min_gene, # Include spots where at least this many genes are detected.
                                   xmin_rect,
                                   xmax_rect,
                                   ymin_rect,
                                   ymax_rect){
  
  if (spaceflow_option == 1){
    dir.create(file.path(paste0(OUT_DIR, "/spaceflow_1")), showWarnings = FALSE)
    
    input_data.filter.temp = input_data.filter[which(input_data.filter$gene_name %in% unique(df_markers$gene)),]
    input_data.filter.temp_matrix = make_spot_gene_matrix(input_data.filter.temp, n_HVG=0)
    
    temp = as.numeric(as.matrix(input_data.filter.temp_matrix))
    spaceflow_matrix = input_data.filter.temp_matrix * 1/min(temp[temp!=0])
    write.table(spaceflow_matrix, paste0(OUT_DIR,"/spaceflow_1/spaceflow_matrix.txt"), sep = "\t", row.names = T, col.names = T, quote = F)
    
    temp = input_data.filter.temp[!duplicated(input_data.filter.temp$spot_id), c("spot_id", "row_i", "col_i")]
    rownames(temp) = temp$spot_id
    temp = temp[rownames(spaceflow_matrix),c("col_i","row_i")]
    write.table(temp, paste0(OUT_DIR,"/spaceflow_1/spaceflow_coord.txt"), sep = "\t", row.names = F, col.names = T, quote = F)
    
    temp = input_data.filter.temp[!duplicated(input_data.filter.temp$spot_id), c("spot_id", feature_label)]
    rownames(temp) = temp$spot_id
    write.table(temp[rownames(spaceflow_matrix),feature_label], paste0(OUT_DIR,"/spaceflow_1/spaceflow_", feature_label, ".txt"), sep = "\t", row.names = F, col.names = T, quote = F)
  } 
  
  else if (spaceflow_option == 2){
    dir.create(file.path(paste0(OUT_DIR, "/spaceflow_2")), showWarnings = FALSE)
    
    input_data.filter.temp = input_data.filter
    
    if (!is.na(min_spot)){
      tab_genes = table(as.factor(input_data.filter$gene_name))
      input_data.filter.temp = input_data.filter.temp[which(input_data.filter.temp$gene_name %in% names(tab_genes[which(tab_genes >= 100)]),]
    }
    
    if (!is.na(min_gene)){
      length(unique(input_data.filter$spot_id))
    }
    
    input_data.filter.temp_matrix = make_spot_gene_matrix(input_data.filter.temp, n_HVG=0)
    
    temp = as.numeric(as.matrix(input_data.filter.temp_matrix))
    spaceflow_matrix = input_data.filter.temp_matrix * 1/min(temp[temp!=0])
    write.table(spaceflow_matrix, paste0(OUT_DIR,"/spaceflow_2/spaceflow_matrix.txt"), sep = "\t", row.names = T, col.names = T, quote = F)
    
    temp = input_data.filter.temp[!duplicated(input_data.filter.temp$spot_id), c("spot_id", "row_i", "col_i")]
    rownames(temp) = temp$spot_id
    temp = temp[rownames(spaceflow_matrix),c("col_i","row_i")]
    write.table(temp, paste0(OUT_DIR,"/spaceflow_2/spaceflow_coord.txt"), sep = "\t", row.names = F, col.names = T, quote = F)
    
    temp = input_data.filter.temp[!duplicated(input_data.filter.temp$spot_id), c("spot_id", feature_label)]
    rownames(temp) = temp$spot_id
    write.table(temp[rownames(spaceflow_matrix), feature_label], paste0(OUT_DIR,"/spaceflow_2/spaceflow_", feature_label, ".txt"), sep = "\t", row.names = F, col.names = T, quote = F)
  }
  
  else if (spaceflow_option == 3){
    dir.create(file.path(paste0(OUT_DIR, "/spaceflow_3")), showWarnings = FALSE)
    
    input_data.filter.temp = input_data[which(input_data$row_i >= xmin_rect & input_data$row_i <= xmax_rect & input_data$col_i >= ymin_rect & input_data$col_i <= ymax_rect),]
    input_data.filter.temp_matrix = make_spot_gene_matrix(input_data.filter.temp, n_HVG=0)
    
    temp = as.numeric(as.matrix(input_data.filter.temp_matrix))
    spaceflow_matrix = input_data.filter.temp_matrix * 1/min(temp[temp!=0])
    write.table(spaceflow_matrix, paste0(OUT_DIR,"/spaceflow_3/spaceflow_matrix.txt"), sep = "\t", row.names = T, col.names = T, quote = F)
    
    temp = input_data.filter.temp[!duplicated(input_data.filter.temp$spot_id), c("spot_id", "row_i", "col_i")]
    rownames(temp) = temp$spot_id
    temp = temp[rownames(spaceflow_matrix),c("col_i","row_i")]
    write.table(temp, paste0(OUT_DIR,"/spaceflow_3/spaceflow_coord.txt"), sep = "\t", row.names = F, col.names = T, quote = F)
    
  }
  
  return()
  
}


#### Read SpaceFlow output #### 
read_spaceflow_output <- function(input_data,
                                  input_data.filter,
                                  df_markers,
                                  feature_label,
                                  spaceflow_option,
                                  OUT_DIR,
                                  xmin_rect,
                                  xmax_rect,
                                  ymin_rect,
                                  ymax_rect){
  
  spaceflow_matrix = data.frame(fread(paste0(OUT_DIR,"/spaceflow_", spaceflow_option, "/spaceflow_matrix.txt"), sep = "\t", nThread = 16))
  rownames(spaceflow_matrix) = spaceflow_matrix$V1
  spaceflow_matrix = spaceflow_matrix[,-1]
  
  if (spaceflow_option == 1){
    input_data.filter.temp = input_data.filter[which(input_data.filter$gene_name %in% unique(df_markers$gene)),]
    
    df_spatial_domains = input_data.filter.temp[!duplicated(input_data.filter.temp$spot_id), c("spot_id", "row_i", "col_i", feature_label)]
    rownames(df_spatial_domains) = df_spatial_domains$spot_id
    df_spatial_domains = df_spatial_domains[rownames(spaceflow_matrix),]
    df_spatial_domains$spatial_domain = factor(read.table(paste0(OUT_DIR,"/spaceflow_", spaceflow_option, "/domains.tsv"))$V1)
  }
  
  else if (spaceflow_option == 2){
    input_data.filter.temp = input_data.filter
    
    df_spatial_domains = input_data.filter.temp[!duplicated(input_data.filter.temp$spot_id), c("spot_id", "row_i", "col_i", feature_label)]
    rownames(df_spatial_domains) = df_spatial_domains$spot_id
    df_spatial_domains = df_spatial_domains[rownames(spaceflow_matrix),]
    df_spatial_domains$spatial_domain = factor(read.table(paste0(OUT_DIR,"/spaceflow_", spaceflow_option, "/domains.tsv"))$V1)
  }
  
  else if (spaceflow_option == 3){
    input_data.filter.temp = input_data[which(input_data$row_i >= xmin_rect & input_data$row_i <= xmax_rect & input_data$col_i >= ymin_rect & input_data$col_i <= ymax_rect),]
    
    df_spatial_domains = input_data.filter.temp[!duplicated(input_data.filter.temp$spot_id), c("spot_id", "row_i", "col_i")]
    df_spatial_domains = merge(df_spatial_domains, spot2feature_roi[,c("spot_id", feature_label)], by = "spot_id", all.x = T)
    df_spatial_domains[which(is.na(df_spatial_domains[,feature_label])), feature_label] = "NA"
    rownames(df_spatial_domains) = df_spatial_domains$spot_id
    df_spatial_domains = df_spatial_domains[rownames(spaceflow_matrix),]
    df_spatial_domains$spatial_domain = factor(read.table(paste0(OUT_DIR,"/spaceflow_", spaceflow_option, "/domains.tsv"))$V1)
  }
  
  return(df_spatial_domains)
  
}


#### Assign feature (cell_type, layer, region) to spatial domains #### 
assign_feature_spatial_domains <- function(df_spatial_domains,
                                           spaceflow_option,
                                           OUT_DIR,
                                           OR_thres=1,
                                           p_value_thres=0.01,
                                           feature # cell_type, layer, region
){
  
  # Chi-square between features and spatial domains
  out_list = list()
  for (i in 0:(length(unique(df_spatial_domains$spatial_domain))-1)){
    for (j in sort(as.character(unique(df_spatial_domains[,feature])))){
      a = nrow(df_spatial_domains[which(df_spatial_domains$spatial_domain == i & df_spatial_domains[,feature] == j),])
      b = nrow(df_spatial_domains[which(df_spatial_domains$spatial_domain == i & df_spatial_domains[,feature] != j),])
      c = nrow(df_spatial_domains[which(df_spatial_domains$spatial_domain != i & df_spatial_domains[,feature] == j),])
      d = nrow(df_spatial_domains[which(df_spatial_domains$spatial_domain != i & df_spatial_domains[,feature] != j),])
      
      temp = as.table(matrix(data=c(a,b,c,d), nrow=2, ncol=2, byrow = T))
      out = chisq.test(temp)
      df = data.frame(feature = j,
                      spatial_domain = i,
                      odds_ratio = (a+1)*(d+1)/((b+1)*(c+1)),
                      p_value = out$p.value)
      colnames(df)[1] = feature
      out_list[[paste0(i,"_",j)]] = df
    }
  }
  
  df = do.call("rbind", out_list)
  
  
  domain.feature_matrix_OR = reshape2::dcast(df, as.formula(paste0(feature, " ~ spatial_domain")), value.var = "odds_ratio", fill = 0)
  rownames(domain.feature_matrix_OR) = domain.feature_matrix_OR[,feature]
  domain.feature_matrix_OR = domain.feature_matrix_OR[,-1]
  
  domain.feature_matrix_pval = reshape2::dcast(df, as.formula(paste0(feature, " ~ spatial_domain")), value.var = "p_value", fill = NA)
  rownames(domain.feature_matrix_pval) = domain.feature_matrix_pval[,feature]
  domain.feature_matrix_pval = domain.feature_matrix_pval[,-1]
  
  temp = round(domain.feature_matrix_OR,2) * ((domain.feature_matrix_OR > OR_thres) * (domain.feature_matrix_pval < p_value_thres))
  temp[which(temp == 0 | is.na(temp), arr.ind = T)] = ""
  
  wb <- createWorkbook()
  addWorksheet(wb, "OR")
  writeData(wb, "OR", domain.feature_matrix_OR, startRow = 1, startCol = 1, rowNames = T, colNames = T)
  addWorksheet(wb, "p_value")
  writeData(wb, "p_value", domain.feature_matrix_pval, startRow = 1, startCol = 1, rowNames = T, colNames = T)
  addWorksheet(wb, "enrichment")
  writeData(wb, "enrichment", temp, startRow = 1, startCol = 1, rowNames = T, colNames = T)
  
  saveWorkbook(wb, file = paste0(OUT_DIR, "/spaceflow_", spaceflow_option, "/chisq_domain_", feature, ".xlsx"), overwrite = TRUE)
  
}

#### Assign layers to spatial domains (OBSOLETE) #### 
# assign_layers_spatial_domains <- function(df_spatial_domains,
#                                           spaceflow_option,
#                                           OUT_DIR,
#                                           OR_thres=1,
#                                           p_value_thres=0.01){
#   
#   # Chi-square between cell types and spatial domains
#   out_list = list()
#   for (i in 0:(length(unique(df_spatial_domains$spatial_domain))-1)){
#     for (j in sort(as.character(unique(df_spatial_domains$layer)))){
#       a = nrow(df_spatial_domains[which(df_spatial_domains$spatial_domain == i & df_spatial_domains$layer == j),])
#       b = nrow(df_spatial_domains[which(df_spatial_domains$spatial_domain == i & df_spatial_domains$layer != j),])
#       c = nrow(df_spatial_domains[which(df_spatial_domains$spatial_domain != i & df_spatial_domains$layer == j),])
#       d = nrow(df_spatial_domains[which(df_spatial_domains$spatial_domain != i & df_spatial_domains$layer != j),])
#       
#       temp = as.table(matrix(data=c(a,b,c,d), nrow=2, ncol=2, byrow = T))
#       out = chisq.test(temp)
#       df = data.frame(layer = j,
#                       spatial_domain = i,
#                       odds_ratio = (a+1)*(d+1)/((b+1)*(c+1)),
#                       p_value = out$p.value)
#       out_list[[paste0(i,"_",j)]] = df
#     }
#   }
#   
#   df = do.call("rbind", out_list)
#   
#   domain.layer_matrix_OR = reshape2::dcast(df, layer ~ spatial_domain, value.var = "odds_ratio", fill = 0)
#   rownames(domain.layer_matrix_OR) = domain.layer_matrix_OR$layer
#   domain.layer_matrix_OR = domain.layer_matrix_OR[,-1]
#   
#   domain.layer_matrix_pval = reshape2::dcast(df, layer ~ spatial_domain, value.var = "p_value", fill = NA)
#   rownames(domain.layer_matrix_pval) = domain.layer_matrix_pval$layer
#   domain.layer_matrix_pval = domain.layer_matrix_pval[,-1]
#   
#   temp = round(domain.layer_matrix_OR,2) * ((domain.layer_matrix_OR > OR_thres) * (domain.layer_matrix_pval < p_value_thres))
#   temp[which(temp == 0 | is.na(temp), arr.ind = T)] = ""
#   
#   wb <- createWorkbook()
#   addWorksheet(wb, "OR")
#   writeData(wb, "OR", domain.layer_matrix_OR, startRow = 1, startCol = 1, rowNames = T, colNames = T)
#   addWorksheet(wb, "p_value")
#   writeData(wb, "p_value", domain.layer_matrix_pval, startRow = 1, startCol = 1, rowNames = T, colNames = T)
#   addWorksheet(wb, "enrichment")
#   writeData(wb, "enrichment", temp, startRow = 1, startCol = 1, rowNames = T, colNames = T)
#   
#   saveWorkbook(wb, file = paste0(OUT_DIR, "/spaceflow_", spaceflow_option, "/chisq_domain_layer.xlsx"), overwrite = TRUE)
#   
# }


#### Dotplot with cluster-feature enrichment and significance #### 
plot_cluster_feature_dotplot <- function(OUT_DIR,
                                         spaceflow_option,
                                         feature_label,
                                         p_value_thres=0.01,
                                         p_value_cutoff=NA){
  
  domain.feature_matrix_OR = read.xlsx(paste0(OUT_DIR, "/spaceflow_", spaceflow_option, "/chisq_domain_", feature_label, ".xlsx"), rowNames = T, sheet = "OR")
  domain.feature_matrix_pval = read.xlsx(paste0(OUT_DIR, "/spaceflow_", spaceflow_option, "/chisq_domain_", feature_label, ".xlsx"), rowNames = T, sheet = "p_value")
  
  # col_order = c(1,6,11,5,7,8,4,2,9,10,3,12)
  # domain.feature_matrix_OR = domain.feature_matrix_OR[,col_order]
  # domain.feature_matrix_pval = domain.feature_matrix_pval[,col_order]
  # colnames(domain.feature_matrix_OR) = c(0:11)
  # colnames(domain.feature_matrix_pval) = c(0:11)
  
  df = reshape2::melt(as.matrix(domain.feature_matrix_OR))
  colnames(df) = c(feature_label, "cluster", "odds_ratio")
  df = cbind(df, reshape2::melt(as.matrix(domain.feature_matrix_pval))[,3])
  colnames(df)[4] = "p_val"
  df[,feature_label] = factor(df[,feature_label], levels = unique(df[,feature_label])[order(unique(df[,feature_label]), decreasing = T)])
  
  # To avoid representing very small p-values
  if (!is.na(p_value_cutoff)){
    df$p_val[df$p_val < p_value_cutoff] = p_value_cutoff
  }
  
  
  p<-ggplot(df, aes_string(x="cluster", y=feature_label)) +
    geom_point(aes(size = odds_ratio, color = -log10(p_val))) +
    scale_x_continuous(breaks = as.numeric(colnames(domain.feature_matrix_OR)), labels = colnames(domain.feature_matrix_OR)) +
    scale_size_continuous(breaks = c(0, 1, 5, 10), labels = c("0","1","5","10"), range = c(0,15)) +
    scale_colour_gradient2(low = "blue", mid = "yellow", high="red", midpoint = -log10(p_value_thres)) +
    labs(size="Odds ratio", colour="-log10(p_val)") +
    xlab("Spatial cluster") + ylab(feature_label) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=10),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          plot.title = element_text(hjust = 0.5),
          legend.position = "bottom")
  
  png(paste0(OUT_DIR, "/spaceflow_", spaceflow_option, "/spatial_domain_", feature_label, "_dotplot.png"), width = 7, height = 3, res = 400, units = "in")
  print(p)
  dev.off()
}


#### Plot spatial domains #### 
plot_spatial_domains <- function(df_spatial_domains,
                                 OUT_DIR,
                                 out_filename="spot_spatial_domain_map.png",
                                 spaceflow_option,
                                 feature_label,
                                 only_cluster_assigned_to_feature,
                                 spatial_domains_updated_names=c(),
                                 cluster_colors=c(),
                                 xmin_rect,
                                 xmax_rect,
                                 ymin_rect,
                                 ymax_rect){
  
  domain.feature_matrix_enrich = read.xlsx(paste0(OUT_DIR, "/spaceflow_", spaceflow_option, "/chisq_domain_", feature_label, ".xlsx"), rowNames = T, sheet = "enrichment")
  
  create_legend_labels <- function(x){
    out = as.character(paste(rownames(domain.feature_matrix_enrich)[which(domain.feature_matrix_enrich[,x] != "")], collapse = ","))
    if (out == ""){
      out = "NA"
    }
    return(out)
  }
  
  ### Default clusters
  if (length(spatial_domains_updated_names) == 0){
    
    legend_labels = paste(as.numeric(colnames(domain.feature_matrix_enrich)),
                          as.character(sapply(colnames(domain.feature_matrix_enrich), create_legend_labels)),
                          sep = "-")
    
    # Remove clusters not assigned to feature
    if (only_cluster_assigned_to_feature == T){
      temp_label = which(sapply(legend_labels, function(x){strsplit(x,"-")[[1]][2]}) != "NA")
      legend_labels = legend_labels[temp_label]
      temp_cluster = sapply(legend_labels, function(x){strsplit(x,"-")[[1]][1]})
    } else {
      temp_cluster = unique(df_spatial_domains$spatial_domain)
    }

    if (spaceflow_option == 1 | spaceflow_option == 2){
      p<-ggplot(data=df_spatial_domains[which(df_spatial_domains$spatial_domain %in% temp_cluster),]) +
        geom_point(aes(x=col_i-min(col_i), y=max(row_i)-row_i, color=spatial_domain), size = 0.5) +
        {if(length(cluster_colors)>0)scale_color_manual(values = cluster_colors, labels=legend_labels)} +
        {if(length(cluster_colors)==0)scale_color_discrete(labels=legend_labels)} +
        scale_x_continuous(breaks=breaks_x,
                           labels=paste0(as.character(labels_x), " um")) +
        scale_y_continuous(breaks=breaks_y,
                           labels=paste0(as.character(labels_y), " um")) +
        xlab("") +
        ylab("") +
        guides(colour = guide_legend(override.aes = list(size=2))) +
        theme_bw() +
        theme(legend.title = element_blank())
      
      p1 <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_rect(fill='transparent'),
                      axis.text.x = element_blank(),
                      axis.text.y = element_blank(),
                      axis.ticks = element_blank(),
                      plot.background = element_rect(fill='transparent', color=NA),
                      legend.position = "none")
      
    } else if (spaceflow_option == 3){
      breaks_3_x = seq(0, xmax_rect-xmin_rect, length=4)
      labels_3_x = round(breaks_3_x * 0.015) + xmin_rect * 0.015
      
      breaks_3_y = seq(0, ymax_rect-ymin_rect, length=4)
      labels_3_y = round(breaks_3_y * 0.015) + ymin_rect * 0.015
      
      p<-ggplot(data=df_spatial_domains[which(df_spatial_domains$spatial_domain %in% temp_cluster),]) +
        geom_point(aes(x=col_i-min(col_i), y=max(row_i)-row_i, color=spatial_domain), size = 0.5) +
        {if(length(cluster_colors)>0)scale_color_manual(values = cluster_colors, labels=legend_labels)} +
        {if(length(cluster_colors)==0)scale_color_discrete(labels=legend_labels)} +
        scale_x_continuous(breaks=breaks_3_x,
                           labels=paste0(as.character(labels_3_x), " um")) +
        scale_y_continuous(breaks=breaks_3_y,
                           labels=paste0(as.character(labels_3_y), " um")) +
        xlab("") +
        ylab("") +
        guides(colour = guide_legend(override.aes = list(size=2))) +
        theme_bw() +
        theme(legend.title = element_blank())

      p1 <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_rect(fill='transparent'),
                      axis.text.x = element_blank(),
                      axis.text.y = element_blank(),
                      axis.ticks = element_blank(),
                      plot.background = element_rect(fill='transparent', color=NA),
                      legend.position = "none")
    }
  
  }
  
  ### Custom cluster order
  else if (length(spatial_domains_updated_names) > 0){
    spatial_domains_updated = c(0:(length(unique(df_spatial_domains$spatial_domain))-1))
    names(spatial_domains_updated) = as.character(spatial_domains_updated_names)
    df_spatial_domains$spatial_domain_updated = factor(as.character(sapply(as.character(df_spatial_domains$spatial_domain), function(x){spatial_domains_updated[x]})), levels = c(0:(length(unique(df_spatial_domains$spatial_domain))-1)))
    
    legend_labels = paste(as.numeric(colnames(domain.feature_matrix_enrich)),
                          as.character(sapply(names(spatial_domains_updated), create_legend_labels)),
                          sep = "-")
    
    if (spaceflow_option == 1 | spaceflow_option == 2){
      p<-ggplot(data=df_spatial_domains) +
        geom_point(aes(x=col_i-min(col_i), y=max(row_i)-row_i, color=spatial_domain_updated), size = 0.5) +
        {if(length(cluster_colors)>0)scale_color_manual(values = cluster_colors, labels=legend_labels)} +
        {if(length(cluster_colors)==0)scale_color_discrete(labels=legend_labels)} +
        scale_x_continuous(breaks=breaks_x,
                           labels=paste0(as.character(labels_x), " um")) +
        scale_y_continuous(breaks=breaks_y,
                           labels=paste0(as.character(labels_y), " um")) +
        xlab("") +
        ylab("") +
        guides(colour = guide_legend(override.aes = list(size=2))) +
        theme_bw() +
        theme(legend.title = element_blank())
      
      p1 <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_rect(fill='transparent'),
                      axis.text.x = element_blank(),
                      axis.text.y = element_blank(),
                      axis.ticks = element_blank(),
                      plot.background = element_rect(fill='transparent', color=NA),
                      legend.position = "none")
      
    } else if (spaceflow_option == 3){
      breaks_3_x = seq(0, xmax_rect-xmin_rect, length=4)
      labels_3_x = round(breaks_3_x * 0.015) + xmin_rect * 0.015
      
      breaks_3_y = seq(0, ymax_rect-ymin_rect, length=4)
      labels_3_y = round(breaks_3_y * 0.015) + ymin_rect * 0.015
      
      p<-ggplot(data=df_spatial_domains) +
        geom_point(aes(x=col_i-min(col_i), y=max(row_i)-row_i, color=spatial_domain_updated), size = 0.5) +
        {if(length(cluster_colors)>0)scale_color_manual(values = cluster_colors, labels=legend_labels)} +
        {if(length(cluster_colors)==0)scale_color_discrete(labels=legend_labels)} +
        scale_x_continuous(breaks=breaks_3_x,
                           labels=paste0(as.character(labels_3_x), " um")) +
        scale_y_continuous(breaks=breaks_3_y,
                           labels=paste0(as.character(labels_3_y), " um")) +
        xlab("") +
        ylab("") +
        guides(colour = guide_legend(override.aes = list(size=2))) +
        theme_bw() +
        theme(legend.title = element_blank())
      
      p1 <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_rect(fill='transparent'),
                      axis.text.x = element_blank(),
                      axis.text.y = element_blank(),
                      axis.ticks = element_blank(),
                      plot.background = element_rect(fill='transparent', color=NA),
                      legend.position = "none")

    }
    
  }
  
  ggsave(p1, width = 3.7 * n_tiles_col_ROI, height = 5 * n_tiles_row_ROI, dpi = 400, units = "in",
         filename=paste0(OUT_DIR, "/spaceflow_", spaceflow_option, "/transparent_", out_filename),
         bg="transparent")
  
  png(paste0(OUT_DIR, "/spaceflow_", spaceflow_option, "/", out_filename), width = 7, height = 6, res = 400, units = "in")
  print(p)
  dev.off()
  
  return(df_spatial_domains)
  
}


#### Assign domain to a spot using nearest neighbor #### 
assign_spatial_domain_nearest_neighbor <- function(input_data,
                                                   df_spatial_domains,
                                                   spatial_domains_updated_names=c()){
  
  input_data.spatial = input_data[, c(1,2,3,4,8,9)]
  input_data.spatial = input_data.spatial[!duplicated(input_data.spatial$spot_id),]
  rownames(input_data.spatial) = input_data.spatial$spot_id
  
  assign_spatial_domain_spot <- function(spot_id){
    if (spot_id %in% df_spatial_domains$spot_id){
      out = as.character(df_spatial_domains[spot_id, "spatial_domain"])
    } else {
      row_i = input_data.spatial[spot_id, "row_i"]
      col_i = input_data.spatial[spot_id, "col_i"]
      
      dist_vect = (row_i - df_spatial_domains$row_i)^2 + (col_i - df_spatial_domains$col_i)^2
      out = as.character(df_spatial_domains[rownames(df_spatial_domains)[which(dist_vect == min(dist_vect))], "spatial_domain"])[1]
    }
    return(out)
  }
  
  out_list = pbmcapply::pbmclapply(rownames(input_data.spatial), assign_spatial_domain_spot, mc.cores = 32)
  input_data.spatial$spatial_domain = unlist(out_list)
  
  if (length(spatial_domains_updated_names) > 0){
    spatial_domains_updated = c(0:(length(unique(df_spatial_domains$spatial_domain))-1))
    names(spatial_domains_updated) = as.character(spatial_domains_updated_names)
    input_data.spatial$spatial_domain_updated = factor(as.character(sapply(as.character(input_data.spatial$spatial_domain), function(x){spatial_domains_updated[x]})), levels = c(0:(length(unique(input_data.spatial$spatial_domain))-1)))
    
  }

  return(input_data.spatial)
}



#### Perform DE analysis between two samples by cell types ####

# To use this function, "input_data.filter" from cell type analysis is necessary and saved to rds file
# DE analysis is sample 1 vs sample 2 (i.e. log2FC > 0 if 1 > 2)
# sample_1 and sample_2 are sample IDs
# label_1 and label_2 are custom labels, for example conditions (normal and AD)
perform_DE_analysis <- function(main_dir,
                                OUT_DIR,
                                sample_1,
                                sample_2,
                                label_1,
                                label_2,
                                ROI_label_1,
                                ROI_label_2,
                                tiles_1=NA,
                                tiles_2=NA,
                                pct_thres=0.05
)

{
  
  dir.create(OUT_DIR, recursive = T, showWarnings = F)
  
  ### cell type-spot stats
  df1 = read.xlsx(paste0(main_dir, "/", sample_1, "/", ROI_label_1, "/cell_types/spot2cell_table.xlsx"))
  df2 = read.xlsx(paste0(main_dir, "/", sample_2, "/", ROI_label_2, "/cell_types/spot2cell_table.xlsx"))
  
  print(wilcox.test(df1$perc_spots[1:(nrow(df1)-1)], df2$perc_spots[1:(nrow(df2)-1)], paired = T))
  
  print("Read input data")
  input_data.filter_1 = readRDS(paste0(main_dir, "/", sample_1, "/", ROI_label_1, "/cell_types/input_data.filter.rds"))
  input_data.filter_1$spot_id = paste0(input_data.filter_1$spot_id, ".", label_1)
  if (!is.na(tiles_1)){
    input_data.filter_1 = input_data.filter_1[which(input_data.filter_1$tile_id %in% tiles_1),]
  }
  print(paste0("Number of spots sample 1: ", length(unique(input_data.filter_1$spot_id))))
  
  input_data.filter_2 = readRDS(paste0(main_dir, "/", sample_2, "/", ROI_label_2, "/cell_types/input_data.filter.rds"))
  input_data.filter_2$spot_id = paste0(input_data.filter_2$spot_id, ".", label_2)
  if (!is.na(tiles_2)){
    input_data.filter_2 = input_data.filter_2[which(input_data.filter_2$tile_id %in% tiles_2),]
  }
  print(paste0("Number of spots sample 2: ", length(unique(input_data.filter_2$spot_id))))
  
  if (!is.na(pct_thres)){
    temp_genes_1 = names(which(table(as.factor(input_data.filter_1$gene_name)) / length(unique(input_data.filter_1$spot_id)) > pct_thres))
    temp_genes_2 = names(which(table(as.factor(input_data.filter_2$gene_name)) / length(unique(input_data.filter_2$spot_id)) > pct_thres))
    temp_genes = union(temp_genes_1, temp_genes_2)
    
    input_data.filter_1 = input_data.filter_1[which(input_data.filter_1$gene_name %in% temp_genes),]
    input_data.filter_2 = input_data.filter_2[which(input_data.filter_2$gene_name %in% temp_genes),]
  }
  
  
  print("Make spot by gene matrices")
  input_data.filter_1_matrix = make_spot_gene_matrix(input_data.filter_1, n_HVG=0)
  input_data.filter_2_matrix = make_spot_gene_matrix(input_data.filter_2, n_HVG=0)
  
  # Transform the two matrices to have the same genes
  union_genes = union(colnames(input_data.filter_1_matrix), colnames(input_data.filter_2_matrix))
  print(paste0("Number of genes under DE analysis: ", length(union_genes)))
  
  temp_mat = matrix(data = 0, nrow = nrow(input_data.filter_1_matrix), ncol = length(setdiff(union_genes, colnames(input_data.filter_1_matrix))))
  rownames(temp_mat) = rownames(input_data.filter_1_matrix)
  colnames(temp_mat) = setdiff(union_genes, colnames(input_data.filter_1_matrix))
  input_data.filter_1_matrix = cbind(input_data.filter_1_matrix, temp_mat)
  
  temp_mat = matrix(data = 0, nrow = nrow(input_data.filter_2_matrix), ncol = length(setdiff(union_genes, colnames(input_data.filter_2_matrix))))
  rownames(temp_mat) = rownames(input_data.filter_2_matrix)
  colnames(temp_mat) = setdiff(union_genes, colnames(input_data.filter_2_matrix))
  input_data.filter_2_matrix = cbind(input_data.filter_2_matrix, temp_mat)
  
  # Normalize data
  input_data.filter_1_matrix_norm = log(10000 * input_data.filter_1_matrix / rowSums(input_data.filter_1_matrix) + 1)
  input_data.filter_2_matrix_norm = log(10000 * input_data.filter_2_matrix / rowSums(input_data.filter_2_matrix) + 1)
  
  
  print("Perform DE analysis")
  out_list = list()
  for (i in sort(unique(input_data.filter_1$cell_type))){
    
    print(i)
    
    spot_1 = input_data.filter_1[which(input_data.filter_1$cell_type == i), "spot_id"]
    spot_2 = input_data.filter_2[which(input_data.filter_2$cell_type == i), "spot_id"]
    
    expr_1 = log(colMeans(expm1(input_data.filter_1_matrix_norm[spot_1,]))+1)
    expr_2 = log(colMeans(expm1(input_data.filter_2_matrix_norm[spot_2,]))+1)
    
    pct_1 = colSums(input_data.filter_1_matrix_norm[spot_1,]>0)/length(spot_1)
    pct_2 = colSums(input_data.filter_2_matrix_norm[spot_2,]>0)/length(spot_2)
    
    out = data.frame(expr_1, expr_2, pct_1, pct_2)
    out = out[order(out$expr_1, decreasing = T),]
    colnames(out) = paste(c("expr","expr","pct","pct"), rep(c(label_1,label_2),2), sep = ".")
    
    # DE analysis
    out$avg_log2FC = out[,paste0("expr.", label_1)] - out[,paste0("expr.", label_2)]
    out$p_val = unlist(pbmcapply::pbmclapply(rownames(out), function(x){wilcox.test(input_data.filter_1_matrix_norm[,x], input_data.filter_2_matrix_norm[,x])$p.value}, mc.cores = 32))
    out$p_val_adj = p.adjust(out$p_val, method = "bonferroni")
    out_list[[i]] = out
  }
  
  wb <- createWorkbook()
  for (i in sort(unique(input_data.filter_1$cell_type))){
    addWorksheet(wb, i)
    writeData(wb, i, out_list[[i]], startRow = 1, startCol = 1, rowNames = T, colNames = T)
  }
  saveWorkbook(wb, file = paste0(OUT_DIR,"/DE_analysis_cell_types.xlsx"), overwrite = TRUE)
  
  return()
  
}


select_DE_genes <- function(OUT_DIR,
                            label_1,
                            label_2,
                            avg_log2FC_thres,
                            p_val_adj_thres,
                            pct_thres)
{
  
  # Select DE genes
  DE_genes_list = list()
  for (i in excel_sheets(paste0(OUT_DIR,"/DE_analysis_cell_types.xlsx"))){
    temp = read.xlsx(paste0(OUT_DIR,"/DE_analysis_cell_types.xlsx"), sheet = i, rowNames = T)
    temp = temp[which(abs(temp$avg_log2FC) > avg_log2FC_thres & temp$p_val_adj < p_val_adj_thres & (temp[,paste0("pct.",label_1)] > pct_thres | temp[,paste0("pct.",label_1)] > pct_thres)),]
    temp = temp[order(temp$avg_log2FC, decreasing = T),]
    DE_genes_list[[i]] = temp
  }
  
  wb <- createWorkbook()
  for (i in excel_sheets(paste0(OUT_DIR,"/DE_analysis_cell_types.xlsx"))){
    addWorksheet(wb, i)
    writeData(wb, i, DE_genes_list[[i]], startRow = 1, startCol = 1, rowNames = T, colNames = T)
  }
  saveWorkbook(wb, file = paste0(OUT_DIR,"/DE_genes_cell_types.xlsx"), overwrite = TRUE)
  
  # DE genes stats
  df = data.frame(cell_type = names(DE_genes_list),
                  n_DE_genes = unlist(lapply(DE_genes_list, function(x){nrow(x)})),
                  n_up_genes = unlist(lapply(DE_genes_list, function(x){nrow(x[which(x$avg_log2FC>0),])})),
                  n_down_genes = unlist(lapply(DE_genes_list, function(x){nrow(x[which(x$avg_log2FC<0),])})))
  write.xlsx(df, paste0(OUT_DIR, "/DE_genes_cell_types_stats.xlsx"), rowNames = F)
  
  return()
}
  

perform_GO_enrich_analysis <- function(OUT_DIR,
                                       label_1,
                                       label_2)
{
  
  # Select DE genes
  DE_genes_list = list()
  for (i in excel_sheets(paste0(OUT_DIR,"/DE_analysis_cell_types.xlsx"))){
    temp = read.xlsx(paste0(OUT_DIR,"/DE_analysis_cell_types.xlsx"), sheet = i, rowNames = T)
    temp = temp[which(abs(temp$avg_log2FC) > avg_log2FC_thres & temp$p_val_adj < p_val_adj_thres & (temp[,paste0("pct.",label_1)] > pct_thres | temp[,paste0("pct.",label_1)] > pct_thres)),]
    temp = temp[order(temp$avg_log2FC, decreasing = T),]
    DE_genes_list[[i]] = temp
  }
  
  
  # Using gene id
  print("Perform pathway enrichment analysis using gene id")
  
  DE_go_list = list()
  for (i in excel_sheets(paste0(OUT_DIR,"/DE_analysis_cell_types.xlsx"))){
    temp_genes = as.character(sapply(gencode.v41.annotation.gene[which(gencode.v41.annotation.gene$gene_name %in% rownames(DE_genes_list[[i]])), "gene_id"], function(x){strsplit(x,"\\.")[[1]][1]}))
    ego <- enrichGO(gene          = temp_genes,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "all",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
    DE_go_list[[i]] = data.frame(ego)
  }
  
  wb <- createWorkbook()
  for (i in names(DE_go_list)){
    if (nrow(DE_go_list[[i]] > 0)){
      addWorksheet(wb, i)
      DE_go_list[[i]] = DE_go_list[[i]][order(DE_go_list[[i]]$Count, decreasing = T),]
      writeData(wb, i, DE_go_list[[i]], startRow = 1, startCol = 1, rowNames = F, colNames = T)
    }
  }
  saveWorkbook(wb, file = paste0(OUT_DIR,"/GO_enrich_gene_id.xlsx"), overwrite = TRUE)
  
  # Using symbols
  print("Perform pathway enrichment analysis using gene symbols")
  
  DE_go_list = list()
  for (i in excel_sheets(paste0(OUT_DIR,"/DE_analysis_cell_types.xlsx"))){
    ego <- enrichGO(gene          = rownames(DE_genes_list[[i]]),
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'SYMBOL',
                    ont           = "all",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
    DE_go_list[[i]] = data.frame(ego)
  }
  
  wb <- createWorkbook()
  for (i in names(DE_go_list)){
    if (nrow(DE_go_list[[i]] > 0)){
      addWorksheet(wb, i)
      DE_go_list[[i]] = DE_go_list[[i]][order(DE_go_list[[i]]$Count, decreasing = T),]
      writeData(wb, i, DE_go_list[[i]], startRow = 1, startCol = 1, rowNames = F, colNames = T)
    }
  }
  saveWorkbook(wb, file = paste0(OUT_DIR,"/GO_enrich_gene_name.xlsx"), overwrite = TRUE)
  
  return()
}




#### Perform DE analysis within the same sample between tiles ####

# DE analysis is tiles_1 vs tiles_2 (i.e. log2FC > 0 if 1 > 2)
# label_1 and label_2 are custom labels, for example conditions (normal and AD)
perform_DE_analysis_between_tiles <- function(input_data,
                                              OUT_DIR,
                                              tiles_1,
                                              tiles_2,
                                              my_genes=NA, # specific subset of genes to be used for DE analysis
                                              pct_thres=0.05,
                                              by_feature=NA, # for example "cell_type" or "region" to perform the DE analysis for each feature separately
                                              label_1,
                                              label_2
)
  
{
  
  dir.create(OUT_DIR, recursive = T, showWarnings = F)
  
  print("Read input data")
  input_data_1 = input_data[which(input_data$tile_id %in% tiles_1),]
  print(paste0("Number of spots tile group 1: ", length(unique(input_data_1$spot_id))))
  
  input_data_2 = input_data[which(input_data$tile_id %in% tiles_2),]
  print(paste0("Number of spots tile group 2: ", length(unique(input_data_2$spot_id))))
  
  # select input genes only
  if (sum(!is.na(my_genes)) > 0){
    input_data_1 = input_data_1[which(input_data_1$gene_name %in% my_genes),]
    input_data_2 = input_data_2[which(input_data_2$gene_name %in% my_genes),]
  }
  
  if (!is.na(pct_thres)){
    temp_genes_1 = names(which(table(as.factor(input_data_1$gene_name)) / length(unique(input_data_1$spot_id)) > pct_thres))
    temp_genes_2 = names(which(table(as.factor(input_data_2$gene_name)) / length(unique(input_data_2$spot_id)) > pct_thres))
    temp_genes = union(temp_genes_1, temp_genes_2)
    
    input_data_1 = input_data_1[which(input_data_1$gene_name %in% temp_genes),]
    input_data_2 = input_data_2[which(input_data_2$gene_name %in% temp_genes),]
  }
  
  print("Make spot by gene matrices")
  input_data_1_matrix = make_spot_gene_matrix(input_data_1, n_HVG=0)
  input_data_2_matrix = make_spot_gene_matrix(input_data_2, n_HVG=0)
  
  # Transform the two matrices to have the same genes
  union_genes = union(colnames(input_data_1_matrix), colnames(input_data_2_matrix))
  print(paste0("Number of genes under DE analysis: ", length(union_genes)))
  
  temp_mat = matrix(data = 0, nrow = nrow(input_data_1_matrix), ncol = length(setdiff(union_genes, colnames(input_data_1_matrix))))
  rownames(temp_mat) = rownames(input_data_1_matrix)
  colnames(temp_mat) = setdiff(union_genes, colnames(input_data_1_matrix))
  input_data_1_matrix = cbind(input_data_1_matrix, temp_mat)
  
  temp_mat = matrix(data = 0, nrow = nrow(input_data_2_matrix), ncol = length(setdiff(union_genes, colnames(input_data_2_matrix))))
  rownames(temp_mat) = rownames(input_data_2_matrix)
  colnames(temp_mat) = setdiff(union_genes, colnames(input_data_2_matrix))
  input_data_2_matrix = cbind(input_data_2_matrix, temp_mat)
  
  # Normalize data
  input_data_1_matrix_norm = log(10000 * input_data_1_matrix / rowSums(input_data_1_matrix) + 1)
  input_data_2_matrix_norm = log(10000 * input_data_2_matrix / rowSums(input_data_2_matrix) + 1)
  
  print("Perform DE analysis")
  out_list = list()
  
  if (!is.na(by_feature)){
    for (i in sort(unique(input_data_1[,by_feature]))){
      
      print(i)
      
      spot_1 = input_data_1[which(input_data_1[,by_feature] == i), "spot_id"]
      spot_2 = input_data_2[which(input_data_2[,by_feature] == i), "spot_id"]
      
      expr_1 = log(colMeans(expm1(input_data_1_matrix_norm[spot_1,]))+1)
      expr_2 = log(colMeans(expm1(input_data_2_matrix_norm[spot_2,]))+1)
      
      pct_1 = colSums(input_data_1_matrix_norm[spot_1,]>0)/length(spot_1)
      pct_2 = colSums(input_data_2_matrix_norm[spot_2,]>0)/length(spot_2)
      
      out = data.frame(expr_1, expr_2, pct_1, pct_2)
      out = out[order(out$expr_1, decreasing = T),]
      colnames(out) = paste(c("expr","expr","pct","pct"), rep(c(label_1,label_2),2), sep = ".")
      
      # DE analysis
      out$avg_log2FC = out[,paste0("expr.", label_1)] - out[,paste0("expr.", label_2)]
      out$p_val = unlist(pbmcapply::pbmclapply(rownames(out), function(x){wilcox.test(input_data_1_matrix_norm[,x], input_data_2_matrix_norm[,x])$p.value}, mc.cores = 32))
      out$p_val_adj = p.adjust(out$p_val, method = "bonferroni")
      out_list[[i]] = out
    }
    
    wb <- createWorkbook()
    for (i in sort(unique(input_data_1[,by_feature]))){
      addWorksheet(wb, i)
      writeData(wb, i, out_list[[i]], startRow = 1, startCol = 1, rowNames = T, colNames = T)
    }
    saveWorkbook(wb, file = paste0(OUT_DIR,"/DE_analysis_", by_feature, ".xlsx"), overwrite = TRUE)
  }
  
  else {
    spot_1 = input_data_1[, "spot_id"]
    spot_2 = input_data_2[, "spot_id"]
    
    expr_1 = log(colMeans(expm1(input_data_1_matrix_norm[spot_1,]))+1)
    expr_2 = log(colMeans(expm1(input_data_2_matrix_norm[spot_2,]))+1)
    
    pct_1 = colSums(input_data_1_matrix_norm[spot_1,]>0)/length(spot_1)
    pct_2 = colSums(input_data_2_matrix_norm[spot_2,]>0)/length(spot_2)
    
    out = data.frame(expr_1, expr_2, pct_1, pct_2)
    out = out[order(out$expr_1, decreasing = T),]
    colnames(out) = paste(c("expr","expr","pct","pct"), rep(c(label_1,label_2),2), sep = ".")
    
    # DE analysis
    out$avg_log2FC = out[,paste0("expr.", label_1)] - out[,paste0("expr.", label_2)]
    out$p_val = unlist(pbmcapply::pbmclapply(rownames(out), function(x){wilcox.test(input_data_1_matrix_norm[,x], input_data_2_matrix_norm[,x])$p.value}, mc.cores = 32))
    out$p_val_adj = p.adjust(out$p_val, method = "bonferroni")
    
    write.xlsx(out, paste0(OUT_DIR,"/DE_analysis.xlsx"), overwrite = TRUE, rowNames = T)
  }
  
  
  return()
  
}


select_DE_genes_between_tiles <- function(OUT_DIR,
                                          out_filename,
                                          input_spreadsheet,
                                          label_1,
                                          label_2,
                                          avg_log2FC_thres,
                                          p_val_thres,
                                          p_val_adj_thres,
                                          pct_thres)
{
  
  # Select DE genes
  DE_genes_list = list()
  for (i in excel_sheets(input_spreadsheet)){
    temp = read.xlsx(input_spreadsheet, sheet = i, rowNames = T)
    temp = temp[which(abs(temp$avg_log2FC) > avg_log2FC_thres  & temp$p_val < p_val_thres & temp$p_val_adj < p_val_adj_thres & (temp[,paste0("pct.",label_1)] > pct_thres | temp[,paste0("pct.",label_1)] > pct_thres)),]
    temp = temp[order(temp$avg_log2FC, decreasing = T),]
    DE_genes_list[[i]] = temp
  }
  
  wb <- createWorkbook()
  for (i in excel_sheets(input_spreadsheet)){
    addWorksheet(wb, i)
    writeData(wb, i, DE_genes_list[[i]], startRow = 1, startCol = 1, rowNames = T, colNames = T)
  }
  saveWorkbook(wb, file = paste0(OUT_DIR,"/", out_filename, ".xlsx"), overwrite = TRUE)
  
  # DE genes stats
  df = data.frame(cell_type = names(DE_genes_list),
                  n_DE_genes = unlist(lapply(DE_genes_list, function(x){nrow(x)})),
                  n_up_genes = unlist(lapply(DE_genes_list, function(x){nrow(x[which(x$avg_log2FC>0),])})),
                  n_down_genes = unlist(lapply(DE_genes_list, function(x){nrow(x[which(x$avg_log2FC<0),])})))
  write.xlsx(df,paste0(OUT_DIR,"/", out_filename, "_stats.xlsx"), rowNames = F)
  
  return(DE_genes_list)
}



#### Giotto analysis for spatial proximity enrichment of cell types and ligand-receptor analysis ####
giotto_analysis <- function(input_data.filter.rds,
                            tiles,
                            pct_thres=NA,
                            OUT_DIR){
  
  print("Read input data")
  input_data.filter = readRDS(input_data.filter.rds)
  input_data.filter = input_data.filter[which(input_data.filter$tile_id %in% tiles),]
  
  if (!is.na(pct_thres)){
    temp_genes = names(which(table(as.factor(input_data.filter$gene_name)) / length(unique(input_data.filter$spot_id)) > pct_thres))
    input_data.filter = input_data.filter[which(input_data.filter$gene_name %in% temp_genes),]
    
    my_out_dir = paste0(OUT_DIR, "/giotto_", pct_thres*100)
  
  } else {
    my_out_dir = paste0(OUT_DIR, "/giotto")
  }
  
  dir.create(my_out_dir, recursive = T, showWarnings = F)
  
  print("Create gene x spot matrix and coordinate matrix as input data for giotto")
  input_data.filter_matrix = t(make_spot_gene_matrix(input_data.filter, n_HVG=0))
  write.table(input_data.filter_matrix, paste0(my_out_dir,"/expression_matrix.txt"), sep = "\t", row.names = T, col.names = T, quote = F)
  
  temp = input_data.filter[!duplicated(input_data.filter$spot_id), c("spot_id", "row_i", "col_i")]
  rownames(temp) = temp$spot_id
  temp = temp[colnames(input_data.filter_matrix),c("col_i","row_i")]
  write.table(temp, paste0(my_out_dir,"/spot_coord.txt"), sep = "\t", row.names = F, col.names = F, quote = F)
  
  my_giotto_object = createGiottoObject(raw_exprs = paste0(my_out_dir,"/expression_matrix.txt"),
                                        spatial_locs = paste0(my_out_dir,"/spot_coord.txt"))
  
  print("Processing Giotto object")
  my_giotto_object <- filterGiotto(gobject = my_giotto_object, 
                                   expression_threshold = 0.5, 
                                   gene_det_in_min_cells = 5, 
                                   min_det_genes_per_cell = 0)
  my_giotto_object <- normalizeGiotto(gobject = my_giotto_object)
  
  cell_types = sort(unique(input_data.filter$cell_type))
  names(cell_types) = cell_types
  
  temp = input_data.filter[!duplicated(input_data.filter$spot_id),]
  rownames(temp) = temp$spot_id
  my_giotto_object@cell_metadata$cell_type = temp[my_giotto_object@cell_metadata$cell_ID, "cell_type"]
  
  my_giotto_object = annotateGiotto(gobject = my_giotto_object, 
                                    annotation_vector = cell_types, 
                                    cluster_column = 'cell_type', 
                                    name = 'cell_types')
  
  print("Create network (required for binSpect methods)")
  my_giotto_object = createSpatialNetwork(gobject = my_giotto_object, minimum_k = 2, name = "Delaunay_network")
  
  # print("Identify genes with a spatial coherent expression profile")
  # km_spatialgenes = binSpect(my_giotto_object, bin_method = 'kmeans')
  
  set.seed(seed = 2841)
  cell_proximities = cellProximityEnrichment(gobject = my_giotto_object,
                                             cluster_column = 'cell_types',
                                             spatial_network_name = 'Delaunay_network',
                                             adjust_method = 'fdr',
                                             number_of_simulations = 1000)
  saveRDS(cell_proximities, paste0(my_out_dir, "/cell_proximities.rds"))
  
  LR_data = fread("/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/hg38_annotation/CellTalkDB_human_lr_pair.txt", nThread = 16)[,c(2,3)]
  
  LR_data[, ligand_det := ifelse(ligand_gene_symbol %in% my_giotto_object@gene_ID, T, F)]
  LR_data[, receptor_det := ifelse(receptor_gene_symbol %in% my_giotto_object@gene_ID, T, F)]
  LR_data_det = LR_data[ligand_det == T & receptor_det == T]
  
  select_ligands = LR_data_det$ligand_gene_symbol
  select_receptors = LR_data_det$receptor_gene_symbol
  
  if (length(select_ligands) > 0 & length(select_receptors) > 0){
    print("Ligand-receptor analysis")
    
    print("Get statistical significance of gene pair expression changes based on expression")
    expr_only_scores = exprCellCellcom(gobject = my_giotto_object,
                                       cluster_column = 'cell_types',
                                       random_iter = 500,
                                       gene_set_1 = select_ligands,
                                       gene_set_2 = select_receptors)
    saveRDS(expr_only_scores, paste0(my_out_dir, "/expr_only_scores_LR.rds"))
    
    print("Get statistical significance of gene pair expression changes upon cell-cell interaction")
    spatial_all_scores = spatCellCellcom(my_giotto_object,
                                         spatial_network_name = 'Delaunay_network',
                                         cluster_column = 'cell_types',
                                         random_iter = 500,
                                         gene_set_1 = select_ligands,
                                         gene_set_2 = select_receptors,
                                         adjust_method = 'fdr',
                                         do_parallel = T,
                                         cores = 32,
                                         verbose = 'none')
    saveRDS(spatial_all_scores, paste0(my_out_dir, "/spatial_all_scores_LR.rds"))
  }

  saveRDS(my_giotto_object, paste0(my_out_dir, "/my_giotto_object.rds"))
  
  ## Combine spatial and expression based cell-cell communication data.tables
  # comb_comm = combCCcom(spatialCC = spatial_all_scores,
  #                       exprCC = expr_only_scores)
  # saveRDS(comb_comm, paste0(my_out_dir, "/comb_comm.rds"))
  
  return()
}


##### Giotto plots
giotto_plot <- function(input_dir){
  
  ### Load data
  my_giotto_object = readRDS(paste0(input_dir, "/my_giotto_object.rds"))
  cell_proximities = readRDS(paste0(input_dir, "/cell_proximities.rds"))
  
  ### Spatial proximity enrichment of cell types
  # heatmap
  p<-cellProximityHeatmap(gobject = my_giotto_object, 
                          CPscore = cell_proximities, 
                          order_cell_types = F, scale = F,
                          color_breaks = c(min(cell_proximities$enrichm_res$enrichm), 0, max(cell_proximities$enrichm_res$enrichm)),
                          color_names = c('blue', 'white', 'red'))
  png(paste0(input_dir,"/cellProximityHeatmap.png"), width = 5, height = 5, res = 300, units = "in")
  print(p)
  dev.off()
  
  # scaled heatmap
  p<-cellProximityHeatmap(gobject = my_giotto_object, 
                          CPscore = cell_proximities, 
                          order_cell_types = F, scale = T,
                          color_breaks = c(-1.5, 0, 1.5), 
                          color_names = c('blue', 'white', 'red'))
  png(paste0(input_dir,"/cellProximityHeatmap_scaled.png"), width = 5, height = 5, res = 300, units = "in")
  print(p)
  dev.off()
  
  # network
  p<-cellProximityNetwork(gobject = my_giotto_object, 
                          CPscore = cell_proximities, 
                          remove_self_edges = T, only_show_enrichment_edges = F)
  png(paste0(input_dir,"/cellProximityNetwork.png"), width = 5, height = 5, res = 300, units = "in")
  print(p)
  dev.off()
  
  # network with self-edges
  p<-cellProximityNetwork(gobject = my_giotto_object, 
                          CPscore = cell_proximities,
                          remove_self_edges = F, self_loop_strength = 0.3,
                          only_show_enrichment_edges = F,
                          rescale_edge_weights = T,
                          node_size = 8,
                          edge_weight_range_depletion = c(1,2),
                          edge_weight_range_enrichment = c(2,5))
  png(paste0(input_dir,"/cellProximityNetwork_SelfEdges.png"), width = 5, height = 5, res = 300, units = "in")
  print(p)
  dev.off()
  
  
  ### Ligand-receptor analysis
  if (file.exists(paste0(input_dir, "/spatial_all_scores_LR.rds"))){
    # expr_only_scores = readRDS(paste0(input_dir, "/expr_only_scores_LR.rds"))
    spatial_all_scores = readRDS(paste0(input_dir, "/spatial_all_scores_LR.rds"))
    # comb_comm = readRDS(paste0(input_dir, "/comb_comm.rds"))
    
    LR_data = fread("/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/hg38_annotation/CellTalkDB_human_lr_pair.txt", nThread = 16)[,c(2,3)]
    
    LR_data[, ligand_det := ifelse(ligand_gene_symbol %in% my_giotto_object@gene_ID, T, F)]
    LR_data[, receptor_det := ifelse(receptor_gene_symbol %in% my_giotto_object@gene_ID, T, F)]
    LR_data_det = LR_data[ligand_det == T & receptor_det == T]
    
    select_ligands = LR_data_det$ligand_gene_symbol
    select_receptors = LR_data_det$receptor_gene_symbol
    
    # select top LR
    selected_spat = spatial_all_scores[p.adj <= 0.5 & abs(log2fc) > 0.1 & lig_nr >= 2 & rec_nr >= 2]
    
    if (nrow(selected_spat) > 0){
      data.table::setorder(selected_spat, -PI)
      write.xlsx(selected_spat, paste0(input_dir, "/ligand_receptor_table.xlsx"))
      
      if (sum(selected_spat$lig_expr & selected_spat$rec_expr) > 0){ # at least a LR pair with non-zero expression
        
        if (length(unique(selected_spat[order(-abs(PI))]$LR_comb)) > 40){
          top_LR_ints = unique(selected_spat[order(-abs(PI))]$LR_comb)[1:40]
        } else {
          top_LR_ints = unique(selected_spat[order(-abs(PI))]$LR_comb)
        }
        
        top_LR_cell_ints = unique(selected_spat[which(selected_spat$LR_comb %in% top_LR_ints),"LR_cell_comb"])$LR_cell_comb
        # top_LR_cell_ints = unique(selected_spat[order(-abs(PI))]$LR_cell_comb)[1:33]
        
        # plotCCcomHeatmap(gobject = my_giotto_object,
        #                  comScores = spatial_all_scores,
        #                  selected_LR = top_LR_ints,
        #                  selected_cell_LR = top_LR_cell_ints,
        #                  show = 'LR_expr')
        
        p<-plotCCcomDotplot(gobject = my_giotto_object,
                            comScores = spatial_all_scores,
                            selected_LR = top_LR_ints,
                            selected_cell_LR = top_LR_cell_ints,
                            cluster_on = 'PI')
        
        png(paste0(input_dir,"/plotCCcomDotplot.png"), width = 5, height = 7, res = 300, units = "in")
        print(p)
        dev.off()
        
        # # top differential activity levels for ligand receptor pairs
        # plotRankSpatvsExpr(gobject = my_giotto_object,
        #                    comb_comm,
        #                    expr_rnk_column = 'exprPI_rnk',
        #                    spat_rnk_column = 'spatPI_rnk',
        #                    midpoint = 10)
        # 
        # ## * recovery ####
        # ## predict maximum differential activity
        # plotRecovery(gobject = my_giotto_object,
        #              comb_comm,
        #              expr_rnk_column = 'exprPI_rnk',
        #              spat_rnk_column = 'spatPI_rnk',
        #              ground_truth = 'spatial')
      }
    }
  }

  return()
  
}


#####Cell type proximity differential analysis
cell_type_proximity_differential_analysis <- function(main_dir,
                                                      sample_1,
                                                      sample_2,
                                                      ROI_label_1,
                                                      ROI_label_2,
                                                      label_1,
                                                      label_2,
                                                      pct_thres=NA)
{
  
  if(!is.na(pct_thres)){
    out_dir = paste0(main_dir, "/DE_analysis/", sample_1, "_vs_", sample_2, "/pct_", pct_thres*100, "/")
    cell_proximities_1 = readRDS(paste0(main_dir, "/", sample_1, "/", ROI_label_1, "/cell_types/giotto_", pct_thres*100,"/cell_proximities.rds"))
    cell_proximities_2 = readRDS(paste0(main_dir, "/", sample_2, "/", ROI_label_2, "/cell_types/giotto_", pct_thres*100,"/cell_proximities.rds"))
  } else {
    out_dir = paste0(main_dir, "/DE_analysis/", sample_1, "_vs_", sample_2, "/")
    cell_proximities_1 = readRDS(paste0(main_dir, "/", sample_1, "/", ROI_label_1, "/cell_types/giotto/cell_proximities.rds"))
    cell_proximities_2 = readRDS(paste0(main_dir, "/", sample_2, "/", ROI_label_2, "/cell_types/giotto/cell_proximities.rds"))
  }
  
  # p_higher_orig: p-value associated with higher observed than expected
  # p_lower_orig: p-value associated with lower observed than expected
  cell_proximities_1 = data.frame(cell_proximities_1$enrichm_res)
  cell_proximities_2 = data.frame(cell_proximities_2$enrichm_res)
  
  # Sort by cell-cell interaction
  temp1 = cell_proximities_1[order(as.character(cell_proximities_1$unified_int)), c("unified_int","enrichm","PI_value")]
  temp2 = cell_proximities_2[order(as.character(cell_proximities_2$unified_int)), c("unified_int","enrichm","PI_value")]
  
  df_out = merge(temp1, temp2, by = "unified_int")[,c(1,2,4,3,5)]
  colnames(df_out)[2:5] = paste0(c(rep("enrichm.", 2), rep("PI_value.",2)), c(label_1,label_2))
  
  df_out$log2FC.enrichm = df_out[,paste0("enrichm.",label_1)] - df_out[,paste0("enrichm.",label_2)]
  df_out$log2FC.PI_value = df_out[,paste0("PI_value.",label_1)] - df_out[,paste0("PI_value.",label_2)]
  df_out = df_out[order(df_out$log2FC.PI_value, decreasing = T),]
  write.xlsx(df_out, paste0(out_dir, "/cell_type_spatial_proximity.xlsx"), rowNames=F)
  
  print(paste0("p_val PI_value: ", t.test(temp1$PI_value, temp2$PI_value, paired = T)$p.value))
  print(paste0("p_val enrichm: ", t.test(temp1$enrichm, temp2$enrichm, paired = T)$p.value))
  
  ### Plot heatmap with enrichment value
  for (j in c("log2FC.enrichm", "log2FC.PI_value")){
    df_heatmap = df_out[,c("unified_int",j)]
    df_heatmap$row = sapply(as.character(df_heatmap$unified_int), function(x){strsplit(x, "--")[[1]][1]})
    df_heatmap$col = sapply(as.character(df_heatmap$unified_int), function(x){strsplit(x, "--")[[1]][2]})
    
    m1 <- matrix(0, length(unique(df_heatmap$row)), length(unique(df_heatmap$col)))
    rownames(m1) = sort(unique(df_heatmap$row))
    colnames(m1) = sort(unique(df_heatmap$row))
    m1[as.matrix(df_heatmap[,c("row","col")])] <- df_heatmap[,j]
    
    temp_matrix = m1 + t(m1) - diag(length(unique(df_heatmap$row)))*m1
    
    col_fun = colorRamp2(breaks = c(quantile(temp_matrix[!is.na(temp_matrix)], 0.02), 
                                    0,
                                    quantile(temp_matrix[!is.na(temp_matrix)], 0.98)), 
                         colors=c("blue","white","red"))
    
    p1<-Heatmap(temp_matrix, cluster_rows = F, cluster_columns = F, na_col = "black", col = col_fun, column_title = "",
                heatmap_legend_param = list(title = ""))
    png(paste0(out_dir, "heatmap_cell_type_proximities_", j, ".png"), width = 4.5, height = 4, res = 200, units = "in")
    print(p1)
    dev.off()
  }
  
  return()
}










