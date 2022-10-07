#!/usr/bin/env Rscript

library("optparse")
library("data.table")

##### Read input argument
option_list = list(
  make_option(c("-i", "--input_file"), type="character", default=NULL, 
              help="L2R1__L1R1_dedup.hifislida2.o", metavar="character"),
  make_option(c("-o", "--out_file"), type="character", default=NULL,
              help="Output file", metavar="character"),
  make_option(c("-M", "--max_size_ROI"), type="character", default=NULL,
              help="Maximum size of the ROI, either on the rows or the columns", metavar="character"),
  make_option(c("-m", "--min_size_ROI"), type="character", default=NULL,
              help="Maximum size of the ROI, either on the rows or the columns", metavar="character"),
  make_option(c("-p", "--p_value"), type="character", default=NULL,
              help="p-value threshold", metavar="character")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input_file)){
  print_help(opt_parser)
  stop("Input file L2R1__L1R1_dedup.hifislida2.o must be supplied!", call.=FALSE)
}

if (is.null(opt$out_file)){
  print_help(opt_parser)
  stop("Output file with full path must be supplied!", call.=FALSE)
}


#####
input_data = data.frame(fread(opt$input_file))
#input_data = data.frame(fread("/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/test_sample_new/lib2/L2R1_mapping/L2R1__L1R1_dedup_1_1.hifislida2.o"))

input_data = input_data[,-1]
colnames(input_data) = c("tile", "n_HiFi_reads", "n_barcodes")

# tile_matrix is the matrix for arranging Tile IDs as their relative position on the surface. We have 11 x 6 tiles (rows x cols)
x = 1101:1111
tile_matrix = cbind(x, x+100, x+200, x+300, x+400, x+500)

mat = matrix(data=NA, nrow=1, ncol=6)

### We estimate the ROI to be rectangular.

# max_size_ROI: the maximum size of a side of the ROI 
# max_size_ROI = 4 indicates that ROI can be at most 4 x 4 tiles
# min_size_ROI: the minimum size of a side of the ROI 
# min_size_ROI = 2 indicates that ROI shall be at least 2 x 2 tiles


max_size_ROI = as.numeric(opt$max_size_ROI) # either on the rows or columns
min_size_ROI = as.numeric(opt$min_size_ROI) # either on the rows or columns

w_matrix = cbind(combn(c(min_size_ROI:max_size_ROI),2), combn(c(min_size_ROI:max_size_ROI),2)[c(2,1),], 
                 matrix(data = c(min_size_ROI:max_size_ROI, min_size_ROI:max_size_ROI), 2, max_size_ROI-min_size_ROI+1, byrow = T))

out_list = list()

for(w in 1:ncol(w_matrix)) {
  w_row = w_matrix[1,w]
  w_col = w_matrix[2,w]
  d_row = w_row - 1
  d_col = w_col - 1
  
  for(i in 1:(11-d_row)) {
    for(j in 1:(6-d_col)) {
      
      tiles_ROI = which(input_data$tile %in% tile_matrix[i:(i+d_row),j:(j+d_col)]) # tiles in putative ROI
      tiles_noROI = which(!(input_data$tile %in% tile_matrix[i:(i+d_row),j:(j+d_col)])) # tiles not in putative ROI
      
      p_value = ks.test(input_data[tiles_ROI,"n_HiFi_reads"],input_data[tiles_noROI,"n_HiFi_reads"])$p.value
      
      out = data.frame(w_row, w_col, i, j, avg_HiFi_reads_ROI = mean(input_data[tiles_ROI,"n_HiFi_reads"]), avg_HiFi_reads_noROI = mean(input_data[tiles_noROI,"n_HiFi_reads"]), p_value)
      out_list[[paste0(w, "_", i, "_", j)]] = out
    }
  }
}

df = do.call("rbind", out_list)
df_1 = df[which(df$avg_HiFi_reads_ROI > df$avg_HiFi_reads_noROI & df$p_value < as.numeric(opt$p_value)),]
df_1$avg_log2FC = log2((df_1$avg_HiFi_reads_ROI + 1) / (df_1$avg_HiFi_reads_noROI + 1))

#df_1 = df_1[order(df_1$avg_HiFi_reads_ROI, decreasing = T),]
df_1 = df_1[order(df_1$avg_log2FC, decreasing = T),]

df_roi = df_1[1,]
roi = tile_matrix[df_roi$i:(df_roi$i+df_roi$w_row-1), df_roi$j:(df_roi$j+df_roi$w_col-1)]

#print(as.numeric(roi))
write.table(as.numeric(roi), opt$out_file, row.names = F, col.names = F, quote = F)



