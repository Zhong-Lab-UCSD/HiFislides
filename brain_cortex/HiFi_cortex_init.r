####### INPUT SAMPLES (RUN THE SECTION RELATED TO THE SAMPLE OF INTEREST)

#### HiFi_cortex_1 ####
SAMPLE_NAME="HiFi_cortex_1"

ROI_label = "ROI_1"
n_tiles_row_ROI = 1
n_tiles_col_ROI = 1

ROI_label = "ROI_2"
n_tiles_row_ROI = 2
n_tiles_col_ROI = 2

ROI_label = "ROI_3"
n_tiles_row_ROI = 3
n_tiles_col_ROI = 3

flowcell_type = "NextSeq"
surface = "1"
first_tile_position = "top-left" # this indicates the position of the first tile either "top-left" or "bottom-right"
n_tiles_row = 14
n_tiles_col = 6

#### HiFi_cortex_AD_1 ####
SAMPLE_NAME="HiFi_cortex_AD_1"

ROI_label = "ROI_3"
n_tiles_row_ROI = 4
n_tiles_col_ROI = 6

flowcell_type = "NextSeq"
surface = "1"
first_tile_position = "top-left" # this indicates the position of the first tile either "top-left" or "bottom-right"
n_tiles_row = 14
n_tiles_col = 6


#### HiFi_cortex_AD_2 ####
SAMPLE_NAME="HiFi_cortex_AD_2"

ROI_label = "ROI_2"
n_tiles_row_ROI = 2
n_tiles_col_ROI = 2

flowcell_type = "NextSeq"
surface = "1"
flowcell_half=".top."
first_tile_position = "top-left" # this indicates the position of the first tile either "top-left" or "bottom-right"
n_tiles_row = 14
n_tiles_col = 6
my_giotto_tiles = 1404
#### HiFi_cortex_3 ####
SAMPLE_NAME="HiFi_cortex_3"

ROI_label = "ROI_1"
n_tiles_row_ROI = 4
n_tiles_col_ROI = 2

flowcell_type = "NextSeq"
surface = "1"
first_tile_position = "top-left" # this indicates the position of the first tile either "top-left" or "bottom-right"
n_tiles_row = 14
n_tiles_col = 6
#### HiFi_cortex_AD_3 ####
SAMPLE_NAME="HiFi_cortex_AD_3"

ROI_label = "ROI_1"
n_tiles_row_ROI = 3
n_tiles_col_ROI = 3

flowcell_type = "NextSeq"
surface = "1"
first_tile_position = "top-left" # this indicates the position of the first tile either "top-left" or "bottom-right"
n_tiles_row = 14
n_tiles_col = 6
#### HiFi_cortex_4 ####
SAMPLE_NAME="HiFi_cortex_4"

# ROI_label = "ROI_1"
# n_tiles_row_ROI = 5
# n_tiles_col_ROI = 3

ROI_label = "ROI_2"
n_tiles_row_ROI = 2
n_tiles_col_ROI = 2

flowcell_type = "NextSeq"
surface = "1"
flowcell_half=".bottom."
first_tile_position = "top-left" # this indicates the position of the first tile either "top-left" or "bottom-right"
n_tiles_row = 14
n_tiles_col = 6
my_giotto_tiles = 1413
#### HiFi_cortex_5 ####
SAMPLE_NAME="HiFi_cortex_5"

ROI_label = "ROI_1"
n_tiles_row_ROI = 4
n_tiles_col_ROI = 2

flowcell_type = "NextSeq"
surface = "1"
first_tile_position = "top-left" # this indicates the position of the first tile either "top-left" or "bottom-right"
n_tiles_row = 11
n_tiles_col = 6
#### END OF INPUT SAMPLES SECTION ####


##### Setup directories
OUT_DIR = paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/result/cortex/", SAMPLE_NAME, "/", ROI_label)
DATA_DIR=paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/IGM/", SAMPLE_NAME)
dir.create(OUT_DIR, recursive = T, showWarnings = F)

##### Create tile matrix and plot number of spots and expressed genes over the flowcell
tile_matrix = create_tile_matrix(flowcell_type,
                                 first_tile_position,
                                 n_tiles_row,
                                 n_tiles_col,
                                 surface)

row_tile_offset = 74500
col_tile_offset = 56500

tile_matrix_expanded = transform(expand.grid(i = seq(nrow(tile_matrix)), j = seq(ncol(tile_matrix))), tile_id = c(tile_matrix))
tile_matrix_expanded$row = (tile_matrix_expanded$i - 1) * row_tile_offset
tile_matrix_expanded$col = (tile_matrix_expanded$j - 1) * col_tile_offset
write.table(tile_matrix_expanded, paste0("/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/barcodes/tile_matrix_expanded_", n_tiles_row, "x", n_tiles_col, ".txt"),
            row.names = F, col.names = T, quote = F, sep = "\t")


# To select tiles from the top half or bottom half
half = "top"

if (half == "bottom"){
  if (n_tiles_row %% 2 == 0){
    custom_tiles = as.numeric(tile_matrix[(n_tiles_row %/% 2 + 1):n_tiles_row,])
  } else {
    custom_tiles = as.numeric(tile_matrix[(n_tiles_row %/% 2 + 1):n_tiles_row,])
  }
} else if (half == "top"){
  if (n_tiles_row %% 2 == 0){
    custom_tiles = as.numeric(tile_matrix[1:(n_tiles_row %/% 2),])
  } else {
    custom_tiles = as.numeric(tile_matrix[1:(n_tiles_row %/% 2 + 1),])
  }
}

plot_flowcell_heatmaps(INTEGRATE_DIR=paste0(DATA_DIR, "/L2R1_L2R2_integrate.N1.top.k80./"),
                       SAMPLE_NAME,
                       custom_tiles=NA,
                       normalization_file="/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/data/barcodes/AAANLCHHV_1_1/tile_barcode_number.txt",
                       scaling_factor=1000000,
                       n_tiles_row,
                       n_tiles_col,
                       n_min_read=NA,
                       n_max_read=NA,
                       n_min_spot=NA,
                       n_max_spot=NA,
                       n_min_gene=NA,
                       n_max_gene=NA)


##### Reading input data and compute ROI stats
ROI_tiles = read.table(paste0(DATA_DIR, "/", ROI_label, "_tiles.txt"))$V1
INPUT_DIR = paste0(DATA_DIR,"/L2R1_L2R2_integrate", flowcell_half, "/")

input_data = read_HiFi_input_data(INPUT_DIR=INPUT_DIR,
                                  SAMPLE_NAME,
                                  ROI_label,
                                  tile_matrix)

calculate_ROI_stats(INPUT_DIR=INPUT_DIR,
                    SAMPLE_NAME,
                    OUT_DIR,
                    input_data)


