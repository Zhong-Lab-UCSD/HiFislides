# /home/rcalandrelli/anaconda3/envs/spaceflow_env/lib/python3.7/site-packages/SpaceFlow

### 1. Import SpaceFlow and squidpy package
import squidpy as sq
import scanpy as sc
from SpaceFlow import SpaceFlow
import numpy as np
import pandas as pd
import anndata as ad
import scipy
import matplotlib
import matplotlib.pyplot as plt

### 2. Load the ST data from squidpy package
adata_example = sq.datasets.seqfish() # example from documentation

SAMPLE_DIR = "/mnt/extraids/SDSC_NFS/rcalandrelli/HiFi/result/brain_organoid/HiFi_organoid_1/fastp_filter_k19_ALL/spaceflow_3/"

df = pd.read_csv(SAMPLE_DIR + "spaceflow_matrix.txt", sep = "\t")
counts = np.float32(df.values)

adata = ad.AnnData(counts)
adata.obs_names = df.index
adata.var_names = df.columns
adata.X = scipy.sparse.csr_matrix(adata.X)

sc.pp.filter_genes(adata, min_cells=3)

# adata.X[1,1]
# print(adata.obs_names[:10])
# print(adata.var_names[:10])

### 3. Create SpaceFlow Object
# adata: the count matrix of gene expression, 2D numpy array of size (# of cells, # of genes), type anndata.AnnData, see https://anndata.readthedocs.io/en/latest/ for more info about anndata.
# spatial_locs: spatial locations of cells (or spots) match to rows of the count matrix, 1D numpy array of size (n_locations,), type numpy.ndarray, optional
adata_coord = pd.read_csv(SAMPLE_DIR + "spaceflow_coord.txt", sep = "\t").values
adata.obsm['spatial'] = adata_coord

adata_cell_type = pd.read_csv(SAMPLE_DIR + "spaceflow_cell_type.txt", sep = "\t").values.tolist()
adata.obs['cell_type'] = pd.Categorical(np.array([item for sublist in adata_cell_type for item in sublist]))

sf = SpaceFlow.SpaceFlow(adata=adata, spatial_locs=adata.obsm['spatial'])

### 4. Preprocessing the ST Data
sf.preprocessing_data(n_top_genes=1000, flavor="seurat")
# sf.preprocessing_data(n_top_genes=adata.shape[1], flavor="seurat")

### 5. Train the deep graph network model
sf.train(spatial_regularization_strength=0.1, z_dim=50, lr=1e-3, epochs=1000, max_patience=50, min_stop=100, random_seed=42, gpu=0, regularization_acceleration=True, edge_subset_sz=1000000)

### 6. Domain segmentation of the ST data
sf.segmentation(domain_label_save_filepath=SAMPLE_DIR + "domains.tsv", n_neighbors=50, resolution=1.0)

### 7. Visualization of the annotation and the identified spatial domains
sf.plot_segmentation(segmentation_figure_save_filepath=SAMPLE_DIR + "domain_segmentation.pdf", colormap="tab20", scatter_sz=1., rsz=4., csz=4., wspace=.4, hspace=.5, left=0.125, right=0.9, bottom=0.1, top=0.9)

sc.pl.spatial(adata, color="cell_type", spot_size=6000, return_fig = True)
plt.savefig(SAMPLE_DIR + "cell_types.pdf")

### 8. Idenfify the spatiotemporal patterns of the ST data through pseudo-Spatiotemporal Map (pSM)
sf.pseudo_Spatiotemporal_Map(pSM_values_save_filepath=SAMPLE_DIR + "pSM_values.tsv", n_neighbors=20, resolution=1.0)

sf.plot_pSM(pSM_figure_save_filepath=SAMPLE_DIR + "pseudo-Spatiotemporal-Map.pdf", colormap="roma", scatter_sz=1., rsz=4., csz=4., wspace=.4, hspace=.5, left=0.125, right=0.9, bottom=0.1, top=0.9)


