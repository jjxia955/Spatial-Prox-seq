import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd
from anndata import read_h5ad
import tangram as tg
import os
import matplotlib.pyplot as plt

print("Starting Tangram deconvolution script...")

# Load single-cell reference
adata_sc = read_h5ad('/project2/tays/Junjie/10x/spatialProx/20250227/sce_downsampled_50k.h5ad')
print("Loaded single-cell AnnData:")
print(adata_sc)

# Check if counts are normalized (example for first cell)
print("Checking counts normalization for first cell:")
if hasattr(adata_sc.X, "toarray"):
    unique_counts = np.unique(adata_sc.X.toarray()[0, :])
else:
    unique_counts = np.unique(adata_sc.X[0, :])
print(unique_counts)

sc.pp.normalize_total(adata_sc)
# Log-transform normalized counts
sc.pp.log1p(adata_sc)

# Load spatial data (Visium)
adata_sp = sc.read_visium('/project2/tays/Junjie/10x/spatialProx/20250227/20250227_A1_Human_tonsil_Proxseq/outs')
print("Loaded spatial AnnData:")
print(adata_sp)

# Rank marker genes per cell type
print("Ranking marker genes for each cell type...")
sc.tl.rank_genes_groups(adata_sc, groupby="annotation_figure_1", use_raw=False)

# Extract top marker genes (up to 100)
markers_df = pd.DataFrame(adata_sc.uns["rank_genes_groups"]["names"]).iloc[0:100, :]
markers = list(np.unique(markers_df.melt().value.values))
print(f"Number of unique marker genes selected: {len(markers)}")

# Preprocess AnnDatas for Tangram
print("Preprocessing AnnDatas with Tangram...")
tg.pp_adatas(adata_sc, adata_sp, genes=markers)

# Check GPU availability before mapping
print("Checking GPU status via nvidia-smi:")
os.system("nvidia-smi")

# Run Tangram cell-to-space mapping
print("Starting Tangram mapping...")
ad_map = tg.map_cells_to_space(
    adata_sc,
    adata_sp,
    #mode="cells",
    mode='clusters',
    cluster_label='annotation_figure_1',
    density_prior='rna_count_based',
    #num_epochs=500,
    #device="cuda:0"
    device='cpu'
)
print("Tangram mapping completed.")

# Save mapped AnnData
output_path = 'A1_cell_map_result.h5ad'
print(f"Saving mapped AnnData to {output_path} ...")
ad_map.write_h5ad(output_path)
print("Saved successfully.")

print("Tangram deconvolution script finished.")

print("Projecting cell annotations...", flush=True)
tg.project_cell_annotations(ad_map, adata_sp, annotation="annotation_figure_1")
annotation_list = list(pd.unique(adata_sc.obs['annotation_figure_1']))
print(f"Annotations to plot: {annotation_list}", flush=True)

print("Plotting cell annotations...", flush=True)
tg.plot_cell_annotation_sc(adata_sp, annotation_list, perc=0.02)
plt.savefig('st_cell_mapping.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved st_cell_mapping.png", flush=True)

print("Plotting training scores...", flush=True)
tg.plot_training_scores(ad_map, bins=20, alpha=.5)
plt.savefig('training_metrics.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved training_metrics.png", flush=True)