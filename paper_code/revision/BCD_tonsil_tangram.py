import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd
from anndata import read_h5ad
import tangram as tg
import os
import matplotlib.pyplot as plt

def run_tangram(adata_sc, adata_sp, sample_label, cluster_label='annotation_figure_1', markers=None, device='cpu'):
    print(f"Preprocessing AnnDatas with Tangram for sample {sample_label}...")
    tg.pp_adatas(adata_sc, adata_sp, genes=markers)

    print(f"Starting Tangram mapping sample {sample_label}...")
    ad_map = tg.map_cells_to_space(
        adata_sc,
        adata_sp,
        mode='clusters',
        cluster_label=cluster_label,
        density_prior='rna_count_based',
        device=device
    )
    output_path = f'{sample_label}_cell_map_result.h5ad'
    print(f"Saving mapped AnnData to {output_path} ...")
    ad_map.write_h5ad(output_path)
    print("Saved successfully.")
    
    print(f"Projecting cell annotations sample {sample_label}...", flush=True)
    tg.project_cell_annotations(ad_map, adata_sp, annotation=cluster_label)
    annotation_list = list(pd.unique(adata_sc.obs[cluster_label]))
    
    print(f"Plotting cell annotations sample {sample_label}...", flush=True)
    tg.plot_cell_annotation_sc(adata_sp, annotation_list, perc=0.02)
    plt.savefig(f'st_cell_mapping_{sample_label}.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved st_cell_mapping_{sample_label}.png", flush=True)
    
    print(f"Plotting training scores sample {sample_label}...", flush=True)
    tg.plot_training_scores(ad_map, bins=20, alpha=.5)
    plt.savefig(f'training_metrics_{sample_label}.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved training_metrics_{sample_label}.png", flush=True)
    
    return ad_map

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
adata_sp_b = sc.read_visium('/project2/tays/Junjie/10x/spatialProx/20250227/20250227_B1_Human_tonsil_Proxseq/outs')
adata_sp_c = sc.read_visium('/project2/tays/Junjie/10x/spatialProx/20250227/20250227_C1_Human_tonsil_Proxseq/outs')
adata_sp_d = sc.read_visium('/project2/tays/Junjie/10x/spatialProx/20250227/20250227_D1_Human_tonsil_Proxseq/outs')
print("Loaded spatial AnnData:")

# Rank marker genes per cell type
print("Ranking marker genes for each cell type...")
sc.tl.rank_genes_groups(adata_sc, groupby="annotation_figure_1", use_raw=False)

# Extract top marker genes (up to 100)
markers_df = pd.DataFrame(adata_sc.uns["rank_genes_groups"]["names"]).iloc[0:100, :]
markers = list(np.unique(markers_df.melt().value.values))
print(f"Number of unique marker genes selected: {len(markers)}")

ad_map_b = run_tangram(adata_sc, adata_sp_b, 'B1', markers=markers, device='cpu')
ad_map_c = run_tangram(adata_sc, adata_sp_c, 'C1', markers=markers, device='cpu')
ad_map_d = run_tangram(adata_sc, adata_sp_d, 'D1', markers=markers, device='cpu')

print("Tangram deconvolution script finished.")