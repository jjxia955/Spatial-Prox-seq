
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import gseapy as gp
from adjustText import adjust_text

def analyze_pla_marker_with_gsea_and_volcano_new(
    adata_pla,
    adata_rna,
    pla_marker,
    cluster_labels=None,
    cluster_key='mrna_annotation',
    gene_set_path='c7.all.v2025.1.Hs.symbols.gmt',
    low_pct=25,
    high_pct=75,
    min_cells=100,
    gsea_output_dir='gsea_output',
    top_n_terms=5,
    fc_thresh=1,
    p_thresh=0.05,
    fc_limit=10,
    volcano_title=None,
    output_prefix=None,
    use_raw_pla=False,
    batch_filter=None,
    verbose=True
):
    EPS = 1e-300  # for log stability

    if batch_filter is not None:
        if isinstance(batch_filter, str):
            batch_filter = [batch_filter]
        adata_pla = adata_pla[adata_pla.obs["batch"].isin(batch_filter)].copy()
        adata_rna = adata_rna[adata_rna.obs["batch"].isin(batch_filter)].copy()
        if verbose:
            print(f"[INFO] Filtering by batch: {batch_filter}")
            print(f"adata_pla shape after filtering: {adata_pla.shape}")
            print(f"adata_rna shape after filtering: {adata_rna.shape}")

    if output_prefix:
        gsea_output_dir = os.path.join(gsea_output_dir, output_prefix)
    os.makedirs(gsea_output_dir, exist_ok=True)

    if pla_marker not in adata_pla.var_names:
        raise ValueError(f"PLA marker {pla_marker} not found in adata_pla.var_names.")

    if isinstance(cluster_labels, str):
        cluster_labels = [cluster_labels]
    is_target_cluster = adata_pla.obs[cluster_key].isin(cluster_labels)
    if is_target_cluster.sum() == 0:
        raise ValueError(f"No cells found for cluster labels: {cluster_labels}")
    adata_cluster = adata_pla[is_target_cluster]

    if use_raw_pla and adata_pla.raw is not None:
        marker_values = adata_pla.raw[adata_cluster.obs_names, pla_marker].X
    else:
        marker_values = adata_cluster[:, pla_marker].X

    if not isinstance(marker_values, np.ndarray):
        marker_values = marker_values.toarray()
    marker_values = marker_values.ravel()

    n_total = len(marker_values)
    n_high = int(n_total * (high_pct / 100))
    n_low = int(n_total * (low_pct / 100))

    if n_high < min_cells or n_low < min_cells:
        raise ValueError(f"Too few cells selected for high ({n_high}) or low ({n_low}) group.")

    sorted_idx = np.argsort(marker_values)
    low_idx = sorted_idx[:n_low]
    high_idx = sorted_idx[-n_high:]

    obs_names = adata_cluster.obs_names.to_numpy()
    low_cells = obs_names[low_idx]
    high_cells = obs_names[high_idx]

    if verbose:
        print(f"[INFO] PLA marker: {pla_marker} | Clusters: {cluster_labels}")
        print(f"High group: {len(high_cells)} cells | Low group: {len(low_cells)} cells")

        high_batches = adata_cluster.obs.loc[high_cells, "batch"].value_counts()
        low_batches = adata_cluster.obs.loc[low_cells, "batch"].value_counts()

        print("[INFO] Batch distribution for high group:")
        print(high_batches)
        print("[INFO] Batch distribution for low group:")
        print(low_batches)

    batch_distribution = {
        "high": high_batches.to_dict(),
        "low": low_batches.to_dict()
    }

    batch_df = pd.DataFrame({'High group': high_batches, 'Low group': low_batches}).fillna(0).astype(int).T
    plt.figure(figsize=(2.5, 2.5))
    batch_df.plot(kind='bar', edgecolor='black', width=0.8)
    plt.ylabel('Spot Count')
    plt.title(f'Batch Distribution: {pla_marker} in {",".join(cluster_labels)}')
    plt.xticks(rotation=0)
    plt.tight_layout()
    if output_prefix:
        plot_path = os.path.join(gsea_output_dir, f"{output_prefix}_batch_distribution.svg")
        plt.savefig(plot_path, dpi=300)
        print(f"✅ Batch distribution plot saved to {plot_path}")
        plt.show()
    else:
        plt.show()

    adata_high = adata_rna[adata_rna.obs_names.isin(high_cells)].copy()
    adata_low = adata_rna[adata_rna.obs_names.isin(low_cells)].copy()

    if adata_rna.raw is not None:
        genes = adata_high.var_names
        adata_high.X = adata_rna.raw[adata_high.obs_names, genes].X.copy()
        adata_low.X = adata_rna.raw[adata_low.obs_names, genes].X.copy()

    adata_high.obs["pla_group"] = "high"
    adata_low.obs["pla_group"] = "low"
    adata_combined = adata_high.concatenate(
        adata_low,
        batch_key="pla_group",
        batch_categories=["high", "low"],
        index_unique=None
    )

    sc.tl.rank_genes_groups(
        adata_combined,
        groupby="pla_group",
        method="wilcoxon",
        use_raw=False
    )
    deg_df = sc.get.rank_genes_groups_df(adata_combined, group='high')
    deg_df_for_gsea = deg_df[deg_df['pvals_adj'] < 0.05]

    if deg_df_for_gsea.empty:
        raise ValueError("No significant DEGs found (FDR < 0.05), cannot perform GSEA.")

    rnk = deg_df_for_gsea[['names', 'logfoldchanges']].dropna()
    rnk = rnk.set_index('names').squeeze().sort_values(ascending=False)

    gsea_result = gp.prerank(
        rnk=rnk,
        gene_sets=gene_set_path,
        outdir=gsea_output_dir,
        min_size=15,
        max_size=500,
        permutation_num=1000,
        seed=42,
        format='svg'
    )

    gsea_df = gsea_result.res2d.copy()
    gsea_df.columns = [col.lower() for col in gsea_df.columns]
    if 'fdr q-val' in gsea_df.columns:
        gsea_df.rename(columns={'fdr q-val': 'fdr'}, inplace=True)
    if "term" in gsea_df.columns:
        gsea_df = gsea_df.set_index("term")

    def plot_gsea_terms(subset_df, direction):
        top_terms = subset_df.reindex(subset_df["nes"].abs().sort_values(ascending=False).head(top_n_terms).index).index.tolist()
        for term in top_terms:
            try:
                result = gp.gseaplot(
                    rank_metric=rnk,
                    term=term,
                    figsize=(5, 5),
                    **gsea_result.results[term]
                )
                if isinstance(result, (list, tuple)) and len(result) == 2:
                    fig, ax = result
                    fig.suptitle(f"{direction.capitalize()} Enrichment: {term} (NES={gsea_df.loc[term, 'nes']:.2f})")
                    fig.tight_layout()
                    save_path = os.path.join(gsea_output_dir, f"{output_prefix}_{direction}_{term}.svg")
                    fig.savefig(save_path, bbox_inches="tight")
                    plt.close(fig)
                    print(f"✅ Saved: {save_path}")
            except KeyError:
                print(f"❌ Term not found in GSEA result: {term}")

    print("=== Top Upregulated Pathways (NES > 0) ===")
    plot_gsea_terms(gsea_df[gsea_df["nes"] > 0], "positive")

    print("\n=== Top Downregulated Pathways (NES < 0) ===")
    plot_gsea_terms(gsea_df[gsea_df["nes"] < 0], "negative")

    volcano_df = deg_df.copy()
    volcano_df['-log10(pval)'] = -np.log10(volcano_df['pvals_adj'] + EPS)
    volcano_df['category'] = 'Not significant'
    volcano_df.loc[
        (volcano_df['logfoldchanges'] > fc_thresh) & (volcano_df['pvals_adj'] < p_thresh),
        'category'
    ] = 'Up'
    volcano_df.loc[
        (volcano_df['logfoldchanges'] < -fc_thresh) & (volcano_df['pvals_adj'] < p_thresh),
        'category'
    ] = 'Down'
    extreme_mask = (
        (volcano_df['pvals_adj'] > p_thresh) &
        (np.abs(volcano_df['logfoldchanges']) > fc_limit)
    )
    plot_df = volcano_df[~extreme_mask].copy()

    plt.figure(figsize=(3.5, 3.5))
    sns.scatterplot(
        data=plot_df,
        x='logfoldchanges',
        y='-log10(pval)',
        hue='category',
        palette={'Up': 'red', 'Down': 'blue', 'Not significant': 'black'},
        edgecolor=None,
        alpha=0.8,
        s=20,
        legend=False
    )
    plt.axhline(-np.log10(p_thresh), linestyle='--', color='gray')
    plt.axvline(fc_thresh, linestyle='--', color='gray')
    plt.axvline(-fc_thresh, linestyle='--', color='gray')

    up_genes = plot_df[(plot_df['category'] == 'Up')].nsmallest(top_n_terms, 'pvals')
    down_genes = plot_df[(plot_df['category'] == 'Down')].nsmallest(top_n_terms, 'pvals')

    texts = []
    for _, row in up_genes.iterrows():
        texts.append(plt.text(row['logfoldchanges'], row['-log10(pval)'], row['names'], fontsize=8, color='red'))
    for _, row in down_genes.iterrows():
        texts.append(plt.text(row['logfoldchanges'], row['-log10(pval)'], row['names'], fontsize=8, color='blue'))

    adjust_text(texts, arrowprops=dict(arrowstyle="-", color='gray', lw=0.5))

    plt.title(volcano_title or f"{pla_marker}")
    plt.xlabel("Log2 Fold Change")
    plt.ylabel("-Log10(P Value)")
    plt.tight_layout()
    if output_prefix:
        volcano_path = os.path.join(gsea_output_dir, f"{output_prefix}_volcano.svg")
        plt.savefig(volcano_path, dpi=300)
        print(f"✅ Volcano plot saved to {volcano_path}")
        plt.show()
    else:
        plt.show()

    return {
        "deg_df": deg_df,
        "gsea_df": gsea_df,
        "gsea_result": gsea_result,
        "batch_distribution": batch_distribution,
        "adata_combined": adata_combined,
        "high_cells": high_cells.tolist(),
        "low_cells": low_cells.tolist()
    }
