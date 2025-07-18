
# SpatialProxSeq: Spatial Proximity Sequencing Analysis Toolkit

Author: Junjie Xia  
Email: jjxia@uchicago.edu

## Overview

`SpatialproxseqClasses.py` provides a computational framework for analyzing *Spatial Prox-Seq* datasets. This tool is designed for proximity ligation assay (PLA) data across spatially barcoded spots, allowing for prediction and quantification of protein interactions, and network-based interpretation of protein interactions in situ.

Key features:
- Calculates **protein abundance** and pairwise **protein proximity counts**
- Calculates **Fisher's exact test**, **Pearson residuals**, and **Fractional Overlap** metrics
- Constructs a **protein interaction network** using statistical thresholds
- Fully compatible with spatial and single-cell PLA datasets

---

## Usage

```python
import SpatialproxseqClasses as SPC

# Load your PLA count matrix
# Rows: PLA pairs (e.g., CD19:CD21), Columns: spatial spots or cells
my_data = pd.read_csv("pla_count_matrix.csv", index_col=0)

# Initialize object
sprox = SPC.sproxseqObject(my_data)

# Derive protein-level summaries
sprox.compute_protein_abundance()
sprox.compute_protein_pair_counts()

#check raw data expression 
sporx.pla_count
sprox.prtoein_count
sprox.protein_pair_count

# Run statistical computations
sprox.compute_fisher_pvalues()
sprox.compute_pearson_resid()
sprox.compute_fractional_overlap(method='average') #protein abundance and protein_pair_count are in need

# Construct interaction network
sprox.construct_protein_network(fisher_ratio_threshold=0.15, fisher_pvalue_cutoff=0.05) #fisher_value and fractional overlap value are in need

# Access results
fisher_pvals = sprox.fisher_value
fractional_overlap = sprox.fractional_overlap
pearson_resid = sprox.pearson_resid
protein_network = sprox.protein_network # the network can be further visualized with networkx
```

---

## Core Methods

| Method | Description |
|--------|-------------|
| `compute_fisher_pvalues()` | Performs Fisher’s exact test across PLA pairs per spot, with FDR correction |
| `compute_fractional_overlap(method='average')` | Computes normalized interaction scores based on PLA counts and protein abundance |
| `compute_pearson_resid()` | Computes Pearson residuals between observed and expected PLA matrices |
| `compute_protein_abundance()` | Sums PLA signals for each individual protein across all pairings |
| `compute_protein_pair_counts()` | Sums counts of symmetric PLA pairs (A:B and B:A) |
| `construct_protein_network()` | Builds a NetworkX graph of significant protein-protein associations |

---

## Output Variables

After running the respective methods, the following attributes become available:

- `sprox.fisher_value`: DataFrame of adjusted Fisher p-values
- `sprox.fractional_overlap`: Normalized PLA pair interaction weights
- `sprox.pearson_resid`: Matrix of Pearson residuals
- `sprox.protein_count`: Protein-level abundance across spots
- `sprox.protein_pair_count`: Total pairwise PLA counts (A:B + B:A)
- `sprox.protein_network`: Undirected NetworkX graph object of protein interactions

---

## Applications

This package is intended for:
- **Spatial proximity interaction profiling** using Spatial Proximity Sequencing
- **Network-based inference of protein complexes**

---
