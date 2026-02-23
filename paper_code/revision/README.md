# Spatial Deconvolution Using Tangram

This folder contains code and results for spatial spot deconvolution using Tangram, mapping single-cell RNA-seq reference profiles onto Visium spatial transcriptomics data from human tonsil samples.

The goal is to infer cell-type composition for each spatial spot, which is later used to interpret Spatial Prox-Seq (PLA) signals.

---

## Directory Structure

```
revision/
├── A1_tonsil_tangram.py # Tangram mapping for sample A1
├── BCD_tonsil_tangram.py # Tangram mapping for samples B/C/D
├── Deconvolution_A1.ipynb # Notebook for A1 analysis and visualization
├── Deconvolution_B1.ipynb # Notebook for B1 analysis and visualization
├── A1_cell_map_result.h5ad # Tangram output (A1)
├── B1_cell_map_result.h5ad # Tangram output (B1)
└── README.md
```
