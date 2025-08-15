# Spatial Proximity Sequencing Maps Developmental Dynamics in the Germinal Center

**Huili Wang, Junjie Xia, et al.**

This folder contains code to reproduce the analyses presented in the paper using **three samples**.  
Analyses include spatial clustering and annotation, differential expression analysis, pseudotime inference, and Gene Ontology (GO) enrichment.

---

## 📂 Repository Structure

### 1. Annotation & Preprocessing
- **A1_annotation.ipynb** – Annotation workflow for sample A1.
- **B1_annotation.ipynb** – Annotation workflow for sample B1.
- **D1_annotation.ipynb** – Annotation workflow for sample D1.

### 2. Cross-sample Consistency
- **A1_B1_D1_code2_consistency_sample_correlation_RNA&protein.ipynb**  
  Computes correlation of RNA and protein expression patterns across the three samples.

### 3. Interaction Quantification
- **A1_B1_D1_code3_fraction_occupation_calculation.ipynb**  
  Calculates fractional occupation scores for detected protein pairs.
- **A1_B1_D1_code4_fisher_test_samples_sprox-seq.ipynb**  
  Performs Fisher’s Exact Test to identify significant protein–protein interactions across samples.

### 4. Germinal Center Subset Analyses
- **A1_B1_D1_code5_LZ_DZ_protein_analysis.ipynb**  
  Compares protein expression between Light Zone (LZ) and Dark Zone (DZ) B cells.
- **A1_B1_D1_code6.1_LZ_DZ_Follicle_subset.ipynb**  
  Extracts follicle subsets for downstream analysis.
- **A1_B1_D1_code6.2_protein_pairs_pass_filter_B_cell_cluster.ipynb**  
  Filters and retains significant protein pairs within B cell clusters.

### 5. Pseudotime Analysis
- **A1_B1_D1_code6.3_pseudotime_analysis.ipynb**  
  Infers developmental trajectories of B cell subsets.

### 6. Gene Set Enrichment Analysis (GSEA)
- **A1_B1_D1_code6.4_GSEA_analysis_CD20_CD32.ipynb**  
  GSEA for CD20–CD32 interactions.
- **A1_B1_D1_code6.4_GSEA_analysis_CD21_CD35.ipynb**  
  GSEA for CD21–CD35 interactions.
- **A1_B1_D1_code6.4_GSEA_analysis_ITGA4_VCAM1.ipynb**  
  GSEA for ITGA4–VCAM1 interactions.

### 7. Differential Expression & GO Enrichment
- **DE_GO_analysis.py**  
  Script for performing differential expression analysis and GO enrichment.

---

## 📊 Analysis Workflow Overview

1. **Annotation**  
   Assigns spatial clusters to known tissue or cell-type identities.

2. **Cross-sample Consistency**  
   Compares RNA and protein-level correlations across three samples.

3. **Interaction Quantification**  
   Detects and quantifies significant protein–protein interactions via Fisher’s Exact Test and fractional occupation.

4. **Zone-specific Analyses**  
   Focuses on LZ vs DZ differences in protein expression and interaction patterns.

5. **Pseudotime Inference**  
   Models developmental progression of B cell populations.

6. **Functional Enrichment**  
   Runs GSEA and GO analysis to interpret biological significance of detected patterns.

---

## 🛠️ Requirements

- Python ≥ 3.8
- Jupyter Notebook
- Key packages:
  scanpy
  anndata
  pandas
  numpy
  matplotlib
  seaborn
  statsmodels
  gseapy
