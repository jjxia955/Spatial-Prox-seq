"""
Author: Junjie Xia
Email: jjxia@uchicago.edu

Computational framework for analyzing spatial prox-seq data, including calculation of 
Fisher Exact Test, Pearson Residual, and Fractional Overlap value.
"""

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import networkx as nx
import matplotlib.pyplot as plt

class sproxseqObject:
    def __init__(self, pla_count: pd.DataFrame, sep: str = ":"):
        """
        Initialize the spatial prox-seq object with PLA count matrix and delimiter.

        Parameters:
        - pla_count: DataFrame of shape (PLA products x spots/samples)
        - sep: Delimiter used in PLA product names (default is ':')
        """
        self.pla_count = pla_count
        self.sep = sep

    def compute_fisher_pvalues(self):
        """
        Compute one-sided Fisher's exact test for each PLA product across cells,
        and apply BH correction per column (spot).

        Result is stored in self.fisher_value.
        """
        fisher_p = pd.DataFrame(np.nan, index=self.pla_count.index, columns=self.pla_count.columns)
        probeA = np.array([s.split(self.sep)[0] for s in fisher_p.index])
        probeB = np.array([s.split(self.sep)[1] for s in fisher_p.index])

        for i in fisher_p.index:
            tempA, tempB = i.split(self.sep)
            x00 = (probeA == tempA) & (probeB == tempB)
            x01 = (probeA == tempA) & (probeB != tempB)
            x10 = (probeA != tempA) & (probeB == tempB)
            x11 = (probeA != tempA) & (probeB != tempB)

            temp01 = self.pla_count.loc[x01,:].sum(axis=0)
            temp10 = self.pla_count.loc[x10,:].sum(axis=0)
            temp11 = self.pla_count.loc[x11,:].sum(axis=0)

            for j in fisher_p.columns:
                x00_val = self.pla_count.loc[x00, j].sum()
                fisher_p.at[i, j] = stats.fisher_exact(
                    [[x00_val, temp01[j]], [temp10[j], temp11[j]]],
                    alternative='greater')[1]

        corrected = pd.DataFrame(np.nan, index=fisher_p.index, columns=fisher_p.columns)
        for col in corrected.columns:
            corrected[col] = multipletests(fisher_p[col], method='fdr_bh')[1]

        self.fisher_value = corrected

    def compute_fractional_overlap(self, method: str = 'average', save_path: str = None):
        """
        Compute fractional overlap score for each PLA pair based on protein abundance.

        If protein abundance or protein pair counts are not computed, they will be computed.

        Parameters:
        - method: 'average' or 'abundance'
        - save_path: Optional path to save triplet format output

        Result is stored in self.fractional_overlap.
        """
        if not hasattr(self, 'protein_count'):
            print("[INFO] protein_count not found. Computing protein abundance...")
            self.compute_protein_abundance()

        if not hasattr(self, 'protein_pair_count'):
            print("[INFO] protein_pair_count not found. Computing protein pair counts...")
            self.compute_protein_pair_counts()

        protein_abundance = self.protein_count
        pla_weights = pd.DataFrame(index=self.protein_pair_count.index, columns=self.protein_pair_count.columns)

        for pair in self.protein_pair_count.index:
            proteinA, proteinB = pair.split(self.sep)
            if proteinA not in protein_abundance.index or proteinB not in protein_abundance.index:
                continue

            weightA = self.protein_pair_count.loc[pair] / protein_abundance.loc[proteinA]
            weightB = self.protein_pair_count.loc[pair] / protein_abundance.loc[proteinB]

            if method == 'abundance':
                total = protein_abundance.loc[proteinA] + protein_abundance.loc[proteinB]
                pla_weights.loc[pair] = weightA * (protein_abundance.loc[proteinA] / total) + \
                                        weightB * (protein_abundance.loc[proteinB] / total)
            else:
                pla_weights.loc[pair] = (weightA + weightB) / 2

        if save_path:
            with open(save_path, "w") as f:
                for cell in pla_weights.columns:
                    for pair in pla_weights.index:
                        proteinA, proteinB = pair.split(self.sep)
                        f.write(f"{proteinA}\t{proteinB}\t{pla_weights.at[pair, cell]}\n")

        self.fractional_overlap = pla_weights

    def compute_pearson_resid(self):
        """
        Compute Pearson residuals for each cell using observed and expected PLA matrix.
        Residuals are flattened to {proteinA:proteinB} format.

        Result is stored in self.pearson_residual.
        """
        probes_A = [p.split(self.sep)[0] for p in self.pla_count.index]
        probes_B = [p.split(self.sep)[1] for p in self.pla_count.index]
        residuals = pd.DataFrame(index=self.pla_count.index, columns=self.pla_count.columns)

        for cell in self.pla_count.columns:
            cell_counts = self.pla_count[cell]
            matrix = pd.DataFrame({
                'ProbeA': probes_A,
                'ProbeB': probes_B,
                'Count': cell_counts.values
            }).pivot_table(index='ProbeA', columns='ProbeB', values='Count', fill_value=0)

            row_sums = matrix.sum(axis=1)
            col_sums = matrix.sum(axis=0)
            grand_total = matrix.values.sum()

            if grand_total == 0:
                residuals[cell] = 0
                continue

            expected = np.outer(row_sums, col_sums) / grand_total
            observed = matrix.values

            with np.errstate(divide='ignore', invalid='ignore'):
                pearson_resid = np.where(expected > 0, (observed - expected) / np.sqrt(expected), 0)

            residual_matrix = pd.DataFrame(pearson_resid, index=matrix.index, columns=matrix.columns)
            residual_flat = residual_matrix.stack()
            residual_flat.index = residual_flat.index.map(lambda x: f"{x[0]}{self.sep}{x[1]}")
            residuals[cell] = residual_flat

        self.pearson_resid = residuals

    def compute_protein_pair_counts(self):
        """
        Compute combined counts of probeA:probeB and probeB:probeA for all possible protein pairs.
        Result is stored in self.protein_pair_count.
        """
        AB1 = np.array([s.split(self.sep)[0] for s in self.pla_count.index])
        AB2 = np.array([s.split(self.sep)[1] for s in self.pla_count.index])
        A_unique = np.unique(AB1)
        B_unique = np.unique(AB2)
        new_index = [f'{A_unique[i]}{self.sep}{B_unique[j]}' for i in range(len(A_unique)) for j in range(i, len(B_unique))]
        self.protein_pairs = new_index
        self.protein_pair_count = pd.DataFrame(0, index=new_index, columns=self.pla_count.columns)

        for i in self.protein_pair_count.index:
            probeA, probeB = i.split(self.sep)
            try:
                entry1 = self.pla_count.loc[f'{probeA}{self.sep}{probeB}', :]
            except KeyError:
                entry1 = 0
            try:
                entry2 = self.pla_count.loc[f'{probeB}{self.sep}{probeA}', :]
            except KeyError:
                entry2 = 0
            self.protein_pair_count.loc[i, :] = entry1 + entry2

    def compute_protein_abundance(self):
        """
        Compute total protein abundance by summing PLA counts from both probeA and probeB.
        Result is stored in self.protein_count.
        """
        AB1 = np.array([s.split(self.sep)[0] for s in self.pla_count.index])
        AB2 = np.array([s.split(self.sep)[1] for s in self.pla_count.index])
        AB_unique = np.unique(np.concatenate((AB1, AB2)))
        AB_unique.sort()
        self.proteins = list(AB_unique)
        self.protein_count = pd.DataFrame(0, index=AB_unique, columns=self.pla_count.columns)

        for protein in self.protein_count.index:
            counts_A = self.pla_count.loc[AB1 == protein, :].sum(axis=0)
            counts_B = self.pla_count.loc[AB2 == protein, :].sum(axis=0)
            self.protein_count.loc[protein, :] = counts_A + counts_B
            
            
    def construct_protein_network(self, fisher_ratio_threshold: float = 0.15, fisher_pvalue_cutoff: float = 0.05):
        """
        Construct a symmetric protein interaction network using average fractional overlap matrix,
        filtered by ratio of significant Fisher test p-values.

        Parameters:
        - fisher_ratio_threshold: minimum fraction of spots with Fisher p-value < cutoff to retain a pair
        - fisher_pvalue_cutoff: cutoff value for Fisher p-values (default 0.05)

        Result is stored in self.protein_network, which can be further visualized with networkx.
        """
        if not hasattr(self, 'fractional_overlap'):
            raise ValueError("fractional_overlap not found. Please run compute_fractional_overlap() first.")
        if not hasattr(self, 'fisher_value'):
            raise ValueError("fisher_value not found. Please run compute_fisher_pvalues() first.")

        avg_weight = pd.DataFrame({'mean': self.fractional_overlap.mean(axis=1)})
        avg_weight['probeA'] = [s.split(self.sep)[0] for s in avg_weight.index]
        avg_weight['probeB'] = [s.split(self.sep)[1] for s in avg_weight.index]
        avg_weight = avg_weight.pivot(index='probeA', columns='probeB', values='mean')
 
        upper_triangle = np.triu(avg_weight.fillna(0))
        avg_weight = pd.DataFrame(upper_triangle, index=avg_weight.index, columns=avg_weight.columns)
        avg_weight = avg_weight + avg_weight.T - np.diag(np.diag(avg_weight))

        fisher_ratio = pd.DataFrame({'fisher': (self.fisher_value < fisher_pvalue_cutoff).sum(axis=1) / self.fisher_value.shape[1]})
        fisher_ratio['probeA'] = [s.split(self.sep)[0] for s in fisher_ratio.index]
        fisher_ratio['probeB'] = [s.split(self.sep)[1] for s in fisher_ratio.index]
        fisher_ratio = fisher_ratio.pivot(index='probeA', columns='probeB', values='fisher')

        common_index = fisher_ratio.index.intersection(avg_weight.index)
        common_columns = fisher_ratio.columns.intersection(avg_weight.columns)
        fisher_common = fisher_ratio.loc[common_index, common_columns]
        weight_common = avg_weight.loc[common_index, common_columns]

        mask = fisher_common <= fisher_ratio_threshold
        weight_common_masked = weight_common.mask(mask, 0)
        avg_weight.update(weight_common_masked)
        avg_weight_sym = avg_weight.combine(avg_weight.T, func=np.fmax)

        my_order = list(avg_weight_sym.index)

        G = nx.Graph()
        G.add_nodes_from(my_order)

        for i, protein1 in enumerate(my_order):
            for j, protein2 in enumerate(my_order):
                weight = avg_weight_sym.loc[protein1, protein2]
                if weight > 0:
                    if protein1 == protein2:
                        color = 'red'
                    elif i < j:
                        color = 'blue'
                    else:
                        continue
                    G.add_edge(protein1, protein2, weight=weight, color=color)

        self.protein_network = G

        print(f"[INFO] Network constructed with {G.number_of_nodes()} proteins and {G.number_of_edges()} interactions.")