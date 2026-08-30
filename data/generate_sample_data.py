"""
Generates a synthetic gene expression dataset for the Gene Expression Data
Explorer workshop app.

Design:
    - 300 genes, named with realistic-looking symbols (mix of real human gene
      symbols and generated ones) so the demo feels authentic.
    - 12 samples: 6 "Control" and 6 "Treatment".
    - Expression values are simulated read counts (log2-normalized, like
      typical RNA-seq TPM/CPM data after log transform).
    - ~15% of genes are seeded as differentially expressed (up or down in
      Treatment) so volcano plots / heatmaps have real signal to show.

Output: data/gene_expression_sample.csv
    Columns: GeneID, GeneSymbol, Control_1..Control_6, Treatment_1..Treatment_6
"""

import numpy as np
import pandas as pd

rng = np.random.default_rng(42)

N_GENES = 300
N_CONTROL = 6
N_TREATMENT = 6

# A pool of real-looking human gene symbols (common well-known genes +
# plausible-looking synthetic ones) so the dataset feels like real biology.
KNOWN_GENES = [
    "TP53", "BRCA1", "BRCA2", "EGFR", "MYC", "KRAS", "PTEN", "RB1", "APC",
    "PIK3CA", "AKT1", "STAT3", "VEGFA", "TNF", "IL6", "IL10", "CDKN2A",
    "MTOR", "NOTCH1", "WNT1", "GATA3", "FOXP3", "CD44", "CD8A", "CD4",
    "ACTB", "GAPDH", "HPRT1", "B2M", "RPL13A", "HSP90AA1", "CASP3",
    "BAX", "BCL2", "MDM2", "CCND1", "CDK4", "E2F1", "SMAD4", "TGFB1",
    "JUN", "FOS", "NFKB1", "HIF1A", "SOD1", "CAT", "GPX1", "NOS2",
    "COL1A1", "FN1", "MMP9", "TIMP1", "VIM", "SNAI1", "CDH1", "ESR1",
    "PGR", "ERBB2", "MKI67", "PCNA",
]

symbols = list(KNOWN_GENES)
i = 1
while len(symbols) < N_GENES:
    symbols.append(f"GENE{i:04d}")
    i += 1
symbols = symbols[:N_GENES]
rng.shuffle(symbols)

gene_ids = [f"ENSG{100000 + idx}" for idx in range(N_GENES)]

# Baseline mean expression per gene (log2 scale, roughly 2-14, like log2(TPM+1))
baseline = rng.uniform(2, 12, size=N_GENES)

control_cols = [f"Control_{i+1}" for i in range(N_CONTROL)]
treatment_cols = [f"Treatment_{i+1}" for i in range(N_TREATMENT)]

data = np.zeros((N_GENES, N_CONTROL + N_TREATMENT))

# Seed differential expression for ~15% of genes
n_de = int(N_GENES * 0.15)
de_indices = rng.choice(N_GENES, size=n_de, replace=False)
de_direction = rng.choice([-1, 1], size=n_de)  # down or up in treatment
de_effect = rng.uniform(1.5, 4.0, size=n_de)  # log2 fold change magnitude

fold_change = np.zeros(N_GENES)
for idx, direction, effect in zip(de_indices, de_direction, de_effect):
    fold_change[idx] = direction * effect

biological_noise_sd = rng.uniform(0.3, 0.9, size=N_GENES)

for g in range(N_GENES):
    ctrl_vals = rng.normal(loc=baseline[g], scale=biological_noise_sd[g], size=N_CONTROL)
    treat_vals = rng.normal(
        loc=baseline[g] + fold_change[g], scale=biological_noise_sd[g], size=N_TREATMENT
    )
    data[g, :N_CONTROL] = ctrl_vals
    data[g, N_CONTROL:] = treat_vals

data = np.clip(data, 0, None)
data = np.round(data, 3)

df = pd.DataFrame(data, columns=control_cols + treatment_cols)
df.insert(0, "GeneSymbol", symbols)
df.insert(0, "GeneID", gene_ids)

out_path = "data/gene_expression_sample.csv"
df.to_csv(out_path, index=False)
print(f"Wrote {len(df)} genes x {N_CONTROL + N_TREATMENT} samples -> {out_path}")
