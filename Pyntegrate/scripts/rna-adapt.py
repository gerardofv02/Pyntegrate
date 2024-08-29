import pandas as pd
import numpy as np

df = pd.read_csv('test_2_options_RNA-seq.csv', delimiter="\t")
print(df)

np.random.seed(42) 
gene_lengths_kb = {gene: np.random.uniform(0.5, 20.0) for gene in df['GeneID']}
total_mapped_reads_million = 50.0

df['FPKM_Group1_Time0_Rep1'] = df.apply(
    lambda row: (row['Group1.Time0.Rep1'] / (gene_lengths_kb[row['GeneID']] * total_mapped_reads_million)), axis=1
)

df['log2foldchange'] = np.log2(df['Group1.Time2.Rep1'] / df['Group1.Time0.Rep1'])

print(df.head())
df.to_csv('test_RNA-seq_w_changes.csv', sep="\t")