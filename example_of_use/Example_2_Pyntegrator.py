"""
In this example we are going to see how can we load DNASE-seq data, make some figures about it to 
analyse them and make some figures analyzing DNASE-seq data with RNA-seq data
"""

import Pyntegrate
import matplotlib.pyplot as plt
import numpy as np


bins = 100

features, atac_signal, tsses,tsses_1kb = Pyntegrate.SeqSignalAnalysis.generate_array_simple_signal(dbPath='data/gencode.vM25.annotation.gtf.db',
                                                                                                   filePath='mosimData/test_2_DNase-seq.bw',
                                                                                                   extensionFile='bigwig',
                                                                                                   genome='mm10',
                                                                                                   bins=bins)

normalized_subtracted_values = Pyntegrate.SeqSignalAnalysis.value_array_simple(atac_signal)
x = np.linspace(-1000, 1000, bins)

fig = Pyntegrate.SeqSignalAnalysis.heatmap_no_sorted(signal_values=normalized_subtracted_values, xAxis=x)
fig = Pyntegrate.SeqSignalAnalysis.heatmap_sorted_by_meanValues(signal_values=normalized_subtracted_values, xAxis=x)
fig = Pyntegrate.SeqSignalAnalysis.heatmap_sorted_by_maxValueIndex(signal_values=normalized_subtracted_values, xAxis=x)

data = Pyntegrate.results_table.DEseq2Results('mosimData/test_RNA-seq_w_changes.csv')
data.data['gene_id'] = data.data['GeneID']
data.reindex_to(tsses, attribute='gene_id')
# print(len(normalized_subtracted), len(data), "\n\n\n\n\n\n")
# print(atac_signal)
print(data)


normalized_subtracted, data = Pyntegrate.SeqSignalAnalysis.chip_or_atac_genes_not_used_with_rna(data,
                                                                                                atac_signal,
                                                                                                id="id",
                                                                                                rna_column_name="GeneID",
                                                                                                delete_no_peaks=True)

# print(len(normalized_subtracted), len(data))


values = Pyntegrate.SeqSignalAnalysis.value_array_simple(normalized_subtracted)

fig = Pyntegrate.SeqSignalAnalysis.atac_or_chip_with_rna(signal_values=values, rna=data, xAxis=x)
plt.show()