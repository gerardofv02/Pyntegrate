"""
In this example we are going to see how can we load CHIP-seq data, make some figures about it to 
analyse them and make some figures analyzing CHIP-seq data with RNA-seq data
"""

import os
import gffutils
import Pyntegrate
bins = 100
features, arrays, tsses, tsses_1kb = Pyntegrate.SeqSignalAnalysis.generate_arrays_features_from_tsses_from_db(dbPath='/home/jerry/Documentos/Python3.10/prueba/data/gencode.vM25.annotation.gtf.db', 
                                                                                                              ipSignalPath='/home/jerry/Documentos/Python3.10/prueba/realData/SRR1204544_sort_nondup.bw', 
                                                                                                              extensionIp='bigwig',
                                                                                                              inputSignalPath='/home/jerry/Documentos/Python3.10/prueba/realData/SRR1204546_sort_nondup.bw',
                                                                                                              extensionInput='bigwig', 
                                                                                                              genome='mm10', 
                                                                                                              bins=bins)

print(len(arrays['ip']), len(arrays['input']))

arrays_ip_values,arrays_input_values = Pyntegrate.SeqSignalAnalysis.values_array(array_ip=arrays['ip'], array_input=arrays['input'])

import numpy as np
import matplotlib.pyplot as plt
x = np.linspace(-1000, 1000, bins)
fig = Pyntegrate.SeqSignalAnalysis.distance_from_tss_chipSeq(arrays_ip=arrays_ip_values, arrays_input=arrays_input_values,xAxes=x, name="MiGRafica")

normalized_subtracted = Pyntegrate.SeqSignalAnalysis.calculate_peaks_with_gene_name(arrays_ip=arrays['ip'], arrays_input=arrays['input'])
normalized_subtracted_values = Pyntegrate.SeqSignalAnalysis.calculate_peaks(arrays_ip=arrays_ip_values, arrays_input=arrays_input_values)

fig = Pyntegrate.SeqSignalAnalysis.heatmap_no_sorted(signal_values=normalized_subtracted_values, xAxis=x)
fig = Pyntegrate.SeqSignalAnalysis.heatmap_sorted_by_meanValues(signal_values=normalized_subtracted_values, xAxis=x)
fig = Pyntegrate.SeqSignalAnalysis.heatmap_sorted_by_maxValueIndex(signal_values=normalized_subtracted_values, xAxis=x)

data = Pyntegrate.results_table.DEseq2Results('/home/jerry/Documentos/Python3.10/prueba/realData/desq2_Sanchez.csv')
data.reindex_to(tsses, attribute='gene_name')
# print(len(normalized_subtracted), len(data), "\n\n\n\n\n\n")

normalized_subtracted, data = Pyntegrate.SeqSignalAnalysis.chip_or_atac_genes_not_used_with_rna(data,
                                                                                                normalized_subtracted,
                                                                                                id="gene_name",
                                                                                                rna_column_name="gene_name",
                                                                                                delete_no_peaks=False)

# print(len(normalized_subtracted), len(data))
data.data['log2foldchange'] = data.data['log2FoldChange']

values = Pyntegrate.SeqSignalAnalysis.value_array_simple(normalized_subtracted)

fig = Pyntegrate.SeqSignalAnalysis.atac_or_chip_with_rna(signal_values=values, rna=data, xAxis=x)
plt.show()