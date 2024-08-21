"""
All signal together
"""

import Pyntegrate
import matplotlib.pyplot as plt
import numpy as np
bins = 100
x = np.linspace(-1000, 1000, bins)

atac = Pyntegrate.genomic_signal('mosimData/SampleDataDNase-seq-5000.bw', 'bigwig')
chip = Pyntegrate.genomic_signal('mosimData/SampleDataChip-seq-5000.bw', 'bigwig')

features, atac_signal, tsses,tsses_1kb = Pyntegrate.SeqSignalAnalysis.generate_array_simple_signal(dbPath='data/gencode.vM25.annotation.gtf.db',
                                                                                                   filePath='mosimData/SampleDataDNase-seq-10000.bw',
                                                                                                   extensionFile='bigwig',
                                                                                                   genome='mm10',
                                                                                                   bins=bins)

features, chip_signal, tsses,tsses_1kb = Pyntegrate.SeqSignalAnalysis.generate_array_simple_signal(dbPath='data/gencode.vM25.annotation.gtf.db',
                                                                                                   filePath='mosimData/SampleDataChIP-seq-10000.bw',
                                                                                                   extensionFile='bigwig',
                                                                                                   genome='mm10',
                                                                                                   bins=bins)




data = Pyntegrate.results_table.DEseq2Results('mosimData/def_lfc_results_grupo1_t0_vs_t2.csv', delimiter=",")
data.data['GeneID'] = data.index
data.data['gene_id'] = data.data['GeneID']
data.data['log2foldchange'] = data.data['log2FoldChange']
data.reindex_to(tsses, attribute='gene_id')


atac_signal_good, data = Pyntegrate.SeqSignalAnalysis.chip_or_atac_genes_not_used_with_rna(data,
                                                                                                atac_signal,
                                                                                                id="id",
                                                                                                rna_column_name="GeneID",
                                                                                                delete_no_peaks=True)



chip_signal_good, data = Pyntegrate.SeqSignalAnalysis.chip_or_atac_genes_not_used_with_rna(data,
                                                                                                chip_signal,
                                                                                                id="id",
                                                                                                rna_column_name="GeneID",
                                                                                                delete_no_peaks=True)

normalized_subtracted_atac , normalized_subtracted_chip = Pyntegrate.SeqSignalAnalysis.chip_genes_not_used_with_atac(atac_signal_good,chip_signal_good,by="id")

print(len(normalized_subtracted_atac),len(normalized_subtracted_chip),len(data))
normalized_subtracted_atac, normalized_subtracted_chip =Pyntegrate.SeqSignalAnalysis.values_array(normalized_subtracted_atac, normalized_subtracted_chip)
fig = Pyntegrate.SeqSignalAnalysis.chip_dnase(normalized_subtracted_chip, normalized_subtracted_atac,xAxes=x,name="", xlabel="", ylabel="")
fig = Pyntegrate.SeqSignalAnalysis.all_signal_together(normalized_subtracted_chip, normalized_subtracted_atac,data, xAxis=x)
fig = Pyntegrate.SeqSignalAnalysis.chip_atac_rna_diverse(normalized_subtracted_chip, normalized_subtracted_atac,data,xAxis=x)

plt.show()