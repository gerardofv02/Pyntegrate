"""
CHIP-DNASE
"""

import Pyntegrate
import matplotlib.pyplot as plt
import numpy as np

bins = 100
x = np.linspace(-1000, 1000, bins)

atac = Pyntegrate.genomic_signal('mosimData/SampleDataDNase-seq-5000.bw', 'bigwig')
chip = Pyntegrate.genomic_signal('mosimData/SampleDataChip-seq-5000.bw', 'bigwig')

features, atac_signal, tsses,tsses_1kb = Pyntegrate.SeqSignalAnalysis.generate_array_simple_signal(dbPath='data/gencode.vM25.annotation.gtf.db',
                                                                                                   filePath='mosimData/SampleDataDNase-seq-5000.bw',
                                                                                                   extensionFile='bigwig',
                                                                                                   genome='mm10',
                                                                                                   bins=bins)

features, chip_signal, tsses,tsses_1kb = Pyntegrate.SeqSignalAnalysis.generate_array_simple_signal(dbPath='data/gencode.vM25.annotation.gtf.db',
                                                                                                   filePath='mosimData/SampleDataChip-seq-5000.bw',
                                                                                                   extensionFile='bigwig',
                                                                                                   genome='mm10',
                                                                                                   bins=bins)

atac_signal_good , chip_signal_good = Pyntegrate.SeqSignalAnalysis.chip_genes_not_used_with_atac(atac_signal,chip_signal)

normalized_subtracted_atac = Pyntegrate.SeqSignalAnalysis.value_array_simple(atac_signal_good)
normalized_subtracted_chip = Pyntegrate.SeqSignalAnalysis.value_array_simple(chip_signal_good)

fig = Pyntegrate.SeqSignalAnalysis.chip_dnase(normalized_subtracted_chip, normalized_subtracted_atac,xAxes=x,name="CHIP-DNASE")

plt.show()
