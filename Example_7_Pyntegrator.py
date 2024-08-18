import Pyntegrate
import os
import numpy as np

bins = 100
x = np.linspace(-1000, 1000, bins)
features, arrays, tsses, tsses_1kb = Pyntegrate.SeqSignalAnalysis.generate_arrays_features_from_tsses_from_db(dbPath='data/gencode.vM25.annotation.gtf.db', 
                                                                                                              ipSignalPath='realData/SRR1204544_sort_nondup.bw', 
                                                                                                              extensionIp='bigwig',
                                                                                                              inputSignalPath='realData/SRR1204546_sort_nondup.bw',
                                                                                                              extensionInput='bigwig', 
                                                                                                              genome='mm10', 
                                                                                                              bins=bins)

normmalized_subtracted = Pyntegrate.SeqSignalAnalysis.calculate_peaks_with_gene_name(arrays['ip'], arrays['input'])

Pyntegrate.SeqSignalAnalysis.array_gene_name_and_annotation(normmalized_subtracted, x)