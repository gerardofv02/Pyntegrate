# import Pyntegrate
# import os

# path = os.path.dirname(os.path.abspath(__file__))
# print(path)

# features, arrays ,tsses , tsses_1kb = Pyntegrate.chipSeqSignalAnalysis.generate_arrays_features_from_tsses_from_db(os.path.join(path, 'data/gencode.vM25.annotation.gtf.db'), os.path.join(path, 'BamFiles/SRR1204544_sort_nondup.bw'),'bigwig',os.path.join(path, 'BamFiles/SRR1204546_sort_nondup.bw'),'bigwig',genome = "mm10")
# comprobacion = []
# genes = []
# misarrays = arrays['ip']
# for x in tsses_1kb:
#     genes.append(x[8].split(";")[3].split('"')[1])

# for x in range(len(genes)):
#     if genes[x] == misarrays[x][0]['gene_name']:
#         comprobacion.append(True)

#     else:
#         comprobacion.append(False)
# print("\nArrays: ",arrays['ip'], "\n Tsses: ", type(tsses_1kb), "\nAll true?: ", all(comprobacion))

import Pyntegrate

Pyntegrate.chipSeqSignalAnalysis.bigwigToBed('test_3_ChIP-seq_3.bw','test_3_ChIP-seq_3.bed')