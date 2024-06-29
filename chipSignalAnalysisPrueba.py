import Pyntegrate
import os

path = os.path.dirname(os.path.abspath(__file__))
print(path)

features, arrays = Pyntegrate.chipSeqSignalAnalysis.generate_arrays_features_from_tsses_from_db(os.path.join(path, 'data/gencode.vM25.annotation.gtf.db'), os.path.join(path, 'bwFiles/SRR1204546_sort_nondup.bw'),os.path.join(path, 'bwFiles/SRR1204544_sort_nondup.bw'),genome = "mm10")

print("END")