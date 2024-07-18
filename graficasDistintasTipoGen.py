import Pyntegrate
import os

# Pyntegrate.chipSeqSignalAnalysis.bigwigToBed("SRR1204544_sort_nondup.bw", "SRR1204544_sort_nondup.bw.bed")
# Pyntegrate.chipSeqSignalAnalysis.bigwigToBed("SRR1204546_sort_nondup.bw", "SRR1204546_sort_nondup.bw.bed")

# Pyntegrate.chipSeqSignalAnalysis.homer_annotate_peaks(file_directory="SRR1204544_sort_nondup.bw.bed",gene="mm10", output=("SRR1204544_sort_nondup.bw.bed.txt"))
# Pyntegrate.chipSeqSignalAnalysis.homer_annotate_peaks(file_directory="SRR1204546_sort_nondup.bw.bed",gene="mm10", output=("SRR1204546_sort_nondup.bw.bed.txt"))

path = os.path.dirname(os.path.abspath(__file__))
print(path)

features, arrays ,tsses , tsses_1kb = Pyntegrate.chipSeqSignalAnalysis.generate_arrays_features_from_tsses_from_db(os.path.join(path, 'data/gencode.vM25.annotation.gtf.db'), os.path.join(path, 'BamFiles/SRR1204544_sort_nondup.bw'),'bigwig',os.path.join(path, 'BamFiles/SRR1204546_sort_nondup.bw'),'bigwig',genome = "mm10")


Pyntegrate.chipSeqSignalAnalysis.array_gene_name_and_annotation(bedFile_ip="SRR1204544_sort_nondup.bw.bed",
                                                                 bedFile_input="SRR1204546_sort_nondup.bw.bed",
                                                                 gene="mm10", arrays_ip=arrays['ip'], arrays_input=arrays['input'])