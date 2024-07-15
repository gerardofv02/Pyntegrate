import Pyntegrate
import os

# Pyntegrate.chipSeqSignalAnalysis.create_homer_tag_directory(name="miPrueba1",file="SRR1204544_sort_nondup.bam")
# directorio_completo = os.path.abspath("miPrueba1")
# print(directorio_completo)
# Pyntegrate.chipSeqSignalAnalysis.homer_peaks(tag_directory="miPrueba1")
directorio_completo = os.path.abspath("Sanchez_et_al_peaks.narrowPeak")
print(directorio_completo)
Pyntegrate.chipSeqSignalAnalysis.homer_annotate_peaks(file_directory=directorio_completo,gene="mm10", output="misPeaks.txt")