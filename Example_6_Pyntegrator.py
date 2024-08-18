import Pyntegrate
import pandas as pd
import matplotlib.pyplot as plt

Pyntegrate.SeqSignalAnalysis.homer_annotate_peaks(file_directory="datos/pruebaHomer/Sanchez_et_al_peaks.narrowPeak",
                                                  gene="mm10",
                                                  output="misPeaks.txt")

df = pd.read_csv("misPeaks.txt",delimiter="\t")

fig = Pyntegrate.plotutils.annotate_peaks(df)

plt.show()