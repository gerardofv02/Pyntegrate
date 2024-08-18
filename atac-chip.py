import Pyntegrate
import os
import gffutils
import pybedtools
from pybedtools.featurefuncs import TSS
from gffutils.helpers import asinterval
import multiprocessing
from matplotlib import pyplot as plt
import numpy as np



# def tss_generator():
#     """
#     Generator function to yield TSS of each annotated transcript
#     """
#     for transcript in db.features_of_type('transcript'):
#         yield TSS(asinterval(transcript), upstream=1, downstream=0)

path = os.path.dirname(os.path.abspath(__file__))
print(path)
bins = 100

# db = gffutils.FeatureDB(os.path.join(path, 'data/gencode.vM25.annotation.gtf.db'))
# print(db)

##CArgamos el archivo
atac = Pyntegrate.genomic_signal('mosimData/SampleDataDNase-seq-5000.bw', 'bigwig')
chip = Pyntegrate.genomic_signal('mosimData/SampleDataChip-seq-5000.bw', 'bigwig')
# tsses = pybedtools.BedTool(tss_generator()).saveas('tsses.gtf')


# tsses_1kb = tsses.slop(b=1000, genome='mm10', output='tsses-1kb.gtf')

tsses = pybedtools.BedTool('tsses.gtf')
tsses_1kb = pybedtools.BedTool('tsses-1kb.gtf')
processes = multiprocessing.cpu_count()


##Cargamos señales dnase, lo procesamos y enseñamos algunas gráficas
atac_signal =  atac.array(
    tsses_1kb,
    bins=100,
    processes=processes
)

chip_signal =  chip.array(
    tsses_1kb,
    bins=100,
    processes=processes
)

print(atac_signal, chip_signal)


# filtered_arrays_atac = [arr for arr in atac_signal if np.mean(arr[0]['values']) != 0]
# filtered_arrays_chip = [arr for arr in chip_signal if np.mean(arr[0]['values']) != 0]

# print("\nfiiltered atac-len: ", len(filtered_arrays_atac))
# print("\nfiiltered chip-len: ", len(filtered_arrays_chip))

atac_signal_good = []
chip_signal_good = []

for x in atac_signal:
    for y in chip_signal:
        if x[0]['gene_name'] == y[0]['gene_name']:
            # print(x[0]['id'])
            atac_signal_good.append(x)
            chip_signal_good.append(y)
            break

print("\nLen de atacsignal : ", len(atac_signal_good))
print("\nlen de chip signal : ", len(chip_signal_good))
normalized_subtracted_semi_atac = Pyntegrate.chipSeqSignalAnalysis.value_array_simple(atac_signal_good)
##Normalizamos

normalized_subtracted_semi_atac /= atac.mapped_read_count() / 1e6

##Suponiendo que aquellos valores donde la media sea 0, no existen picos con lo cual los elimino

##QUitamos aquellos sin picos (o eso creo que es así)
# filtered_arrays_atac = [arr for arr in normalized_subtracted_semi_atac if np.mean(arr) != 0]

##COnvertimos en array de numpy para poder hacer la ordenación
normalized_subtracted_complete_atac = np.array(normalized_subtracted_semi_atac)
# print(len(normalized_subtracted_complete_atac))





##Cargamos chip-seq, lo procesamos y enseñamos algunas gráficas


normalized_subtracted_semi_chip = Pyntegrate.chipSeqSignalAnalysis.value_array_simple(chip_signal_good)

normalized_subtracted_semi_chip /= chip.mapped_read_count() / 1e6

##Suponiendo que aquellos valores donde la media sea 0, no existen picos con lo cual los elimino

# filtered_arrays_chip = [arr for arr in normalized_subtracted_semi_chip if np.mean(arr) != 0]

normalized_subtracted_complete_chip = np.array(normalized_subtracted_semi_chip)
# print(len(normalized_subtracted_complete_chip))

fig = Pyntegrate.chipSeqSignalAnalysis.chip_dnase(normalized_subtracted_complete_chip, normalized_subtracted_complete_atac,bins)

plt.show()