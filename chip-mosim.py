import Pyntegrate
import os
import gffutils
import pybedtools
from pybedtools.featurefuncs import TSS
from gffutils.helpers import asinterval
import multiprocessing
from matplotlib import pyplot as plt
import numpy as np



def tss_generator():
    """
    Generator function to yield TSS of each annotated transcript
    """
    for transcript in db.features_of_type('transcript'):
        yield TSS(asinterval(transcript), upstream=1, downstream=0)

path = os.path.dirname(os.path.abspath(__file__))
print(path)
bins = 100
# features, arrays ,tsses , tsses_1kb = Pyntegrate.chipSeqSignalAnalysis.generate_arrays_features_from_tsses_from_db(os.path.join(path, 'data/gencode.vM25.annotation.gtf.db'), os.path.join(path, 'BamFiles/SRR1204544_sort_nondup.bw'),'bigwig',os.path.join(path, 'BamFiles/SRR1204546_sort_nondup.bw'),'bigwig',genome = "mm10", bins=bins)

# print("\nArrays: ",len(arrays['ip']), len(arrays['input']))
# # # db = gffutils.create_db(data=path+"/Homo_sapiens.GRCh38.109.gtf",dbfn=path+"/data/Homo_sapiens.GRCh38.109.gtf.db")
# # db = gffutils.create_db(data=path+"/prueba.gtf",dbfn=path+"/data/prueba.gtf.db")

db = gffutils.FeatureDB(os.path.join(path, 'data/gencode.vM25.annotation.gtf.db'))
print(db)
chip = Pyntegrate.genomic_signal('mosimData/test_3_ChIP-seq_3.bw', 'bigwig')
tsses = pybedtools.BedTool(tss_generator()).saveas('tsses.gtf')


tsses_1kb = tsses.slop(b=1000, genome='mm10', output='tsses-1kb.gtf')
processes = multiprocessing.cpu_count()


##Cargamos chip-seq, lo procesamos y enseñamos algunas gráficas
chip_signal =  chip.array(
    tsses_1kb,
    bins=100,
    processes=processes
)

print(chip_signal)
chip_signal_good = []


chip_signal_good = np.array(chip_signal_good)
normalized_subtracted_semi = Pyntegrate.chipSeqSignalAnalysis.value_array_simple(chip_signal)

normalized_subtracted_semi /= chip.mapped_read_count() / 1e6

##Suponiendo que aquellos valores donde la media sea 0, no existen picos con lo cual los elimino

filtered_arrays = [arr for arr in normalized_subtracted_semi if np.mean(arr) != 0]

normalized_subtracted_complete = np.array(filtered_arrays)
print(normalized_subtracted_semi)
x = np.linspace(-1000, 1000, bins)

fig = Pyntegrate.chipSeqSignalAnalysis.chip_dnase(normalized_subtracted_complete, normalized_subtracted_complete,bins)


fig = Pyntegrate.plotutils.imshow(

    # The array to plot; here, we've subtracted input from IP.
    normalized_subtracted_complete,

    # X-axis to use
    x=x,

    # Change the default figure size to something smaller for this example
    figsize=(3, 7),

    # Make the colorbar limits go from 5th to 99th percentile.
    # `percentile=True` means treat vmin/vmax as percentiles rather than
    # actual values.
    percentile=True,
    vmin=5,
    vmax=99,

    # Style for the average line plot (black line)
    line_kwargs=dict(color='k', label='All'),

    # Style for the +/- 95% CI band surrounding the
    # average line (transparent black)
    fill_kwargs=dict(color='k', alpha=0.3),
)



fig = Pyntegrate.plotutils.imshow(

    # These are the same arguments as above.
    normalized_subtracted_complete,
    x=x,
    figsize=(3, 7),
    vmin=5, vmax=99,  percentile=True,
    line_kwargs=dict(color='k', label='All'),
    fill_kwargs=dict(color='k', alpha=0.3),

    # This is new: sort by mean signal
    sort_by=normalized_subtracted_complete.mean(axis=1)
)

# fig = Pyntegrate.plotutils.imshow(

#     # These are the same arguments as above.
#     normalized_subtracted_complete,
#     x=x,
#     figsize=(3, 7),
#     vmin=5, vmax=99,  percentile=True,
#     line_kwargs=dict(color='k', label='All'),
#     fill_kwargs=dict(color='k', alpha=0.3),

#     # This is new: sort by mean signal
#     sort_by=np.argmax(normalized_subtracted_complete, axis=1)
# )



# fig = Pyntegrate.plotutils.imshow(
#     normalized_subtracted_complete,
#     x=x,
#     figsize=(3, 7),
#     vmin=5, vmax=99,  percentile=True,
#     line_kwargs=dict(color='k', label='All'),
#     fill_kwargs=dict(color='k', alpha=0.3),
#     sort_by=normalized_subtracted_complete.mean(axis=1)
# )

plt.show()