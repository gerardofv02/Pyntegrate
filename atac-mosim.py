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

db = gffutils.FeatureDB(os.path.join(path, 'data/gencode.vM25.annotation.gtf.db'))
print(db)

##CArgamos el archivo
atac = Pyntegrate.genomic_signal('test_2_DNase-seq.bw', 'bigwig')
tsses = pybedtools.BedTool(tss_generator()).saveas('tsses.gtf')


tsses_1kb = tsses.slop(b=1000, genome='mm10', output='tsses-1kb.gtf')
processes = multiprocessing.cpu_count()


##Cargamos señales dnase, lo procesamos y enseñamos algunas gráficas
atac_signal =  atac.array(
    tsses_1kb,
    bins=100,
    processes=processes
)

print(atac_signal)
atac_signal_good = []

##Descartamos primeras posiciones de array (ya que devuelve arrays de 1 posición)
for x in atac_signal:
    atac_signal_good.append(x[0])

atac_signal_good = np.array(atac_signal_good)
normalized_subtracted_semi = Pyntegrate.chipSeqSignalAnalysis.value_array_simple(atac_signal)
##Normalizamos

normalized_subtracted_semi /= atac.mapped_read_count() / 1e6

##Suponiendo que aquellos valores donde la media sea 0, no existen picos con lo cual los elimino

##QUitamos aquellos sin picos (o eso creo que es así)
filtered_arrays = [arr for arr in normalized_subtracted_semi if np.mean(arr) != 0]

##COnvertimos en array de numpy para poder hacer la ordenación
normalized_subtracted_complete = np.array(filtered_arrays)
print(normalized_subtracted_semi)
x = np.linspace(-1000, 1000, bins)



##Primera figura sin ordenación
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



##Segunda y cuarta figura con ordenación de que se muestran la media de los subintervalos

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

##TErcera figurra que se muestra la señal más altas de los subintervalos

fig = Pyntegrate.plotutils.imshow(

    # These are the same arguments as above.
    normalized_subtracted_complete,
    x=x,
    figsize=(3, 7),
    vmin=5, vmax=99,  percentile=True,
    line_kwargs=dict(color='k', label='All'),
    fill_kwargs=dict(color='k', alpha=0.3),

    # This is new: sort by mean signal
    sort_by=np.argmax(normalized_subtracted_complete, axis=1)
)



fig = Pyntegrate.plotutils.imshow(
    normalized_subtracted_complete,
    x=x,
    figsize=(3, 7),
    vmin=5, vmax=99,  percentile=True,
    line_kwargs=dict(color='k', label='All'),
    fill_kwargs=dict(color='k', alpha=0.3),
    sort_by=normalized_subtracted_complete.mean(axis=1)
)

plt.show()