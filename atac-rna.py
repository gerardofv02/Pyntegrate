# import pyBigWig
# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt




# # Abrir el archivo BigWig existente
# bw = pyBigWig.open('test_2_DNase-seq.bw')

# # Imprimir un resumen del archivo BigWig
# print("Resumen del archivo BigWig:")
# print("Chromosomas y sus longitudes:")
# for chrom, length in bw.chroms().items():
#     print(f"{chrom}: {length}")

# print("\nIntervalos en los primeros 100,000 bp de cada cromosoma:")
# for chrom in bw.chroms().keys():
#     intervals = bw.intervals(chrom)
#     print(f"{chrom}: {intervals}")

# # Cerrar el archivo BigWig
# bw.close()

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


data = Pyntegrate.results_table.ResultsTable('test_RNA-seq_w_changes.csv')
###Primero eliminamos aquellos genes que no existan en rna pero si en atac
print(data)
print(type(data.data['GeneID']))
print("ENSMUSG00000035179" in data.data['GeneID'].values)

db = gffutils.FeatureDB(os.path.join(path, 'data/gencode.vM25.annotation.gtf.db'))
print(db)
atac = Pyntegrate.genomic_signal('test_2_DNase-seq.bw', 'bigwig')
tsses = pybedtools.BedTool(tss_generator()).saveas('tsses.gtf')


tsses_1kb = tsses.slop(b=1000, genome='mm10', output='tsses-1kb.gtf')
processes = multiprocessing.cpu_count()


##Cargamos chip-seq, lo procesamos y enseñamos algunas gráficas
atac_signal =  atac.array(
    tsses_1kb,
    bins=100,
    processes=processes
)

# print(atac_signal)
#Quitando lso q no tienen picos
#########################################
atac_signal_semigood = []
for x in atac_signal:
    if np.mean(x[0]['values']) != 0:
        x[0]['id'] = x[0]['id'].split('.')[0]
        atac_signal_semigood.append(x)
#####################################3



###AHora eliminamos aquiellos genes que no existen en rna pero si en atac
atac_signal_good = []
ids_selected = []
for x in atac_signal_semigood:
    # print(x)
    if x[0]['id'] in  data.data['GeneID'].values:
        if x[0]['id'] not in ids_selected:
            atac_signal_good.append(x)
            ids_selected.append(x[0]['id'])


print(len(atac_signal_good))



#########################################

###Ahora eliminas aquelos genes que no existen en atac pero si en RNA
indexes = []
i = 0
for x in data.data['GeneID'].values:
    if x in ids_selected:
        indexes.append(i)
    i += 1
    
data_filtered = data.iloc[indexes]

data = data_filtered
print("Data lenght after filtered:",len(data))

##############################################################



normalized_subtracted = Pyntegrate.chipSeqSignalAnalysis.value_array_simple(atac_signal_good)
normalized_subtracted_2 = np.array(normalized_subtracted)
normalized_subtracted_2 /= atac.mapped_read_count() / 1e6
print(len(normalized_subtracted_2))

x = np.linspace(-1000, 1000, bins)


fig = Pyntegrate.plotutils.imshow(

    # The array to plot; here, we've subtracted input from IP.
    normalized_subtracted_2,

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
    normalized_subtracted_2,
    x=x,
    figsize=(3, 7),
    vmin=5, vmax=99,  percentile=True,
    line_kwargs=dict(color='k', label='All'),
    fill_kwargs=dict(color='k', alpha=0.3),

    # This is new: sort by mean signal
    sort_by=normalized_subtracted.mean(axis=1)
)

fig = Pyntegrate.plotutils.imshow(

    # These are the same arguments as above.
    normalized_subtracted_2,
    x=x,
    figsize=(3, 7),
    vmin=5, vmax=99,  percentile=True,
    line_kwargs=dict(color='k', label='All'),
    fill_kwargs=dict(color='k', alpha=0.3),

    # This is new: sort by mean signal
    sort_by=np.argmax(normalized_subtracted, axis=1)
)



fig = Pyntegrate.plotutils.imshow(
    normalized_subtracted_2,
    x=x,
    figsize=(3, 7),
    vmin=5, vmax=99,  percentile=True,
    line_kwargs=dict(color='k', label='All'),
    fill_kwargs=dict(color='k', alpha=0.3),
    sort_by=normalized_subtracted.mean(axis=1)
)

# plt.show()

print(len(normalized_subtracted_2))
# print(atac_signal_good)



fig = Pyntegrate.plotutils.imshow(
    # Same as before...
    normalized_subtracted_2,
    x=x,
    figsize=(3, 7),
    vmin=5, vmax=99,  percentile=True,
    line_kwargs=dict(color='k', label='All'),
    fill_kwargs=dict(color='k', alpha=0.3),
    sort_by=normalized_subtracted_2.mean(axis=1),
    
    
    # Default was (3,1); here we add another number 
    height_ratios=(3, 1, 1)
)

# `fig.gs` contains the `matplotlib.gridspec.GridSpec` object,
# so we can now create the new axes.
bottom_axes = plt.subplot(fig.gs[2, 0])


# Signal over TSSs of transcripts that were activated upon knockdown.
Pyntegrate.plotutils.ci_plot(
    x,
    normalized_subtracted_2[(data.log2foldchange > 1).values, :],
    line_kwargs=dict(color='#fe9829', label='up'),
    fill_kwargs=dict(color='#fe9829', alpha=0.3),
    ax=bottom_axes)

# Signal over TSSs of transcripts that were repressed upon knockdown
Pyntegrate.plotutils.ci_plot(
    x,
    normalized_subtracted_2[(data.log2foldchange < -1).values, :],
    line_kwargs=dict(color='#8e3104', label='down'),
    fill_kwargs=dict(color='#8e3104', alpha=0.3),
    ax=bottom_axes)

# Signal over TSSs tof transcripts that did not change upon knockdown
Pyntegrate.plotutils.ci_plot(
    x,
    normalized_subtracted_2[((data.log2foldchange >= -1) & (data.log2foldchange <= 1)).values, :],
    line_kwargs=dict(color='.5', label='unchanged'),
    fill_kwargs=dict(color='.5', alpha=0.3),
    ax=bottom_axes)

# Clean up redundant x tick labels, and add axes labels
fig.line_axes.set_xticklabels([])
fig.array_axes.set_xticklabels([])
fig.line_axes.set_ylabel('Average\nenrichement')
fig.array_axes.set_ylabel('Transcripts on chr17')
bottom_axes.set_ylabel('Average\nenrichment')
bottom_axes.set_xlabel('Distance from TSS (bp)')
fig.cax.set_ylabel('Enrichment')

# Add the vertical lines for TSS position to all axes
for ax in [fig.line_axes, fig.array_axes, bottom_axes]:
    ax.axvline(0, linestyle=':', color='k')

# Nice legend
bottom_axes.legend(loc='best', frameon=False, fontsize=8, labelspacing=.3, handletextpad=0.2)
fig.subplots_adjust(left=0.3, right=0.8, bottom=0.05)


plt.show()
