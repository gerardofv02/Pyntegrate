import Pyntegrate
import os
import gffutils
import pybedtools
from pybedtools.featurefuncs import TSS
from gffutils.helpers import asinterval
import multiprocessing


def tss_generator():
    """
    Generator function to yield TSS of each annotated transcript
    """
    for transcript in db.features_of_type('transcript'):
        yield TSS(asinterval(transcript), upstream=1, downstream=0)

path = os.path.dirname(os.path.abspath(__file__))
print(path)

features, arrays = Pyntegrate.chipSeqSignalAnalysis.generate_arrays_features_from_tsses_from_db(os.path.join(path, 'data/gencode.vM25.annotation.gtf.db'), os.path.join(path, 'bwFiles/SRR1204546_sort_nondup.bw'),os.path.join(path, 'bwFiles/SRR1204544_sort_nondup.bw'),genome = "mm10")

# # # db = gffutils.create_db(data=path+"/Homo_sapiens.GRCh38.109.gtf",dbfn=path+"/data/Homo_sapiens.GRCh38.109.gtf.db")
# # db = gffutils.create_db(data=path+"/prueba.gtf",dbfn=path+"/data/prueba.gtf.db")

# db = gffutils.FeatureDB(os.path.join(path, 'data/gencode.vM25.annotation.gtf.db'))
# print(db)

# # fn = "gencode.vM25.annotation.gtf"

# # db = open(fn).read()






# tsses = pybedtools.BedTool(tss_generator()).saveas('tsses.gtf')


# tsses_1kb = tsses.slop(b=1000, genome='mm10', output='tsses-1kb.gtf')






# ip_signal = Pyntegrate.genomic_signal(
#     os.path.join(path, 'bwFiles/SRR1204546_sort_nondup.bw'),
#     'bigwig')


# input_signal = Pyntegrate.genomic_signal(
#     os.path.join(path, 'bwFiles/SRR1204544_sort_nondup.bw'),
#     'bigwig')

# print(ip_signal,input_signal)
# processes = multiprocessing.cpu_count()

# if not os.path.exists('example.npz'):

#     # The signal is the IP ChIP-seq BAM file.
#     ip_array = ip_signal.array(

#         # Look at signal over these windows
#         tsses_1kb,

#         # Bin signal into this many bins per window
#         bins=100,

#         # Use multiple CPUs. Dramatically speeds up run time.
#         processes=processes)

#     print(ip_array[0][:10])
#     # Do the same thing for input.
#     input_array = input_signal.array(
#         tsses_1kb,
#         bins=100,
#         processes=processes)

#     print(input_array[:10])



#     # Normalize to library size. The values in the array
#     # will be in units of "reads per million mapped reads"
#     ip_array /= ip_signal.mapped_read_count() / 1e6
#     input_array /= input_signal.mapped_read_count() / 1e6


#     # Cache to disk. The data will be saved as "example.npz" and "example.features".
#     Pyntegrate.persistence.save_features_and_arrays(
#         features=tsses,
#         arrays={'ip': ip_array, 'input': input_array},
#         prefix='example',
#         link_features=True,
#         overwrite=True)

# features, arrays = Pyntegrate.persistence.load_features_and_arrays(prefix='example')




import numpy as np
x = np.linspace(-1000, 1000, 100)

# Import plotting tools
from matplotlib import pyplot as plt


# Create a figure and axes
fig = plt.figure()
ax = fig.add_subplot(111)


# Plot the IP:
ax.plot(
    # use the x-axis values we created
    x,

    # axis=0 takes the column-wise mean, so with
    # 100 columns we'll have 100 means to plot
    arrays['ip'].mean(axis=0),

    # Make it red
    color='r',

    # Label to show up in legend
    label='IP')


# Do the same thing with the input
ax.plot(
    x,
    arrays['input'].mean(axis=0),
    color='k',
    label='input')


# Add a vertical line at the TSS, at position 0
ax.axvline(0, linestyle=':', color='k')


# Add labels and legend
ax.set_xlabel('Distance from TSS (bp)')
ax.set_ylabel('Average read coverage (per million mapped reads)')
ax.legend(loc='best')

plt.show()


# Calcular la diferencia después de ajustar los tamaños
print(arrays["ip"])

normalized_subtracted = arrays["ip"] - arrays["input"]
# Tweak some font settings so the results look nicer
# plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10


# the metaseq.plotutils.imshow function does a lot of work,
# we just have to give it the right arguments:
fig = Pyntegrate.plotutils.imshow(

    # The array to plot; here, we've subtracted input from IP.
    normalized_subtracted,

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
    normalized_subtracted,
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
    normalized_subtracted,
    x=x,
    figsize=(3, 7),
    vmin=5, vmax=99,  percentile=True,
    line_kwargs=dict(color='k', label='All'),
    fill_kwargs=dict(color='k', alpha=0.3),

    # This is new: sort by mean signal
    sort_by=np.argmax(normalized_subtracted, axis=1)
)



fig = Pyntegrate.plotutils.imshow(
    normalized_subtracted,
    x=x,
    figsize=(3, 7),
    vmin=5, vmax=99,  percentile=True,
    line_kwargs=dict(color='k', label='All'),
    fill_kwargs=dict(color='k', alpha=0.3),
    sort_by=normalized_subtracted.mean(axis=1)
)

plt.show()
# "line_axes" is our handle for working on the lower axes.
# Add some nicer labels.
fig.line_axes.set_ylabel('Average enrichment')
fig.line_axes.set_xlabel('Distance from TSS (bp)')

# "array_axes" is our handle for working on the upper array axes.
# Add a nicer axis label
fig.array_axes.set_ylabel('Transcripts on chr17')

# Remove the x tick labels, since they're redundant
# with the line axes
fig.array_axes.set_xticklabels([])

# Add a vertical line to indicate zero in both the array axes
# and the line axes
fig.array_axes.axvline(0, linestyle=':', color='k')
fig.line_axes.axvline(0, linestyle=':', color='k')

fig.cax.set_ylabel("Enrichment")

print(tsses_1kb[0])
data = Pyntegrate.results_table.DEseq2ResultsPrueba("./desq2_Sanchez.csv")
print(data)

fig = Pyntegrate.plotutils.imshow(
    # Same as before...
    normalized_subtracted,
    x=x,
    figsize=(3, 7),
    vmin=5, vmax=99,  percentile=True,
    line_kwargs=dict(color='k', label='All'),
    fill_kwargs=dict(color='k', alpha=0.3),
    sort_by=normalized_subtracted.mean(axis=1),


    # Default was (3,1); here we add another number
    height_ratios=(3, 1, 1)
)

# `fig.gs` contains the `matplotlib.gridspec.GridSpec` object,
# so we can now create the new axes.

bottom_axes = plt.subplot(fig.gs[2, 0])
upregulated = (data.log2FoldChange > 1)
index = upregulated.values
# upregulated_chipseq_signal = normalized_subtracted[index, :]
# subset = normalized_subtracted[(data.log2FoldChange > 1).values, :]

##BAM FILES
# bottom_axes = plt.subplot(fig.gs[2, 0])
# print(bottom_axes)

# upregulated = (data.log2FoldChange > 1)
# print(upregulated)

# index =np.where(upregulated)[0]
# print(index)
# upregulated_chipseq_signal = np.take(normalized_subtracted, index, axis=0)
# print(upregulated_chipseq_signal)
# valid_indices = ~data.log2FoldChange.isna() & ((data.log2FoldChange > 1).values)
# valid_indices_minus = ~data.log2FoldChange.isna() & ((data.log2FoldChange < -1).values)
# valid_indices_minus_mas = ~data.log2FoldChange.isna() & ((data.log2FoldChange >= -1).values) &((data.log2FoldChange <= 1).values)
# subset = np.take(normalized_subtracted, valid_indices, axis=0)
# print(subset)
Pyntegrate.plotutils.ci_plot(
    x,
    normalized_subtracted[(data.log2FoldChange > 1).values, :],
    line_kwargs=dict(color='#fe9829', label='up'),
    fill_kwargs=dict(color='#fe9829', alpha=0.3),
    ax=bottom_axes)

# Signal over TSSs of transcripts that were repressed upon knockdown
Pyntegrate.plotutils.ci_plot(
    x,
    normalized_subtracted[(data.log2FoldChange < -1).values, :],
    line_kwargs=dict(color='#8e3104', label='down'),
    fill_kwargs=dict(color='#8e3104', alpha=0.3),
    ax=bottom_axes)

# Signal over TSSs tof transcripts that did not change upon knockdown
Pyntegrate.plotutils.ci_plot(
    x,
    normalized_subtracted[((data.log2FoldChange >= -1) & (data.log2FoldChange <= 1)).values, :],
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
# plt.show()