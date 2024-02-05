data_dir = 'metaseq-example/data'

import os
import gffutils
import pybedtools
from pybedtools.featurefuncs import TSS
from gffutils.helpers import asinterval
import Pyntegrator
import multiprocessing
import numpy as np
from matplotlib import pyplot as plt
from Pyntegrator.results_table import ResultsTable


db = gffutils.FeatureDB(os.path.join(data_dir, 'Homo_sapiens.GRCh37.66_chr17.gtf.db'))


def tss_generator():
    """
    Generator function to yield TSS of each annotated transcript
    """
    for transcript in db.features_of_type('transcript'):
        yield TSS(asinterval(transcript), upstream=1, downstream=0)


# A BedTool made out of a generator, and saved to file.
tsses = pybedtools.BedTool(tss_generator()).saveas('tsses.gtf')
tsses_1kb = tsses.slop(b=1000, genome='hg19', output='tsses-1kb.gtf')

ip_signal = Pyntegrator.genomic_signal(
    os.path.join(data_dir, 'wgEncodeHaibTfbsK562Atf3V0416101AlnRep1_chr17.bam'),
    'bam')

print("\n\n\n\n\n\n\n\n\n",ip_signal)

input_signal = Pyntegrator.genomic_signal(
    os.path.join(data_dir, 'wgEncodeHaibTfbsK562RxlchV0416101AlnRep1_chr17.bam'),
    'bam')


processes = multiprocessing.cpu_count()

if not os.path.exists('example.npz'):

    # The signal is the IP ChIP-seq BAM file.
    ip_array = ip_signal.array(

        # Look at signal over these windows
        tsses_1kb,

        # Bin signal into this many bins per window
        bins=100,

        # Use multiple CPUs. Dramatically speeds up run time.
        processes=processes)
    print("\n\n\n\n\n\n\n", ip_array)

    # Do the same thing for input.
    input_array = input_signal.array(
        tsses_1kb,
        bins=100,
        processes=processes)

    print("\n\n\n\n\n\n\n\n", input_array)

    # Normalize to library size. The values in the array
    # will be in units of "reads per million mapped reads"

    print(ip_signal.mapped_read_count())
    
    mapped_count = ip_signal.mapped_read_count()

    ip_array /= mapped_count / 1e6



    input_array /= input_signal.mapped_read_count() / 1e6

    # Cache to disk. The data will be saved as "example.npz" and "example.features".
    Pyntegrator.persistence.save_features_and_arrays(
        features=tsses,
        arrays={'ip': ip_array, 'input': input_array},
        prefix='example',
        link_features=True,
        overwrite=True)

features, arrays = Pyntegrator.persistence.load_features_and_arrays(prefix='example')

# How many features?
assert len(features) == 5708

# This ought to be exactly the same as the number of features in `tsses_1kb.gtf`
assert len(features) == len(tsses_1kb) == 5708

# This shows that `arrays` acts like a dictionary
assert sorted(arrays.keys()) == ['input', 'ip']

# This shows that the IP and input arrays have one row per feature, and one column per bin
assert arrays['ip'].shape == (5708, 100) == arrays['input'].shape

x = np.linspace(-1000, 1000, 100)

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

#plt.show()

normalized_subtracted = arrays['ip'] - arrays['input']

# Tweak some font settings so the results look nicer
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10

# the Pyntegrator.plotutils.imshow function does a lot of work,
# we just have to give it the right arguments:
fig = Pyntegrator.plotutils.imshow(

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

#plt.show()


fig = Pyntegrator.plotutils.imshow(

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

#plt.show()

fig = Pyntegrator.plotutils.imshow(

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

#plt.show()

fig = Pyntegrator.plotutils.imshow(
    normalized_subtracted,
    x=x,
    figsize=(3, 7),
    vmin=5, vmax=99,  percentile=True,
    line_kwargs=dict(color='k', label='All'),
    fill_kwargs=dict(color='k', alpha=0.3),
    sort_by=normalized_subtracted.mean(axis=1)
)

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

#plt.show()

control = ResultsTable(
    os.path.join(data_dir, 'GSM847565_SL2585.table'),
    import_kwargs=dict(index_col=0))

knockdown = ResultsTable(
    os.path.join(data_dir, 'GSM847566_SL2592.table'),
    import_kwargs=dict(index_col=0))


print(len(control.data))
print(control.data.head())



# ---------------------------------------------------------
# Re-align the ResultsTables to match the GTF file
control = control.reindex_to(tsses, attribute='transcript_id')
knockdown = knockdown.reindex_to(tsses, attribute='transcript_id')

print (len(control))
control.data.head()


original_control = ResultsTable(
    os.path.join(data_dir, 'GSM847565_SL2585.table'),
    import_kwargs=dict(index_col=0))

print('ENST00000419929' in original_control.data.index)

print(len(control.data),len(knockdown.data),len(tsses_1kb))

# Everything should be the same length
assert len(control.data) == len(knockdown.data) == len(tsses_1kb) == 5708

# Spot-check some values to make sure the GTF file and the DataFrame match up.
assert tsses[0]['transcript_id'] == control.data.index[0]
assert tsses[100]['transcript_id'] == control.data.index[100]
assert tsses[5000]['transcript_id'] == control.data.index[5000]


# Join the dataframes and create a new pandas.DataFrame.
data = control.data.join(knockdown.data, lsuffix='_control', rsuffix='_knockdown')

# Add a log2 fold change variable
data['log2foldchange'] = np.log2(data.fpkm_knockdown / data.fpkm_control)
print(data.head())

# ---------------------------------------------------------
# How many transcripts on chr17 changed expression?

print( "up:", sum(data.log2foldchange > 1))
print( "down:", sum(data.log2foldchange < -1))



fig = Pyntegrator.plotutils.imshow(
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

# plt.show()

# This is a pandas.Series, True where the log2foldchange was >1
upregulated = (data.log2foldchange > 1)
print(upregulated)

# This gets us the underlying boolean NumPy array which we
# can use to subset the array
index = upregulated.values
print(index)

# This is the subset of the array where the TSS of the transcript
# went up in the ATF3 knockdown
upregulated_chipseq_signal = normalized_subtracted[index, :]
print(upregulated_chipseq_signal)
# We can combine the above steps into the following:
subset = normalized_subtracted[(data.log2foldchange > 1).values, :]

# Signal over TSSs of transcripts that were activated upon knockdown.
Pyntegrator.plotutils.ci_plot(
    x,
    normalized_subtracted[(data.log2foldchange > 1).values, :],
    line_kwargs=dict(color='#fe9829', label='up'),
    fill_kwargs=dict(color='#fe9829', alpha=0.3),
    ax=bottom_axes)

# Signal over TSSs of transcripts that were repressed upon knockdown
Pyntegrator.plotutils.ci_plot(
    x,
    normalized_subtracted[(data.log2foldchange < -1).values, :],
    line_kwargs=dict(color='#8e3104', label='down'),
    fill_kwargs=dict(color='#8e3104', alpha=0.3),
    ax=bottom_axes)

# Signal over TSSs tof transcripts that did not change upon knockdown
Pyntegrator.plotutils.ci_plot(
    x,
    normalized_subtracted[((data.log2foldchange >= -1) & (data.log2foldchange <= 1)).values, :],
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

# fig.savefig('demo.png')
# fig.savefig('demo.svg')



# K-means input data should be normalized (mean=0, stddev=1)
from sklearn import preprocessing
X_scaled = preprocessing.scale(normalized_subtracted)

k = 4

ind, breaks = Pyntegrator.plotutils.new_clustered_sortind(

    # The array to cluster
    X_scaled,

    # Within each cluster, how the rows should be sorted
    row_key=np.mean,

    # How each cluster should be sorted
    cluster_key=np.median,

    # Number of clusters
    k=k)

# Plot the heatmap again
fig = Pyntegrator.plotutils.imshow(
    normalized_subtracted,
    x=x,
    figsize=(3, 9),
    vmin=5, vmax=99,  percentile=True,
    line_kwargs=dict(color='k', label='All'),
    fill_kwargs=dict(color='k', alpha=0.3),

    # A little tricky: `sort_by` expects values to sort by
    # (say, expression values). But we've pre-calculated
    # our actual sort index based on clusters, so we transform
    # it like this
    sort_by=np.argsort(ind),

    # This adds a "strip" axes on the right side, useful
    # for adding extra information. We'll add cluster color
    # codes here.
    strip=True,
)

# De-clutter by hiding labels
plt.setp(
    fig.strip_axes.get_yticklabels()
    + fig.strip_axes.get_xticklabels()
    + fig.array_axes.get_xticklabels(),
    visible=False)

#
fig.line_axes.set_ylabel('Average\nenrichement')
fig.array_axes.set_ylabel('Transcripts on chr17')
fig.strip_axes.yaxis.set_label_position('right')
fig.strip_axes.set_ylabel('Cluster')
fig.cax.set_ylabel('Enrichment')

# Make colors
import matplotlib
cmap = matplotlib.cm.Spectral
colors = cmap(np.arange(k) / float(k))

# This figure will contain average signal for each cluster
fig2 = plt.figure(figsize=(10,3))


last_break = 0
cluster_number = 1
n_panel_rows = 1
n_panel_cols = k
for color, this_break in zip(colors, breaks):
    if cluster_number == 1:
        sharex = None
        sharey = None
    else:
        sharex = fig2.axes[0]
        sharey = fig2.axes[0]

    ax = fig2.add_subplot(
        n_panel_rows,
        n_panel_cols,
        cluster_number,
        sharex=sharex,
        sharey=sharey)

    # The y position is somewhat tricky: the array was
    # displayed using matplotlib.imshow with the argument
    # `origin="lower"`, which means the row in the plot at y=0
    # corresponds to the last row in the array (index=-1).
    # But the  breaks are in array coordinates. So we convert
    # them by subtracting from the total array size.
    xpos = 0
    width = 1
    ypos = len(normalized_subtracted) - this_break
    height = this_break - last_break
    rect = matplotlib.patches.Rectangle(
        (xpos, ypos), width=width, height=height, color=color)
    fig.strip_axes.add_patch(rect)
    fig.array_axes.axhline(ypos, color=color, linewidth=2)

    chunk = normalized_subtracted[last_break:this_break]

    Pyntegrator.plotutils.ci_plot(
        x,
        chunk,
        ax=ax,
        line_kwargs=dict(color=color),
        fill_kwargs=dict(color=color, alpha=0.3),
        )
    ax.axvline(0, color='k', linestyle=':')
    ax.set_title('cluster %s\n(N=%s)' % (cluster_number, len(chunk)))
    cluster_number += 1
    last_break = this_break



# Convert to ResultsTable so we can take advantage of its
# `scatter` method
rt = ResultsTable(data)

# Get the up/down regulated
up = rt.log2foldchange > 1
dn = rt.log2foldchange < -1

# Go back to the ChIP-seq data and create a boolean array
# that is True only for the top TSSes with the strongest
# mean signal
tss_means = normalized_subtracted.mean(axis=1)
strongest_signal = np.zeros(len(tss_means)) == 1
strongest_signal[np.argsort(tss_means)[-25:]] = True

rt.scatter(
    x='fpkm_control',
    y='fpkm_knockdown',
    xfunc=np.log1p,
    yfunc=np.log1p,
    genes_to_highlight=[
        (up, dict(color='#da3b3a', alpha=0.8)),
        (dn, dict(color='#00748e', alpha=0.8)),
        (strongest_signal, dict(color='k', s=50, alpha=1)),

    ],
    general_kwargs=dict(marker='.', color='0.5', alpha=0.2, s=5),
    one_to_one=dict(color='r', linestyle=':')
)


# Perhaps a better analysis would be to plot average
# ChIP-seq signal vs log2foldchange directly. In an imaginary
# world where biology is simple, we might expect TSSes with stronger
# log2foldchange upon knockdown to have stronger ChIP-seq signal
# in the control.
#
# To take advantage of the `scatter` method of ResultsTable objects,
# we simply add the TSS signal means as another variable in the
# dataframe. Then we can refer to it by name in `scatter`.
#
# We'll also use the same colors and genes to highlight from
# above.

rt.data['tss_means'] = tss_means
rt.scatter(
    x='log2foldchange',
    y='tss_means',
    genes_to_highlight=[
        (up, dict(color='#da3b3a', alpha=0.8)),
        (dn, dict(color='#00748e', alpha=0.8)),
        (strongest_signal, dict(color='k', s=50, alpha=1)),
    ],
    general_kwargs=dict(marker='.', color='0.5', alpha=0.2, s=5),
    yfunc=np.log2)
plt.show()