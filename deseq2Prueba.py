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

features, arrays ,tsses , tsses_1kb = Pyntegrate.chipSeqSignalAnalysis.generate_arrays_features_from_tsses_from_db(os.path.join(path, 'data/gencode.vM25.annotation.gtf.db'), os.path.join(path, 'BamFiles/SRR1204544_sort_nondup.bw'),'bigwig',os.path.join(path, 'BamFiles/SRR1204546_sort_nondup.bw'),'bigwig',genome = "mm10")

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


#ESto hay que verlo si ponerlo dentro o fuera de las funciones. Por ahora va dentro
#########################################
# import numpy as np
# x = np.linspace(-1000, 1000, 100)
#########################################

# Import plotting tools

#Haciendo gráfica de ip e input distancia del tss
##################################################################
x = np.linspace(-1000, 1000, 100)


# # Create a figure and axes
# fig = plt.figure()
# ax = fig.add_subplot(111)


# # Plot the IP:
# ax.plot(
#     # use the x-axis values we created
#     x,

#     # axis=0 takes the column-wise mean, so with
#     # 100 columns we'll have 100 means to plot
#     arrays['ip'].mean(axis=0),

#     # Make it red
#     color='r',

#     # Label to show up in legend
#     label='IP')


# # Do the same thing with the input
# ax.plot(
#     x,
#     arrays['input'].mean(axis=0),
#     color='k',
#     label='input')


# # Add a vertical line at the TSS, at position 0
# ax.axvline(0, linestyle=':', color='k')


# # Add labels and legend
# ax.set_xlabel('Distance from TSS (bp)')
# ax.set_ylabel('Average read coverage (per million mapped reads)')
# ax.legend(loc='best')

fig = Pyntegrate.chipSeqSignalAnalysis.distance_from_tss_chipSeq(arrays_ip = arrays["ip"], arrays_input=arrays["input"])
plt.show()
##################################################################



#Calculando picos
################################
# normalized_subtracted = arrays["ip"] - arrays["input"]
normalized_subtracted= Pyntegrate.chipSeqSignalAnalysis.calculate_peaks(arrays_ip=arrays["ip"],arrays_input= arrays["input"])
################################
print("Type of normaliz subtracted before all", type(normalized_subtracted))
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


# "line_axes" is our handle for working on the lower axes.
# Add some nicer labels.
fig.line_axes.set_ylabel('Average enrichment')
fig.line_axes.set_xlabel('Distance from TSS (bp)')

# "array_axes" is our handle for working on the upper array axes.
# Add a nicer axis label
fig.array_axes.set_ylabel('Transcripts')

# Remove the x tick labels, since they're redundant
# with the line axes
fig.array_axes.set_xticklabels([])

# Add a vertical line to indicate zero in both the array axes
# and the line axes
fig.array_axes.axvline(0, linestyle=':', color='k')
fig.line_axes.axvline(0, linestyle=':', color='k')

fig.cax.set_ylabel("Enrichment")




data = Pyntegrate.results_table.DEseq2ResultsPrueba("./desq2_Sanchez.csv",  import_kwargs=dict(index_col=1))
print("\ndata: ",data,"len data: ", len(data) )
data.reindex_to(tsses, attribute='gene_name')
print("\ndata: ",data,"len data: ", len(data) )
# data2 = data.data
import pandas as pd
df2 = pd.read_csv("./tsses.gtf", sep="\t" ,header= None)
print(df2[8].values)
gene_used = []
indexs = []
values_not_data = []
for idx,value in enumerate(df2.values):
    if value[8].split(";")[3].split(" ")[2].split('"')[1] in data.index:
        indexs.append(idx)
        gene_used.append(value[8].split(";")[3].split(" ")[2].split('"')[1])


print(len(indexs))
print("\n Gene_used: ", gene_used)

for value_data in data.index:
    if value_data not in gene_used:
        values_not_data.append(value_data)

print("Len of values not data",len(values_not_data))
normalized_subtracted_good= []
##For put good normalized subtracted
for i in range(len(normalized_subtracted)):
    if i in indexs:
        normalized_subtracted_good.append(normalized_subtracted[i])
##FOr deseq2
data = data.drop(values_not_data)

normalized_subtracted= np.array(normalized_subtracted_good) 
print("Normalized subtract len: ", len(normalized_subtracted))

# data2['log2FoldChange'] = data["fpkm"]
# print("Len data2: ",len(data2), "DAta2: ",data2)
# print("\n\nlen data: ",len(data), "Data: ", data)

fig = Pyntegrate.plotutils.imshow(
    # Same as before...
    normalized_subtracted,
    x=x,
    figsize=(4, 7),
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

print("\n\n\n\n\n, Bottom axes: ", bottom_axes)
upregulated = (data.log2FoldChange > 1)
print("\n\n\n\nUpregulated: ", upregulated, "\nLen upregulated:  ", len(upregulated))
index = upregulated.values
print("\n\n\n Index: ", index, "\nlen de index:", len(index))
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
# # # # # normalized_subtracted = normalized_subtracted[:14543]
print("\n\nnormalized, datafoldchange > 1", normalized_subtracted[(data.log2FoldChange > 1).values, :], "Len: ", len(normalized_subtracted[(data.log2FoldChange > 1).values, :]))

# Signal over TSSs of transcripts that were activated upon knockdown.
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
bottom_axes.legend(loc=0,bbox_to_anchor=(1, 1), frameon=False, fontsize=8, labelspacing=.3, handletextpad=0.2)##aki
fig.subplots_adjust(left=0.3, right=0.8, bottom=0.05)
# plt.show()
plt.show()