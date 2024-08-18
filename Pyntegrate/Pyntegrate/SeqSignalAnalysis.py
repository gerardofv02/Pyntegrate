from ._genomic_signal import *
from .persistence import *
from .array_helpers import *
from .plotutils import *
from .results_table import *
import os
import gffutils
import pybedtools
from pybedtools.featurefuncs import TSS
from gffutils.helpers import asinterval
import multiprocessing
import numpy as np
import sys
from matplotlib import pyplot as plt
import subprocess
from .plotutils import *

def tss_generator(db):
    """
    Generator function to yield TSS of each annotated transcript
    """
    for transcript in db.features_of_type('transcript'):
        yield TSS(asinterval(transcript), upstream=1, downstream=0)

def generate_arrays_features_from_tsses_from_db(dbPath, ipSignalPath, extensionIp,inputSignalPath,extensionInput, genome, bins):

    """
    Function to process ChIP-seq data, specifically to analyze the protein enrichment signal in certain genomic regions, such as transcription start sites (TSS)

    Params:
        - dbPath: String: The path where database of CHIP-seq data is stored (It must be a db generated by the CHIP-seq data with pybedtools library)
        - ipSignalPath: String: The path where IP signal data is stored
        - extensionIp: String: Type of the file of the IP signal data file. Example:('bigwig','bam','bed',...)
        - inputSignalPath: String: The path where INPUT signal data is stored
        - extensionInput: String: Type of the file of the INPUT signal data file example:('bigwig','bam','bed',...)
        - genome: String: The genome of the CHIP-seq data we want to analyze. Example: ('hg19','hg38',...)
        - bins:Int: The amount of subintervals you want to have on each interval
    """

    db = gffutils.FeatureDB(dbPath)
    tsses = pybedtools.BedTool(tss_generator(db)).saveas('tsses.gtf')

    remove_duplicates("tsses.gtf")
    tsses = pybedtools.BedTool("tsses.gtf")
    tsses_1kb = tsses.slop(b=1000, genome=genome, output='tsses-1kb.gtf')

    ip_signal = genomic_signal(ipSignalPath,extensionIp)
    input_signal = genomic_signal(inputSignalPath,extensionInput)




    print(ip_signal,input_signal)
    processes = multiprocessing.cpu_count()
    if not os.path.exists('example.npz'):

        ip_array = ip_signal.array(
            tsses_1kb,
            bins=bins,
            processes=processes)
        input_array = input_signal.array(
            tsses_1kb,
            bins=bins,
            processes=processes)
        

        ip_array_good = []
        input_array_good = []
        for i in ip_array:
            i[0]['values'] /= ip_signal.mapped_read_count() / 1e6
            ip_array_good.append(i[0])

        for y in input_array:
            y[0]['values'] /= input_signal.mapped_read_count() / 1e6
            input_array_good.append(y[0])
        save_features_and_arrays(
            features=tsses,
            arrays={'ip': ip_array_good, 'input': input_array_good},
            prefix='example',
            link_features=True,
            overwrite=True)
        
    features, arrays = load_features_and_arrays(prefix='example')
    return features ,arrays ,tsses, tsses_1kb


def generate_array_simple_signal(dbPath,filePath,extensionFile,genome,bins):

    """
    Function to get featurees, tsses, tsses_1kb and the array of values with gene_names of a simple signal

    Params: 
        - dbPath: String: The path where database of CHIP-seq data is stored (It must be a db generated by the CHIP-seq data with pybedtools library)
        - filePath: String: The file path where is the signal data stored
        - extensionFile: String: Type of the file of the signal data file example:('bigwig','bam','bed',...)
        - genome: String: The genome of the CHIP-seq data we want to analyze. Example: ('hg19','hg38',...)
        - bins: Int: The amount of subintervals you want to have on each interval
    """

    db = gffutils.FeatureDB(dbPath)
    tsses = pybedtools.BedTool(tss_generator(db)).saveas('tsses.gtf')

    remove_duplicates("tsses.gtf")
    tsses = pybedtools.BedTool("tsses.gtf")
    tsses_1kb = tsses.slop(b=1000, genome=genome, output='tsses-1kb.gtf')

    signal = genomic_signal(filePath,extensionFile)
    processes = multiprocessing.cpu_count()
    array = signal.array(
    tsses_1kb,
    bins=bins,
    processes=processes)

    array_good = []
    for i in array:
        i[0]['values'] /= signal.mapped_read_count() / 1e6
        array_good.append(i[0])
    features = tsses
    return features ,array_good ,tsses, tsses_1kb
    


def value_array_simple(array):
    """
    This function convert a simple array of object with values and gene info into an array of only values.
    THis works when an open file with IP or INPUT data is open, we can arrify it with the array function of
    the object and it returns an array of object with the values an gene information. So this funcion helps 
    to just have an array of values.

    Params:
        - array: Array get by the array function of BaseSignal object
    """
    array_final = []
    for x in array:
    #    print(x, type(x))
       array_final.append(x['values'])
    array_final = np.array(array_final)

    return array_final

def values_array(array_ip, array_input):

    """
    This function makes the same as the one before but with two arrays because normally it is used to IP and INPUT

    Params:    
        - arrays_ip: IP array get by the array function of the BaseSignal object
        - arrays_input: INPUT array get by array function of the BaseSignal object
    """
    arrays_ip = []
    arrays_input = []
    for x in array_ip:
       arrays_ip.append( x['values'])
    for y in array_input:
       arrays_input.append( y['values'])

    arrays_ip = np.array(arrays_ip)
    arrays_input = np.array(arrays_input)

    return arrays_ip, arrays_input


def calculate_peaks_with_gene_name(arrays_ip, arrays_input):

    """
    Function to calculate the peaks between IP signal and INPUT signal having an array of objects.

    Params:
        - arrays_ip: Array: The array of the IP signal data (dict object with more information as gene_name)
        - arrays_input: Array: The array of the INPUT signal data (dict object with more information as gene_name)
    """

    normalized_subtracted = np.array([])
    try_boolean = []
    for index,x in enumerate(arrays_ip):

        normalized_subtracted = np.append(normalized_subtracted,{ 'values': arrays_ip[index]['values']-arrays_input[index]['values'],
                                                                  'gene_name': arrays_ip[index]['gene_name'],
                                                                  'chr': arrays_ip[index]['chr'],
                                                                  'start': arrays_ip[index]['start'],
                                                                  'end': arrays_ip[index]['end'],
                                                                  'strand': arrays_ip[index]['strand'],
                                                                  'id': arrays_ip[index]['id']})


    return normalized_subtracted


def calculate_peaks(arrays_ip, arrays_input):

    """
    Function to calculate the peaks between the IP signal and INPUT signal.
    Params:
        - arrays_ip: Array: The array of the IP signal data (just values)
        - arrays_input: Array: The array of the INPUT signal data (just values)
    """

    if(len(arrays_ip) != len(arrays_input)):
        sys.stderr.write("Length of ip is different from input array")
    peaks = arrays_ip - arrays_input
    return peaks


def distance_from_tss_chipSeq(arrays_ip, arrays_input,xAxes, name=""):

    """
    Function to show in a graphic the distance of IP signal data and INPUT signal data from the TSS.
    Params:
        - arrays_ip: Array: The array of the IP signal data
        - arrays_input: Array: The array of the INPUT signal data
        - xAxes: 
    """

    x = xAxes
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(
        x,
        arrays_ip.mean(axis=0),
        color='r',
        label='IP')

    ax.plot(
        x,
        arrays_input.mean(axis=0),
        color='k',
        label='input')

    ax.axvline(0, linestyle=':', color='k')

    ax.set_xlabel('Distance from TSS (bp)')
    ax.set_ylabel('Average read coverage (per million mapped reads)')
    ax.legend(loc='best')
    ax.set_title(name)
    return fig


##Hacer que homer lo instale el propio usuario
def create_homer_tag_directory(name,file):

    """
    Function to create a HOMER tagDirectory from a signal file.
    Params:
        - name: String: Name of the tagDirectory
        - file: String: PAth opf the file you want to create the HOMER TagDirectory 
    """

    subprocess.run(["makeTagDirectory",name,file])
    return


def homer_peaks(tag_directory,style="factor",output="auto"):

    """
    Function to get the peaks enrichment and the peaks positions into a text file
    Params:
        - tag_directory: String: Path of the HOMER tagDIrectory already created from a signal file.
        - style: String: Style you want to use to calculate the enrichment and the positions
        - output: String: Name of the output text file
    """

    print(output,style)

    subprocess.run(["findPeaks", tag_directory, "-style",style,"-o",output])
    return

def homer_annotate_peaks(file_directory,gene, output):

    """
    Function to get a file with the annotation types of the signal file introduced
    Params:
        - file_directory: String: PAth of the signal file to get the types of annotations. It must be a bed or narrowPeak file
        - gene: String: The gene you want to obtain the annotations (mm10,hg19,...) It must be installed from HOMER before the use
        - output: String: The name of the output text file
    """

    with open(output, "w") as out_file, open("logs.log", "w") as logs_file:
        subprocess.run(["annotatePeaks.pl", file_directory, gene], stdout=out_file, stderr=logs_file)



# def chip_genes_not_used_with_rna(tsses, data, normalized_subtracted, id="gene_name", findBy="gene_name"):

#     """
#     Function to delete the genes not used in both data (chip and rna). Need to use this function to analyse both them together because
#     having more genes can have problems with array length. Important to be in the same order normalized_subtracted and tsses because
#     index are used

#     Params:
#         -tsses: Pybedtool class with the transcriptions
#         -data: DEseq2ResultsPrueba class where there is rna-seq data
#         -normalized_subtracted: Array where have the chip-seq peaks calculated
#         -id: String that indicates the condition filter (gene_name, gene_id,...)
#     """
#     df2 = tsses.to_dataframe()
#     print(df2)
#     gene_used = []
#     gene_not_used_chip = []
#     indexs = []
#     values_not_data = []
#     miidx = -1
#     for idx,value in enumerate(df2.values):
#         ##Hacer dinámico
#         # if miidx != -1:
#         #     continue
#         # else:
#         #     for idx, x in enumerate(value[8].split(";")):
#         #         if findBy in x:
#         #             miidx = idx 
#         #             break

#         # if value[8].split(";")[miidx].split(" ")[1].split('"')[0] in data[id]:
#         #     indexs.append(idx)
#         #     gene_used.append(value[8].split(";")[miidx].split(" ")[1].split('"')[0])
#         # else:
#         #     gene_not_used_chip.append(value[8].split(";")[miidx].split(" ")[1].split('"')[0])
            
#         if value[8].split(";")[3].split(" ")[2].split('"')[1] in data[id]:
#             indexs.append(idx)
#             gene_used.append(value[8].split(";")[3].split(" ")[2].split('"')[1])
#         else:
#             gene_not_used_chip.append(value[8].split(";")[3].split(" ")[2].split('"')[1])


#     with open("genes_not_used_chip_seq.log", "w") as chip_log:
#         for gene in gene_not_used_chip:
#             chip_log.write(f"{gene}\n")

#     for value_data in data.index:
#         if value_data not in gene_used:
#             values_not_data.append(value_data)

#     with open("genes_not_used_rna_seq.log", "w") as rna_log:
#         for gene in values_not_data:
#             rna_log.write(f"{gene}\n")

#     print("\nGenes_not_used_rna-seq: ",values_not_data)
#     normalized_subtracted_good= []
#     for i in range(len(normalized_subtracted)):
#         if i in indexs:
#             normalized_subtracted_good.append(normalized_subtracted[i])
#     data = data.drop(values_not_data)

#     normalized_subtracted= np.array(normalized_subtracted_good)
#     print("Normalized subtract len: ", len(normalized_subtracted))

#     return normalized_subtracted,data

def chip_genes_not_used_with_atac(atac_signal,chip_signal, by):

    """
    Function to eliminate form chip and atac signal, genes that are not coincident between them.
    Params: 
        - atac_signal: Array of dict objects: The array of dict objects with the information of the atac signal
        - chip_signal: Array of dict objects: The array of dict objects with the information of the chip signal
        - by: String: The parameter you want to use to delete those genes (gene_name, gene_id,...)
    """

    atac_signal_good = []
    chip_signal_good = []
    for x in atac_signal:
        for y in chip_signal:
            if x[by] == y[by]:
                atac_signal_good.append(x)
                chip_signal_good.append(y)
                break
    print(" Termine lo alrgo")
    return atac_signal_good, chip_signal_good


def chip_or_atac_genes_not_used_with_rna( data, signal, id="gene_name", rna_column_name="gene_name", delete_no_peaks=False):

    """
    Function to delete the genes not used in both data (chip and rna). Need to use this function to analyse both them together because
    having more genes can have problems with array length. Important to be in the same order normalized_subtracted and tsses because
    index are used

    Params:
        -data: DEseq2ResultsPrueba class: where there is rna-seq data
        -signal: Array of dict objects: where have the chip-seq or atac-seq peaks calculated and with names/ids
        -id: String: indicates the condition filter (gene_name, gene_id,...)
        -rna_column_name: String: The name of the column that wants to match with the id condition
        -delete_no_peaks: Boolean: If true, it will delete the ones values which mean is 0 (no peaks), and if False, not deleting them
    """
    ## Quitando lso q no tienen picos
    #########################################
    if delete_no_peaks == True:
        signal_semigood = []
        for x in signal:
            if np.mean(x['values']) != 0:
                if id =="id":
                    x[id] = x[id].split('.')[0]
                    signal_semigood.append(x)
                else:
                    signal_good.append(x)
    else:
        signal_semigood = signal
    #####################################3
    ###AHora eliminamos aquiellos genes que no existen en rna pero si en atac
    signal_good = []
    ids_selected = []
    for x in signal_semigood:
        # print(x)
        ##Column -> GeneID
        ##SIngal -> id
        # print(id)
        # print(rna_column_name)
        if x[id] in  data.data[rna_column_name].values:
            if x[id] not in ids_selected:
                signal_good.append(x)
                ids_selected.append(x[id])

    print("Len signal",len(signal_good))
    #########################################

    ###Ahora eliminas aquelos genes que no existen en atac pero si en RNA
    indexes = []
    i = 0
    for x in data.data[rna_column_name].values:
        if x in ids_selected:
            indexes.append(i)
        i += 1
        
    data_filtered = data.iloc[indexes]

    data = DEseq2Results(data_filtered)
    print("Data type: ", type(data), "\n data filtered_type", type(data_filtered))
    print("Data lenght after filtered:",len(data))

    with open("genes_not_used_signal.log", "w") as signal_log:
        for gene in ids_selected:
            signal_log.write(f"{gene}\n")

    ##############################################################

    return signal_good,data


    # gene_used = []
    # gene_not_used_chip = []
    # values_not_data = []
    # indexs = []
    # print(normalized_subtracted)
    # for idx, x in enumerate(normalized_subtracted):
    #     # print(x[id], data[id])
    #     ##No se si esta comprobación se debería de hacer
        
    #     if("." in x[id]):
    #         print(x[id])
    #         r = x[id].split(".")[0]
    #     if r in data[id]:
    #         gene_used.append(r)
    #         indexs.append(idx)
    #     else:
    #         gene_not_used_chip.append(r)

        # print(len(normalized_subtracted))

    # with open("genes_not_used_chip_seq.log", "w") as chip_log:
    #     for gene in gene_not_used_chip:
    #         chip_log.write(f"{gene}\n")

    # for value_data in data[id]:
    #     if value_data not in gene_used:
    #         values_not_data.append(value_data)

    # with open("genes_not_used_rna_seq.log", "w") as rna_log:
    #     for gene in values_not_data:
    #         rna_log.write(f"{gene}\n")

    # # print("\nGenes_not_used_rna-seq: ",values_not_data)
    # normalized_subtracted_good= []
    # for i in range(len(normalized_subtracted)):
    #     if i in indexs:
    #         normalized_subtracted_good.append(normalized_subtracted[i])
    # data = data.drop(values_not_data)

    # normalized_subtracted= np.array(normalized_subtracted_good)
    # print("Normalized subtract len: ", len(normalized_subtracted))

    # return signal_good,data


def bigwigToBed(bwFile, bedNameFile):

    """
    Function that converts a BigWig file intro a BEd file. 
    Params:
        -bwFile: String: Path to the BigWig file
        -bedNameFile: String: Name to use in the Bed file created
    """
    bwFile = pyBigWig.open(bwFile)
    chroms = bwFile.chroms()

    bed_data = []
    i = 0
    for chrom in chroms:
        intervals = bwFile.intervals(chrom)
        for start, end, score  in intervals:
            bed_data.append([chrom, start, end, i, score, ""])
            i += 1

    bed_df = pd.DataFrame(bed_data, columns=["chrom", "start", "end", "ID", "value", ""])
    bed_df.to_csv(bedNameFile, sep='\t', header=True, index=False)

def array_gene_name_and_annotation(normalized_subtracted,xAxes):

    """
    Function to analyse CHIP-seq or ATAC-seq data per annotation
    Params:
        - normalized_subtracted: Array of objects with gene information of CHIP-seq or ATAC-seq data
        - xAxes: 

    """

    bed_data = []
    i = 0
    for e in normalized_subtracted:
        i +=1
        bed_data.append([e['chr'], e['start'], e['end'], i ])
    
    df_bed = pd.DataFrame(bed_data, columns=['chr', 'start', 'end', 'ID'])
    df_bed.to_csv("normalized_subtracted.bed",sep='\t', header=True, index=False)

    homer_annotate_peaks(file_directory="normalized_subtracted.bed",gene="mm10", output="normalized_subtracted.bed.txt")
    df_homer = pd.read_csv("normalized_subtracted.bed.txt",delimiter="\t")

    print(df_homer)
    df_homer["Annotation"] = df_homer["Annotation"].str.replace(r'\(.*', '', regex=True)
    print(df_homer)
    unique_values_homer = df_homer["Annotation"].unique()

    array_gene_annotation = []
    for idx, row in df_homer.iterrows():
        array_gene_annotation.append({"gene_name": row['Gene Name'], "annotation": row['Annotation']})

    df_homer_total = {key:[] for key in unique_values_homer}

    for x in normalized_subtracted:
        for y in array_gene_annotation:
            if x['gene_name'] == y['gene_name']:
                df_homer_total[y['annotation']].append(x['values'])
                break
    x = xAxes

    
    for key in df_homer_total:

        if len(df_homer_total[key]) == 0:
            continue
        
        else:
            df_homer_total[key] = np.array(df_homer_total[key])

            heatmap_no_sorted(df_homer_total[key],x)

            heatmap_sorted_by_meanValues(df_homer_total[key],x)

            heatmap_sorted_by_maxValueIndex(df_homer_total[key],x)

            print("AHora viene: ", key)

            plt.show()
            


    print(normalized_subtracted,df_homer_total)

    






def chip_dnase(chip_array, dnase_array,xAxes,name="", xlabel="", ylabel=""):

    """
    Function to show in a graphic the distance of IP signal data and INPUT signal data from the TSS.
    Params:
        - arrays_ip: Array: The array of the IP signal data
        - arrays_input: Array: The array of the INPUT signal data
    """
    x = xAxes

    fig = plt.figure()
    ax = fig.add_subplot(111)
    # Eje principal para DNase-seq
    line1, = ax.plot(
        x,
        dnase_array.mean(axis=0),
        color='k',
        label='DNase-seq'
    )
    ax.set_xlabel('Distance from the TSS')
    ax.set_ylabel('DNase-seq signal intensity')
    ax.axvline(0, linestyle=':', color='k')
    ax.set_title(name) 
    # Eje secundario para ChIP-seq
    ax2 = ax.twinx()
    line2, = ax2.plot(
        x,
        chip_array.mean(axis=0),
        color='r',
        label='ChIP-seq'
    )
    ax2.set_ylabel('ChIP-seq Signal intensity', color='r')
    ax2.tick_params(axis='y', labelcolor='r')

    lines = [line1, line2]
    labels = [line.get_label() for line in lines]
    legend = ax.legend(lines, labels, loc='best', title='Legend')

    return fig


def heatmap_no_sorted(signal_values, xAxis):

    """
    Function to do a heatmap to analyse ATAC-seq or CHIP-seq data without sorting the data
    Params: 
        - signal_values: Array: Array where are stored the signal values
        - xAxis: 
    """

    fig = imshow(
    signal_values,
    x=xAxis,
    figsize=(3, 7),
    percentile=True,
    vmin=5,
    vmax=99,
    line_kwargs=dict(color='k', label='All'),
    fill_kwargs=dict(color='k', alpha=0.3),
    )
    return fig

def heatmap_sorted_by_meanValues(signal_values, xAxis):

    """
    Function to do a heatmap to analyse ATAC-seq or CHIP-seq data sorting the data by the mean of the subintervals of each interval
    Params: 
        - signal_values: Array: Array where are stored the signal values
        - xAxis: 
    """

    fig = imshow(
    signal_values,
    x=xAxis,
    figsize=(3, 7),
    vmin=5, vmax=99,  percentile=True,
    line_kwargs=dict(color='k', label='All'),
    fill_kwargs=dict(color='k', alpha=0.3),
    sort_by=signal_values.mean(axis=1)
    )
    return fig

def heatmap_sorted_by_maxValueIndex(signal_values, xAxis):

    """
    Function to do a heatmap to analyse ATAC-seq or CHIP-seq data sorting the data by the index of the max value of each subinterval
    Params: 
        - signal_values: Array: Array where are stored the signal values
        - xAxis: 
    """

    fig = imshow(
    signal_values,
    x=xAxis,
    figsize=(3, 7),
    vmin=5, vmax=99,  percentile=True,
    line_kwargs=dict(color='k', label='All'),
    fill_kwargs=dict(color='k', alpha=0.3),
    sort_by=np.argmax(signal_values, axis=1),
    )
    return fig

def atac_or_chip_with_rna(signal_values, rna, xAxis,signalName=""):

    """
    Function to analyse and make a figure of CHIP or ATAC data with RNA data
    Params:
        - signal_values: Array: Array where are stored the signal values
        - rna: ResultsTable object: Object with the RNA data. It must have a column called 'log2foldchange' on it to work
        - xAxis:
    """

    x = xAxis
    fig = imshow(
        signal_values,
        x=x,
        figsize=(4, 7),
        vmin=5, vmax=99,  percentile=True,
        line_kwargs=dict(color='k', label='All'),
        fill_kwargs=dict(color='k', alpha=0.3),
        sort_by=signal_values.mean(axis=1),
        height_ratios=(3, 1, 1),
    )
    bottom_axes = plt.subplot(fig.gs[2, 0])

    ci_plot(
        x,
        signal_values[(rna.log2foldchange > 1).values, :],
        line_kwargs=dict(color='#fe9829', label='up'),
        fill_kwargs=dict(color='#fe9829', alpha=0.3),
        ax=bottom_axes)
    ci_plot(
        x,
        signal_values[(rna.log2foldchange < -1).values, :],
        line_kwargs=dict(color='#8e3104', label='down'),
        fill_kwargs=dict(color='#8e3104', alpha=0.3),
        ax=bottom_axes)
    ci_plot(
        x,
        signal_values[((rna.log2foldchange >= -1) & (rna.log2foldchange <= 1)).values, :],
        line_kwargs=dict(color='.5', label='unchanged'),
        fill_kwargs=dict(color='.5', alpha=0.3),
        ax=bottom_axes)
    
    fig.line_axes.set_xticklabels([])
    fig.array_axes.set_xticklabels([])
    fig.line_axes.set_ylabel('Average\nenrichement')
    fig.array_axes.set_ylabel(signalName + ' transcripts')
    bottom_axes.set_ylabel('Average\nenrichment')
    bottom_axes.set_xlabel('Distance from TSS (bp)')
    fig.cax.set_ylabel('Enrichment')

    for ax in [fig.line_axes, fig.array_axes, bottom_axes]:
        ax.axvline(0, linestyle=':', color='k')

    bottom_axes.legend(loc='best', frameon=False, fontsize=8, labelspacing=.3, handletextpad=0.2, bbox_to_anchor=(1.1, 0.5))
    fig.subplots_adjust(left=0.3, right=0.8, bottom=0.05)


    return fig



def all_signal_together(chip_signal_values, atac_signal_values,rna, xAxis,name=""):
    x = xAxis

    fig = imshow(
        chip_signal_values,
        x=x,
        figsize=(4, 7),
        vmin=5, vmax=99,  percentile=True,
        line_kwargs=dict(color='k', label='All'),
        fill_kwargs=dict(color='k', alpha=0.3),
        sort_by=chip_signal_values.mean(axis=1),
        height_ratios=(3, 1, 1)
    )


    bottom_axes = plt.subplot(fig.gs[2, 0])

    ax2 = bottom_axes.twinx()

    ci_plot(
        x,
        chip_signal_values[(rna.log2foldchange > 1).values, :],
        line_kwargs=dict(color='#fe9829', label='Chip-up'),
        fill_kwargs=dict(color='#fe9829', alpha=0.1),
        ax=bottom_axes)


    ci_plot(
        x,
        chip_signal_values[(rna.log2foldchange < -1).values, :],
        line_kwargs=dict(color='#8e3104', label='Chip-down'),
        fill_kwargs=dict(color='#8e3104', alpha=0.1),
        ax=bottom_axes)


    ci_plot(
        x,
        chip_signal_values[((rna.log2foldchange >= -1) & (rna.log2foldchange <= 1)).values, :],
        line_kwargs=dict(color='.5', label='Chip-unchanged'),
        fill_kwargs=dict(color='.5', alpha=0.1),
        ax=bottom_axes)




    ci_plot(
        x,
        atac_signal_values[(rna.log2foldchange > 1).values, :],
        line_kwargs=dict(color='#6c4675', label='Atac-up'),
        fill_kwargs=dict(color='#6c4675', alpha=0.1),
        ax=ax2)


    ci_plot(
        x,
        atac_signal_values[(rna.log2foldchange < -1).values, :],
        line_kwargs=dict(color='#987234', label='Atac-down'),
        fill_kwargs=dict(color='#987234', alpha=0.1),
        ax=ax2)

 
    ci_plot(
        x,
        atac_signal_values[((rna.log2foldchange >= -1) & (rna.log2foldchange <= 1)).values, :],
        line_kwargs=dict(color='k', label='Atac-unchanged'),
        fill_kwargs=dict(color='k', alpha=0.1),
        ax=ax2)

    fig.line_axes.set_xticklabels([])
    fig.array_axes.set_xticklabels([])
    fig.line_axes.set_ylabel('Average\nenrichement')
    fig.array_axes.set_ylabel('Chip-seq transcripts')
    bottom_axes.set_ylabel('Average\nenrichment')
    bottom_axes.set_xlabel('Distance from TSS (bp)')
    fig.cax.set_ylabel('Enrichment')


    for ax in [fig.line_axes, fig.array_axes, bottom_axes]:
        ax.axvline(0, linestyle=':', color='k')


    bottom_axes.legend(loc='best', frameon=False, fontsize=8, labelspacing=.3, handletextpad=0.2,bbox_to_anchor=(1.1, 1))
    ax2.legend(loc='best', frameon=False, fontsize=8, labelspacing=.3, handletextpad=0.2,bbox_to_anchor=(1.1, 0.5))
    fig.subplots_adjust(left=0.3, right=0.8, bottom=0.05)


    return fig



# def chip_dnase_rna_diverse(chip_array, dnase_array,rna,xAxes,name="", xlabel="", ylabel=""):

#     """
#     Function to show in a graphic the distance of IP signal data and INPUT signal data from the TSS.
#     Params:
#         - arrays_ip: Array: The array of the IP signal data
#         - arrays_input: Array: The array of the INPUT signal data
#     """

#     x = xAxes


#     fig = imshow(
#         # Same as before...
#         chip_array,
#         x=x,
#         figsize=(4, 7),
#         vmin=5, vmax=99,  percentile=True,
#         line_kwargs=dict(color='k', label='All'),
#         fill_kwargs=dict(color='k', alpha=0.3),
#         sort_by=chip_array.mean(axis=1),
#         height_ratios=(3, 1, 1)
#     )

#     fig.line_axes.set_xticklabels([])
#     fig.array_axes.set_xticklabels([])
#     fig.line_axes.set_ylabel('Average\nenrichement (ChIP-seq)')
#     fig.array_axes.set_ylabel('Chip-seq transcripts')
#     fig.cax.set_ylabel('Enrichment (ChIP-seq)')
#     fig.subplots_adjust(left=0.3, right=0.8, bottom=0.05)

#     for ax1 in [fig.line_axes, fig.array_axes]:
#         ax1.axvline(0, linestyle=':', color='k')

#     fig2 = imshow(
#         # Same as before...
#         dnase_array,
#         x=x,
#         figsize=(4, 7),
#         vmin=5, vmax=99,  percentile=True,
#         line_kwargs=dict(color='k', label='All'),
#         fill_kwargs=dict(color='k', alpha=0.3),
#         sort_by=dnase_array.mean(axis=1),
#         height_ratios=(3, 1, 1)
#     )
#     fig2.line_axes.set_xticklabels([])
#     fig2.array_axes.set_xticklabels([])
#     fig2.line_axes.set_ylabel('Average\nenrichement )(Dnase-seq)')
#     fig2.array_axes.set_ylabel('DNase-seq transcripts')
#     fig2.cax.set_ylabel('Enrichment (Dnase-seq)')
#     fig2.subplots_adjust(left=0.3, right=0.8, bottom=0.05)

#     for ax2 in [fig2.line_axes, fig2.array_axes]:
#         ax2.axvline(0, linestyle=':', color='k')


#     fig3 = plt.figure()
#     ax = fig3.add_subplot(111)


#     # Eje principal para DNase-seq
#     line1, = ax.plot(
#         x,
#         dnase_array.mean(axis=0),
#         color='k',
#         label='DNase-seq'
#     )
#     ax.set_xlabel('Distance from the TSS')
#     ax.set_ylabel('DNase-seq signal intensity')
#     ax.axvline(0, linestyle=':', color='k')
#     ax.set_title(name)  # Asegúrate de que 'name' esté definido

#     # Eje secundario para ChIP-seq
#     ax4 = ax.twinx()
#     line2, = ax4.plot(
#         x,
#         chip_array.mean(axis=0),
#         color='r',
#         label='ChIP-seq'
#     )
#     ax4.set_ylabel('ChIP-seq Signal intensity', color='r')
#     ax4.tick_params(axis='y', labelcolor='r')

#     ci_plot(
#         x,
#         chip_array[(rna.log2foldchange > 1).values, :],
#         line_kwargs=dict(color='#fe9829', label='Chip-up'),
#         fill_kwargs=dict(color='#fe9829', alpha=0.1),
#         ax=ax4)

#     # Signal over TSSs of transcripts that were repressed upon knockdown
#     ci_plot(
#         x,
#         chip_array[(rna.log2foldchange < -1).values, :],
#         line_kwargs=dict(color='#8e3104', label='Chip-down'),
#         fill_kwargs=dict(color='#8e3104', alpha=0.1),
#         ax=ax4)

#     # Signal over TSSs tof transcripts that did not change upon knockdown
#     ci_plot(
#         x,
#         chip_array[((rna.log2foldchange >= -1) & (rna.log2foldchange <= 1)).values, :],
#         line_kwargs=dict(color='.5', label='Chip-unchanged'),
#         fill_kwargs=dict(color='.5', alpha=0.1),
#         ax=ax4)
#             # Signal over TSSs of transcripts that were activated upon knockdown.



#     ci_plot(
#         x,
#         dnase_array[(rna.log2foldchange > 1).values, :],
#         line_kwargs=dict(color='#6c4675', label='Atac-up'),
#         fill_kwargs=dict(color='#6c4675', alpha=0.1),
#         ax=ax)

#     # Signal over TSSs of transcripts that were repressed upon knockdown
#     ci_plot(
#         x,
#         dnase_array[(rna.log2foldchange < -1).values, :],
#         line_kwargs=dict(color='#987234', label='Atac-down'),
#         fill_kwargs=dict(color='#987234', alpha=0.1),
#         ax=ax)

#     # Signal over TSSs tof transcripts that did not change upon knockdown
#     ci_plot(
#         x,
#         dnase_array[((rna.log2foldchange >= -1) & (rna.log2foldchange <= 1)).values, :],
#         line_kwargs=dict(color='.8', label='Atac-unchanged'),
#         fill_kwargs=dict(color='.8', alpha=0.1),
#         ax=ax)
    
    
#     lines = [line1, line2]
#     labels = [line.get_label() for line in lines]
#     legend = ax.legend(lines, labels, loc='best', title='Legend')
#     # Añadir leyendas para ambos ejes en un solo lugar

    


#     return fig

def chip_atac_rna_diverse(chip, atac, rna, xAxis):
    fig = atac_or_chip_with_rna(rna=rna,signal_values=chip,xAxis=xAxis,signalName="ChIP-seq")
    fig = atac_or_chip_with_rna(rna=rna, signal_values=atac, xAxis=xAxis,signalName="DNase-seq")
    return fig