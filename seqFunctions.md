# Signal Functions

Once the signals from our files have been calculated, we can proceed to explain the variety of functions available in the Pyntegrate library and how to use them. Most of these functions are found within the SeqSignalAnalysis module. Additionally, each of these functions has its own documentation. To view it, you can use the following Python script:
```python
help(Pyntegrate.FilesFunction.Function)
```
The function that calculates the signals for the data types returns an array of intervals, where each element is an object with the following fields:

- Values: Signal values
- gene_name: Name of the gene in the interval
- gene_id: ID of the gene in the interval
- strand: Indicates the strand from which the interval was obtained
- chr: Chromosome of the interval
- start: Start of the interval
- stop: End of the interval

To begin, we'll explain how to clean the signal data to enable joint analysis between different types. Afterward, we’ll describe the functions used to generate graphs and perform the analysis.

## Data Cleaning Functions

### ChIP-seq/RNA-seq or DNase-seq/RNA-seq

This function is designed to return the intersecting genes between these types. Since ChIP-seq and DNase-seq share similar data processing, the same function can be used for both to find the intersecting genes.

This function is called **chip_or_atac_genes_not_used_with_rna** and has the following input parameters:

- data: DEseq2ResultsPrueba class or similar: Stores the RNA-seq data
- signal: Array of objects: Stores the array of objects mentioned earlier, containing the signal information
- id: String: Specifies the value to filter by (id, gene_name, etc.)
- rna_column_name: String: Name of the RNA column for filtering
- delete_no_peaks: Boolean: If True, removes genes from the signal where the mean signal value is 0

### Signal Values

In some functions, it is necessary to only add the signal values from ChIP-seq and ATAC-seq. This function returns an array containing only the signal values, without the object.

This function is called **value_array_simple** and has the following input parameter:

- array: Array of objects: Stores the array of objects mentioned earlier, containing the signal information

###  ChIP-seq/DNase-seq

This function is designed to return the intersecting genes between ChIP-seq and DNase-seq.

This function is called **chip_genes_not_used_with_atac** and has the following input parameters:

- atac_signal: Array of objects: Stores the array of objects mentioned earlier, containing the DNase signal information
- chip_signal: Array of objects: Stores the array of objects mentioned earlier, containing the ChIP signal information
- by: String: Specifies the value to filter by (id, gene_name, etc.)

### Calculating ChIP-seq Peaks

When you have both an IP and INPUT file in ChIP-seq, a specific function in Pyntegrate can calculate the peaks.
This function is called **calculate_peaks_with_gene_name** and has the following input parameters:

- arrays_ip: Array of objects: Stores the array of objects mentioned earlier, containing the IP signal information
- arrays_input: Array of objects: Stores the array of objects mentioned earlier, containing the INPUT signal information

### Bigwig to Bed Converter

This function was created for internal use in other functions, but it can also be used directly. It converts a bigwig file to bed format, which is necessary for some HOMER functions.
This function is called **bigwigToBed** and has the following input parameters:

- bwFile: String: Path to the bigwig file
- bedNameFile: String: Name of the bed file to create

## Data Analysis and Graph Generation Functions

### IP, INPUT Distance from TSS Graph

This function returns a graph showing the distance of the IP and INPUT signals from TSS.
The function is named **distance_from_tss_chipSeq** and has the following input parameters:

- arrays_ip: Array: Array containing the IP signal values
- arrays_input: Array: Array containing the INPUT signal values
- xAxes: Sequence of equally spaced numbers (numpy.linespace)

### Heatmaps by Annotation Type

This function generates heatmaps for different types of DNA regions. It is useful for analyzing different types, such as promoters, exons, etc. The function is named **annotation_type_heatmap** and has the following input parameters:

- normalized_subtracted: Array of objects: Stores the array of objects mentioned earlier, containing the signal information (either ChIP or DNase)
- xAxes: Sequence of equally spaced numbers (numpy.linespace)

### ChIP-DNASE

This function was created to allow joint analysis of ChIP-seq and DNase-seq data. It generates a graph where the y-axis shows the values for both ChIP and DNase with their respective scales, and the x-axis shows the distance to TSS.
The function is named **chip_dnase** and has the following input parameters:

- chip_array: Array of objects: Stores the array of objects mentioned earlier, containing the ChIP signal information
- dnase_array: Array of objects: Stores the array of objects mentioned earlier, containing the DNase signal information
- xAxes: Sequence of equally spaced numbers (numpy.linespace)

### Heatmaps

These heatmaps are used to analyze ChIP-seq and DNase-seq. There are three types of heatmaps, which vary in the order in which values are presented in the heatmap. These functions are called **heatmap_no_sorted**, **heatmap_sorted_by_meanValues**, and **heatmap_sorted_by_maxValueIndex** and have the following input parameters:

- signal_values: Array of objects: Stores the array of objects mentioned earlier, containing the signal information (either ChIP or DNase)
- xAxis: Sequence of equally spaced numbers (numpy.linespace)

###  ChIP-seq/RNA-seq or DNase-seq/RNA-seq

This function graphically represents the values of ChIP-seq or DNase-seq with RNA-seq. It displays three graphs showing the heatmap with the ChIP or DNase signal, the enrichment of this signal, and the combination with RNA-seq.
This function is called **atac_or_chip_with_rna** and has the following input parameters:

- signal_values: Array of objects: Stores the array of objects mentioned earlier, containing the signal information (either ChIP or DNase)
- rna: DEseq2ResultsPrueba class or similar: Stores the RNA-seq data
- xAxis: Sequence of equally spaced numbers (numpy.linespace)

### All Signals Together

This function was created to analyze all data types combined. To do this, a 3D graph is generated.
The function is named **all_signal_together_new** and has the following input parameters:

- chip_signal_values: Array of objects: Stores the array of objects mentioned earlier, containing the ChIP-seq signal information
- atac_signal_values: Array of objects: Stores the array of objects mentioned earlier, containing the DNase-seq signal information
- rna: DEseq2ResultsPrueba class or similar: Stores the RNA-seq data
