# Pyntegrate

Pyntegrate is a library for genomic data analysis, built on the foundation of another library called [metaseq](https://github.com/daler/metaseq)  created by Ryan Dale. The update of this library is driven by the need to analyze a wider variety of data types. The types of data supported by Pyntegrate are:

1.  ChIP-seq
2.  RNA-seq
3.  DNase-seq/ATAC-seq

With these types of data accepted by the library, various functions have been created to enable the generation of graphs with them and to allow for cross-analysis.

This library also utilizes software called [HOMER](http://homer.ucsd.edu/homer/).It is especially used for identifying ChIP-seq peak enrichments and distinguishing between different types of DNA regions based on intervals.


## Pyntegrate Installation Instructions

To install Pyntegrate, you need to have Git installed.
If you are using Linux, you should run the following commands:
```console
sudo apt update
sudo apt install git
```

Once Git is installed, is needed to have some essentials tools installed. For that, the following commando should be executed:
```console
sudo apt-get install build-essential
```

Once the essentials tools are installed, the following command should be executed:

```console
pip install git+https://github.com/gerardofv02/Pyntegrate.git
```
For the full fucntionality, is required to install bedtools. For that, there are two options:
- Using conda:
  ```console
  conda install -c bioconda bedtools
  ```
- Not using conda:
  ```console
  sudo apt-get install bedtools
  ```

## HOMER Installation

To use the HOMER-dependent functions (which are optional for using the library), you should follow the installation instructions available at [this link](http://homer.ucsd.edu/homer/introduction/install.html)

## Getting Started

To start using Pyntegrate, you need genomic data files with the following extensions:

- bigwig, bam: For full functionality (ChIP-seq, ATAC-seq)
- bed, bigbed, gff, gtf, vcf: These provide most functionality but do not support read counting from your file (ChIP-seq, ATAC-seq)
- csv: For full functionality (RNA-seq)

Once you have this data, to obtain signal values in ChIP-seq and ATAC-seq and load them as objects within your code, you should run the following function:

```python
features, signal, tsses, tsses_1kb = Pyntegrate.SeqSignalAnalysis.generate_array_simple_signal(dbPath='/path/to/transcript/data',
                                                                                                   filePath='/path/to/signal/data',
                                                                                                   extensionFile='name-of-extension-file',
                                                                                                   genome='genome',
                                                                                                   bins="number-of-subintervals-want-to-generate")
```
To obtain RNA-seq values, you should create the following object, which will be stored as that object within the variable (for example, of type DEseq2):

```python
data = Pyntegrate.results_table.DEseq2Results('/path/to/RNA.csv')
```

Once you have these variables, a wide range of functions can be performed with the signal data. These functions are better described in [this link](./seqFunctions.md)

## Usage Examples

Below are some usage examples that may be very useful for understanding the functionality and extensibility of the library:

### [Example 1](./example_of_use/Example_1_Pyntegrate.ipynb)

In this example, ChIP-seq analysis will be performed initially. This will include peak calculation and heatmap generation to enable analysis. Later, RNA-seq will be integrated and analyzed together with ChIP-seq. An example of the generated graph is as follows:

![CHIP-RNA](./images_examples/CHIP-RNA.png)

This example is fully available in this repository. The link is: [Example_1_Pyntegrate](./example_of_use/Example_1_Pyntegrate.ipynb)

### [Example 2](./example_of_use/Example_2_Pyntegrate.ipynb)

In this example, DNase-seq analysis will be performed initially. Heatmaps will be generated based on signal values for subsequent joint analysis with RNA-seq. An example of the generated graph (similar to the ChIP-seq and RNA-seq combination) is as follows:

![DNase-RNA](./images_examples/DNase-RNA.png)

This example is fully available in this repository. The link is: [Example_2_Pyntegrate](./example_of_use/Example_2_Pyntegrate.ipynb)

### [Example 3](./example_of_use/Example_3_Pyntegrate.ipynb)

In this example, ChIP-seq will be analyzed together with DNase-seq. The separate analysis of each of these types can be found in the two previous examples. For the combined analysis of these types, a graph will be generated (using a Pyntegrate function) as shown below:

![CHIP-DNASE](./images_examples/CHIP-DNASE.png)

This example is fully available in this repository. The link is: [Example_3_Pyntegrate](./example_of_use/Example_3_Pyntegrate.ipynb)

### [Example 4](./example_of_use/Example_4_Pyntegrate.ipynb)

In this example, all data types (ChIP-seq, RNA-seq, and DNase-seq) will be analyzed together to get a complete view of the process. Two different functions have been generated for this purpose: one for analyzing pairs of data types separately, and another for combining all data types and analyzing them together using a 3D graph. An example of this graph is shown below:

![CHIP-DNASE-RNA](./images_examples/CHIP-DNASE-RNA.png)

This example is fully available in this repository. The link is: [Example_4_Pyntegrate](./example_of_use/Example_4_Pyntegrate.ipynb)

### [Example 5](./example_of_use/Example_5_Pyntegrate.ipynb)

In this example, a volcano plot will be generated to analyze RNA-seq data. This plot will show which genes are repressed, expressed, and normal. Additionally, if desired, the gene names can be displayed on the plot. An example of this graph is shown below:

![Volcano plot](./images_examples/Volcanoplot.png)

This example is fully available in this repository. The link is: [Example_5_Pyntegrate](./example_of_use/Example_5_Pyntegrate.ipynb)

### [Example 6](./example_of_use/Example_6_Pyntegrate.py)

In this example, a circular diagram will be generated. By passing a file created by the HOMER function that shows annotation types, you will be able to visualize the different percentages of various annotation types in that file. An example of this graph is shown below:

![Pie Chart annotations](./images_examples/PieChartAnnotations.png)

This example is fully available in this repository. The link is: [Example_6_Pyntegrate](./example_of_use/Example_6_Pyntegrate.py)

### [Example 7](./example_of_use/Example_7_Pyntegrate.py)

In this example, heatmaps will be generated repetitively by varying the data. The data from which the heatmaps are generated are annotation types. That is, different heatmaps will be generated according to the type of annotation, especially for analyzing types separately. An example of the generated heatmaps is shown below (in this case, it is for Promoter):

![Promoter heatmaps](./images_examples/Promoter.png)

This example is fully available in this repository. The link is: [Example_7_Pyntegrate](./example_of_use/Example_7_Pyntegrate.py)
