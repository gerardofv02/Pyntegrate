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
chip = Pyntegrate.genomic_signal('test_3_ChIP-seq_3.bw', 'bigwig')
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

for x in chip_signal:
    chip_signal_good.append(x[0])

chip_signal_good = np.array(chip_signal_good)
normalized_subtracted_semi = Pyntegrate.chipSeqSignalAnalysis.value_array_simple(chip_signal)

normalized_subtracted_semi /= chip.mapped_read_count() / 1e6

##Suponiendo que aquellos valores donde la media sea 0, no existen picos con lo cual los elimino

filtered_arrays = [arr for arr in normalized_subtracted_semi if np.mean(arr) != 0]

normalized_subtracted_complete = np.array(filtered_arrays)
print(normalized_subtracted_semi)
x = np.linspace(-1000, 1000, bins)


##Cargamos el rna-seq

rna = Pyntegrate.results_table.ResultsTable('test_3_RNA-seq3.csv')

rna.data['gene_name'] = rna.data['Unnamed: 0']

rna.reindex_to(tsses, attribute='gene_name')
rna.data['id'] = rna['gene_name']

rna.delete_genes_no_peaks_chip(chip_signal_good)

norma_sub_bueno, data = Pyntegrate.chipSeqSignalAnalysis.genes_not_used_with_rna_new( rna, chip_signal_good, id="id", findBy="gene_id")

##No estoy encontrando datos conjuntos para hacer análisis juntos con estos datos pero si que puedo abrirlos bien y realizar bien las funciones. De todos modso he creado otra fuinción para eliminar genes que no estén en ambos lados y es mejor

# print(norma_sub_bueno,rna)
# print(rna)

