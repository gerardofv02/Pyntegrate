import Pyntegrate
import multiprocessing
from matplotlib import pyplot as plt
import matplotlib
import numpy as np
import Pyntegrate
import pybedtools


# ##Supported_formats()
# at = Pyntegrate._genomic_signal.supported_formats()
# print(at)

###Genomic_signal(bam)
# gs = Pyntegrate.genomic_signal(Pyntegrate.example_filename('x.bam'), 'bam')
# print(gs)

# ##Genomic_signal(bed)
# gs = Pyntegrate.genomic_signal(Pyntegrate.example_filename('gdc.bed'), 'bed')
# print(gs)

##ver para hacer local_coverage

# ##Genomic_signal(bigbed) ERROR REVISAR
# gs = Pyntegrate.genomic_signal(Pyntegrate.example_filename('gdc.bigbed'), 'bigbed')
# print(gs)

##Genomic_signal(bigwig)
gs = Pyntegrate.genomic_signal(Pyntegrate.example_filename('gdc.bigwig'), 'bigwig')
print(gs)

# ##Genomic_signal(gtf(objecto bed))
# gs = Pyntegrate.genomic_signal("Homo_sapiens.GRCh38.109.gtf","gtf")
# print(gs)


"""Probando este archivo se ha posiso también comprobar que el archivo de array helpers también funcionan debido a que es el único archivo que usa estas funciones y con lo cual si funciona aqui, funciona alli"""