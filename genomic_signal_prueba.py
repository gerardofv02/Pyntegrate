import Pyntegrate


###Supported_formats()
# at = Pyntegrate._genomic_signal.supported_formats()
# print(at)

###Genomic_signal(bam)
# gs = Pyntegrate.genomic_signal(Pyntegrate.example_filename('x.bam'), 'bam')
# print(gs)

###Genomic_signal(bed)
# gs = Pyntegrate.genomic_signal(Pyntegrate.example_filename('gdc.bed'), 'bed')
# print(gs)

###Genomic_signal(bigbed) ERROR REVISAR
# gs = Pyntegrate.genomic_signal(Pyntegrate.example_filename('gdc.bigbed'), 'bigbed')
# print(gs)

###Genomic_signal(bigwig)
gs = Pyntegrate.genomic_signal(Pyntegrate.example_filename('gdc.bigwig'), 'bigwig')
print(gs)