import Pyntegrate
from Pyntegrate.arrayify import Binner

##Genomic_signal(bigwig)
gs = Pyntegrate.genomic_signal(Pyntegrate.example_filename('gdc.bigwig'), 'bigwig')
print(gs)

x = Binner(genome="hg38",windowsize=1000)

x.to_npz(gs)