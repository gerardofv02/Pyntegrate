import pysam

# Ruta al archivo BAM
bam_file = "./BamFiles/SRR1204546_sort_nondup.bam"

# Generar el índice
pysam.index(bam_file)