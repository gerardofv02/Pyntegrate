import Pyntegrate
import os
import gffutils
import pybedtools
from pybedtools.featurefuncs import TSS
from gffutils.helpers import asinterval
import multiprocessing


def tss_generator():
    """
    Generator function to yield TSS of each annotated transcript
    """
    for transcript in db.features_of_type('transcript'):
        yield TSS(asinterval(transcript), upstream=1, downstream=0)

path = os.path.dirname(os.path.abspath(__file__))
print(path)

# db = gffutils.create_db(data=path+"/Homo_sapiens.GRCh38.109.gtf",dbfn=path+"/data/Homo_sapiens.GRCh38.109.gtf.db")

db = gffutils.FeatureDB(os.path.join(path, 'data/Homo_sapiens.GRCh38.109.gtf.db'))

tsses = pybedtools.BedTool(tss_generator()).saveas('tsses.gtf')

tsses_1kb = tsses.slop(b=1000, genome='hg38', output='tsses-1kb.gtf')

ip_signal = Pyntegrate.genomic_signal(
    os.path.join(path, 'SRR3134987_nodup.bw'),
    'bigwig')

input_signal = Pyntegrate.genomic_signal(
    os.path.join(path, 'SRR3134986_nodup.bw'),
    'bigwig')

print(ip_signal,input_signal)
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

    # Do the same thing for input.
    input_array = input_signal.array(
        tsses_1kb,
        bins=100,
        processes=processes)

    # Normalize to library size. The values in the array
    # will be in units of "reads per million mapped reads"
    ip_array /= ip_signal.mapped_read_count() / 1e6
    input_array /= input_signal.mapped_read_count() / 1e6

    # Cache to disk. The data will be saved as "example.npz" and "example.features".
    Pyntegrate.persistence.save_features_and_arrays(
        features=tsses,
        arrays={'ip': ip_array, 'input': input_array},
        prefix='example',
        link_features=True,
        overwrite=True)
# t = Pyntegrate.results_table.DEseq2ResultsPrueba("./desq2_Sanchez.csv")
# print(t.disenriched())

