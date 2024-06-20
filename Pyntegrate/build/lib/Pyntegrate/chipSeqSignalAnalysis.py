from ._genomic_signal import *
from .persistence import *
from .array_helpers import *
import os
import gffutils
import pybedtools
from pybedtools.featurefuncs import TSS
from gffutils.helpers import asinterval
import multiprocessing
import numpy as np

def tss_generator(db):
    """
    Generator function to yield TSS of each annotated transcript
    """
    for transcript in db.features_of_type('transcript'):
        yield TSS(asinterval(transcript), upstream=1, downstream=0)

def generate_arrays_features_from_tsses_from_db(dbPath, ipSignalPath, extensionIp,inputSignalPath,extensionInput, genome):
    
    ##Primer paso 
    db = gffutils.FeatureDB(dbPath) #Primer paso
    tsses = pybedtools.BedTool(tss_generator(db)).saveas('tsses.gtf')

    remove_duplicates("tsses.gtf")
    tsses = pybedtools.BedTool("tsses.gtf")
    tsses_1kb = tsses.slop(b=1000, genome=genome, output='tsses-1kb.gtf')
    


    ##Segundo paso
    ip_signal = genomic_signal(ipSignalPath,extensionIp)
    input_signal = genomic_signal(inputSignalPath,extensionInput)
    
    
    

    print(ip_signal,input_signal)
    processes = multiprocessing.cpu_count()
    if not os.path.exists('example.npz'):

        ##Tercer paso

        # The signal is the IP ChIP-seq BAM file.
        ip_array = ip_signal.array(

            # Look at signal over these windows
            tsses_1kb,

            # Bin signal into this many bins per window
            bins=100,

            # Use multiple CPUs. Dramatically speeds up run time.
            processes=processes)

        # print(ip_array[0][:10])
        # Do the same thing for input.

        ##Cuartopaso
        input_array = input_signal.array(
            tsses_1kb,
            bins=100,
            processes=processes)

        # print(input_array[:10])



        # Normalize to library size. The values in the array
        # will be in units of "reads per million mapped reads"
        ##SExto paso

        ip_array /= ip_signal.mapped_read_count() / 1e6
        input_array /= input_signal.mapped_read_count() / 1e6


        print("Ip_array: " , ip_array)
        print("Input_array: " , input_array)

        ##SExto paso
        # Cache to disk. The data will be saved as "example.npz" and "example.features".
        save_features_and_arrays(
            features=tsses,
            arrays={'ip': ip_array, 'input': input_array},
            prefix='example',
            link_features=True,
            overwrite=True)

    features, arrays = load_features_and_arrays(prefix='example')
    return features ,arrays ,tsses, tsses_1kb