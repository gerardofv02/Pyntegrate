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

db = gffutils.FeatureDB(os.path.join(path, 'data/gencode.vM25.annotation.gtf.db'))
print(db)

##CArgamos el archivo
mipath = os.path.join(path, 'mosimData/random_atac_seq.narrowPeak')
print(mipath)
atac = Pyntegrate.genomic_signal(mipath, 'narrowpeak')
print(atac)
print(type(atac))
tsses = pybedtools.BedTool(tss_generator()).saveas('tsses.gtf')


tsses_1kb = tsses.slop(b=1000, genome='mm10', output='tsses-1kb.gtf')
processes = multiprocessing.cpu_count()

atac_signal =  atac.array(
    tsses_1kb,
    bins=100,
    processes=processes
)

print(atac_signal)
print(len(atac_signal))

normalized_subtracted_semi = Pyntegrate.chipSeqSignalAnalysis.value_array_simple(atac_signal)

print(normalized_subtracted_semi, len(normalized_subtracted_semi))


filtered_arrays = [arr for arr in normalized_subtracted_semi if np.mean(arr) != 0]

##COnvertimos en array de numpy para poder hacer la ordenación
normalized_subtracted_complete = np.array(filtered_arrays)
print(len(normalized_subtracted_complete))
x = np.linspace(-1000, 1000, bins)



##Primera figura sin ordenación
fig = Pyntegrate.plotutils.imshow(

    # The array to plot; here, we've subtracted input from IP.
    normalized_subtracted_complete,

    # X-axis to use
    x=x,

    # Change the default figure size to something smaller for this example
    figsize=(3, 7),

    # Make the colorbar limits go from 5th to 99th percentile.
    # `percentile=True` means treat vmin/vmax as percentiles rather than
    # actual values.
    percentile=True,
    vmin=5,
    vmax=99,

    # Style for the average line plot (black line)
    line_kwargs=dict(color='k', label='All'),

    # Style for the +/- 95% CI band surrounding the
    # average line (transparent black)
    fill_kwargs=dict(color='k', alpha=0.3),
)



##Segunda y cuarta figura con ordenación de que se muestran la media de los subintervalos

fig = Pyntegrate.plotutils.imshow(

    # These are the same arguments as above.
    normalized_subtracted_complete,
    x=x,
    figsize=(3, 7),
    vmin=5, vmax=99,  percentile=True,
    line_kwargs=dict(color='k', label='All'),
    fill_kwargs=dict(color='k', alpha=0.3),

    # This is new: sort by mean signal
    sort_by=normalized_subtracted_complete.mean(axis=1)
)

##TErcera figurra que se muestra la señal más altas de los subintervalos

fig = Pyntegrate.plotutils.imshow(

    # These are the same arguments as above.
    normalized_subtracted_complete,
    x=x,
    figsize=(3, 7),
    vmin=5, vmax=99,  percentile=True,
    line_kwargs=dict(color='k', label='All'),
    fill_kwargs=dict(color='k', alpha=0.3),

    # This is new: sort by mean signal
    sort_by=np.argmax(normalized_subtracted_complete, axis=1)
)



fig = Pyntegrate.plotutils.imshow(
    normalized_subtracted_complete,
    x=x,
    figsize=(3, 7),
    vmin=5, vmax=99,  percentile=True,
    line_kwargs=dict(color='k', label='All'),
    fill_kwargs=dict(color='k', alpha=0.3),
    sort_by=normalized_subtracted_complete.mean(axis=1)
)

plt.show()











# import pandas as pd
# import os

# def check_narrowpeak_format(filepath):
#     try:
#         # Lee el archivo usando pandas para verificar el formato
#         df = pd.read_csv(filepath, sep='\t', header=None)

#         # Imprime las primeras filas y los tipos de datos de cada columna
#         print("Primeras filas del archivo:")
#         print(df.head())

#         print("\nTipos de datos de cada columna:")
#         print(df.dtypes)

#         # Verifica que las columnas de coordenadas sean enteras
#         if not pd.api.types.is_integer_dtype(df[1]) or not pd.api.types.is_integer_dtype(df[2]):
#             raise ValueError("Las coordenadas cromosómicas deben ser enteras.")
        
#         # Verificar si todas las filas tienen exactamente 10 columnas
#         if not (df.shape[1] == 10):
#             raise ValueError("El archivo debe tener exactamente 10 columnas.")
        
#         print("El archivo está en el formato correcto.")
#     except Exception as e:
#         print(f"Error en el formato del archivo: {e}")

# # Verificar el archivo narrowPeak
# path = os.path.dirname(os.path.abspath(__file__))
# mipath = os.path.join(path, 'mosimData/related_atac_seq.narrowPeak')
# check_narrowpeak_format(mipath)



# import pybedtools

# # Define un intervalo de prueba
# test_interval = pybedtools.Interval('chr1', 10000, 10500)

# # Carga el archivo narrowPeak como un BedTool
# bedtool = pybedtools.BedTool('mosimData/related_atac_seq.narrowPeak')

# # Intenta hacer una intersección simple
# try:
#     result = bedtool.intersect(test_interval, wa=True)
#     print("Intersección exitosa, número de intervalos encontrados:", len(result))
#     for interval in result:
#         print(interval)
# except Exception as e:
#     print(f"Error al intentar la intersección: {e}")



# import pybedtools
# pybedtools.set_tempdir('tmp/')
# import pandas as pd

# # Cargar el archivo narrowPeak
# narrowpeak_df = pd.read_csv('mosimData/related_atac_seq.narrowPeak', sep='\t', header=None)

# # Verificar los cromosomas y los rangos
# print("Cromosomas disponibles:", narrowpeak_df[0].unique())
# print("\nRango de posiciones por cromosoma:")
# for chrom in narrowpeak_df[0].unique():
#     chrom_data = narrowpeak_df[narrowpeak_df[0] == chrom]
#     min_pos = chrom_data[1].min()
#     max_pos = chrom_data[2].max()
#     print(f"{chrom}: {min_pos}-{max_pos}")

# import pybedtools

# # Define un intervalo válido
# valid_chrom = 'chr10'
# valid_start = 3972940  # Asegúrate de que el inicio está dentro del rango
# valid_stop = 128977854   # Asegúrate de que el fin está dentro del rango
# # Crear el intervalo como un BedTool desde la cadena
# interval_str = f"{valid_chrom}\t{valid_start}\t{valid_stop}"
# query_interval = pybedtools.BedTool(interval_str, from_string=True)

# # Cargar el archivo narrowPeak como un BedTool
# bedtool = pybedtools.BedTool('mosimData/related_atac_seq.narrowPeak')

# # Intentar hacer una intersección simple
# try:
#     result = bedtool.intersect(query_interval, wa=True)
#     print("Intersección exitosa, número de intervalos encontrados:", len(result))
#     for i, interval in enumerate(result):
#         if i < 10:  # Mostrar solo los primeros 10 para no abrumar
#             print(interval)
# except Exception as e:
#     print(f"Error al intentar la intersección: {e}")


