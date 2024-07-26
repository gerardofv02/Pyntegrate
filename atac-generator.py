# from ncls import NCLS
# import numpy as np
# import matplotlib.pyplot as plt

# # Generar intervalos aleatorios
# num_intervals = 1000
# genome_size = 1000000
# starts = np.random.randint(0, genome_size, num_intervals)
# ends = starts + np.random.randint(100, 500, num_intervals)
# ids = np.arange(num_intervals)

# # Crear un NCLS
# ncls = NCLS(starts, ends, ids)

# # Simular lecturas ATAC-seq
# num_reads = 5000
# read_starts = np.random.randint(0, genome_size, num_reads)
# read_ends = read_starts + np.random.randint(50, 200, num_reads)
# read_ids = np.arange(num_reads)

# # Encontrar intersecciones
# overlaps = ncls.all_overlaps_both(read_starts, read_ends, read_ids)

# print(overlaps)

# # Visualizar los intervalos y las lecturas
# plt.hist(starts, bins=100, alpha=0.5, label='Intervals')
# plt.hist(read_starts, bins=100, alpha=0.5, label='Reads')
# plt.xlabel('Genomic position')
# plt.ylabel('Count')
# plt.legend()
# plt.title('Simulated ATAC-seq intervals and reads')
# plt.show()

# # Imprimir algunos resultados de las intersecciones
# for overlap in overlaps[:10]:
#     interval_start = overlap
#     print(f"Interval {interval_start}")


















# Importar bibliotecas necesarias
import ncls
import numpy as np
 
# Definir el genoma de referencia (Ejemplo: longitud de los cromosomas)
genome_size = 3e9  # Tamaño aproximado del genoma humano
 
# Definir el número de picos de accesibilidad y su longitud
num_peaks = 10000
peak_length = 200  # Longitud promedio de los picos de ATAC-seq
 
# Generar posiciones aleatorias para los picos
peak_positions = np.random.randint(0, genome_size, num_peaks)
 
# Crear un archivo BED simulado con los picos de ATAC-seq
with open("simulated_atacseq.bed", "w") as bed_file:
    for pos in peak_positions:
        bed_file.write(f"chr1\t{pos}\t{pos + peak_length}\n")
 
# Usar NCLS para analizar los datos simulados
# Aquí deberías usar las funciones específicas de NCLS para cargar y analizar los datos
# Ejemplo (dependiendo de la implementación de NCLS):
ncls_data = NCLS.load_bed("simulated_atacseq.bed")
# Realizar análisis adicional según sea necesario






# import numpy as np
 
# # Definir los parámetros para la simulación
# genome_size = 3e9  # tamaño aproximado del genoma humano
# num_peaks_control = 10000
# num_peaks_experimental = 15000
# peak_length = 500  # longitud promedio de los picos de ATAC-seq
 
# # Función para generar picos
# def generate_peaks(num_peaks):
#     positions = np.random.randint(0, genome_size - peak_length, size=num_peaks)
#     return [("chr1", pos, pos + peak_length) for pos in positions]
 
# # Generar datos para control y experimental
# control_peaks = generate_peaks(num_peaks_control)
# experimental_peaks = generate_peaks(num_peaks_experimental)
 
# # Guardar los datos en formato BED
# def save_peaks(peaks, filename):
#     with open(filename, "w") as file:
#         for chr, start, end in peaks:
#             file.write(f"{chr}\t{start}\t{end}\n")
 
# save_peaks(control_peaks, "control_peaks.bed")
# save_peaks(experimental_peaks, "experimental_peaks.bed")