import pyBigWig
import numpy as np
import pandas as pd

# Parámetros para la simulación de ATAC-seq
fragment_length_mean = 150  # Longitud promedio del fragmento en ATAC-seq
fragment_length_std = 30    # Desviación estándar de la longitud del fragmento

# Leer el archivo bigWig de DNase-seq
dnase_bw = pyBigWig.open('test_3_DNase_4.bw')

# Obtener todos los cromosomas
chromosomes = dnase_bw.chroms().keys()

def simulate_atac_fragments(chrom, values, length_mean, length_std):
    simulated_fragments = []
    
    # Encontrar picos significativos
    peaks = []
    start = None
    for i, value in enumerate(values):
        if value > 0:
            if start is None:
                start = i
        else:
            if start is not None:
                peaks.append((start, i))
                start = None
    if start is not None:
        peaks.append((start, len(values)))
    
    # Generar fragmentos simulados
    for start, end in peaks:
        num_fragments = np.random.randint(3, 10)
        
        for _ in range(num_fragments):
            frag_length = int(np.random.normal(length_mean, length_std))
            frag_length = max(frag_length, 50)  # Asegurarse de que el fragmento no sea demasiado corto
            
            # Asegurarse de que el rango sea válido
            if end - start > frag_length:
                frag_start = np.random.randint(start, end - frag_length)
                frag_end = frag_start + frag_length
                
                simulated_fragments.append([chrom, frag_start, frag_end])
    
    return simulated_fragments

all_simulated_fragments = []

# Procesar cada cromosoma
for chrom in chromosomes:
    values = dnase_bw.values(chrom, 0, dnase_bw.chroms()[chrom])
    fragments = simulate_atac_fragments(chrom, values, fragment_length_mean, fragment_length_std)
    all_simulated_fragments.extend(fragments)

# Convertir a DataFrame
atac_simulated_df = pd.DataFrame(all_simulated_fragments, columns=['chrom', 'start', 'end'])

# Guardar los datos simulados en un archivo BED
atac_simulated_df.to_csv('simulated_atac_peaks.bed', sep='\t', header=False, index=False)

# Cerrar el archivo bigWig
dnase_bw.close()
