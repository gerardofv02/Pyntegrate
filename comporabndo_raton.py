import pyBigWig

def verificar_raton(bigwig_file):
    try:
        # Abrir el archivo BigWig
        bw = pyBigWig.open(bigwig_file)

        # Obtener la información del encabezado
        header = bw.chroms()
        
        # Imprimir la información del encabezado
        print("Información del archivo BigWig:")
        for chrom in header:
            print(f"Chromosome: {chrom}, Length: {header[chrom]}")
        
        # Buscar si hay referencias a Mus musculus (Esto depende del contenido del archivo)
        # Si el archivo contiene nombres de cromosomas específicos de ratón, podrías buscar "chr" seguido de números
        mus_musculus_present = any("chr" in chrom for chrom in header)
        
        if mus_musculus_present:
            print("El archivo BigWig contiene datos de Mus musculus.")
        else:
            print("El archivo BigWig no contiene datos de Mus musculus o no se puede determinar con la información disponible.")
        
        # Cerrar el archivo BigWig
        bw.close()
    except Exception as e:
        print(f"Error al abrir el archivo BigWig: {e}")

# Ruta al archivo BigWig
bigwig_file = "./BamFiles/SRR3134986_nodup.bw"
verificar_raton(bigwig_file)
