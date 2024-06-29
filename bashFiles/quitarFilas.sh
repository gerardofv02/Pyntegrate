#!/bin/bash

# Nombre del archivo GTF
archivo="Homo_sapiens.GRCh38.109.gtf"
# Número de línea a partir del cual eliminar
linea_a_eliminar=316886

# Obtener el número total de líneas en el archivo
num_lineas=$(wc -l < "$archivo")

# Verificar si el número de línea a eliminar es válido
if [ "$linea_a_eliminar" -le "$num_lineas" ]; then
    # Crear un archivo temporal con las líneas anteriores a la línea a eliminar
    head -n "$((linea_a_eliminar - 1))" "$archivo" > temporal.gtf
    
    # Reemplazar el archivo original con el archivo temporal
    mv temporal.gtf "$archivo"
    echo "Se eliminaron todas las líneas a partir de la línea $linea_a_eliminar."
else
    echo "El número de línea a eliminar es mayor que el número total de líneas en el archivo."
fi