#!/bin/bash

# Nombre del archivo de entrada
input_file="SRR1204544_sort_nondup.bw.bed"

# Nombre del archivo temporal
temp_file="temp.bed"

# Eliminar todo a partir de la línea 50000
head -n 49999 "$input_file" > "$temp_file"

# Reemplazar el archivo original con el archivo temporal
mv "$temp_file" "$input_file"