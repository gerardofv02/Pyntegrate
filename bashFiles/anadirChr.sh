#!/bin/bash

# Verifica si se proporcionó un archivo como argumento
if [ $# -ne 1 ]; then
    echo "Uso: $0 archivo.gtf"
    exit 1
fi

# Verifica si el archivo existe
if [ ! -f "$1" ]; then
    echo "El archivo $1 no existe."
    exit 1
fi

# Crea una copia de respaldo del archivo original
cp "$1" "$1.bak"

# Agrega "chr" al principio de cada línea y guarda el resultado en un nuevo archivo
awk '{print "chr"$0}' "$1.bak" > "$1"

echo "Se ha agregado 'chr' al principio de cada línea en $1."
