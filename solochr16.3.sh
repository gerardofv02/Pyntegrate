#!/bin/bash

# Nombre del archivo de entrada y salida
archivo="tsses.gtf"

# Verificar si el archivo existe
if [ ! -f "$archivo" ]; then
    echo "El archivo $archivo no existe."
    exit 1
fi

# Crear un archivo temporal para almacenar las líneas filtradas
archivo_temp=$(mktemp)

# Filtrar las líneas que comienzan con 'ch16' y guardarlas en el archivo temporal
grep '^chr16' "$archivo" > "$archivo_temp"

# Reemplazar el archivo original con el archivo filtrado
mv "$archivo_temp" "$archivo"

echo "Se han eliminado todas las líneas que no comienzan con 'ch16' del archivo $archivo."
