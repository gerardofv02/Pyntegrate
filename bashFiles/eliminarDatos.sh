#!/bin/bash

input_file="completo.gtf"
output_file="gencode.vM25.annotation.gtf"

# Verificar que el archivo de entrada existe
if [[ ! -f "$input_file" ]]; then
    echo "El archivo $input_file no existe."
    exit 1
fi

# Crear o limpiar el archivo de salida
> "$output_file"

# Crear un diccionario para almacenar los conteos por tipo
declare -A type_count

# Leer el archivo línea por línea
while IFS=$'\t' read -r -a fields; do
    # Verificar si la línea está vacía
    if [[ -z "${fields[0]}" ]]; then
        continue
    fi

    # Obtener los primeros 4 caracteres del primer campo (columna del cromosoma/contig)
    type=${fields[0]:0:5}
    echo "Tipo detectado: $type"

    # Incrementar el contador para este tipo si no existe
    if [[ -z "${type_count[$type]}" ]]; then
        type_count[$type]=0
    fi

    # Si hay menos de 5 filas para este tipo, agregar la línea al archivo de salida
    if [[ ${type_count[$type]} -lt 50000 ]]; then
        echo -e "${fields[*]}" >> "$output_file"
        echo "Añadiendo línea al archivo de salida: ${fields[*]}"
        ((type_count[$type]++))
    fi
done < "$input_file"

echo "Proceso completado. Los datos filtrados se guardaron en $output_file."
