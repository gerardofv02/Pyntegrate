#!/bin/bash

# Verificar si se proporcionó el archivo de entrada
if [ $# -ne 1 ]; then
    echo "Uso: $0 <archivo_de_entrada>"
    exit 1
fi

archivo_entrada="$1"

# Verificar si el archivo de entrada existe
if [ ! -f "$archivo_entrada" ]; then
    echo "El archivo $archivo_entrada no existe."
    exit 1
fi

# Archivo temporal para almacenar líneas únicas
archivo_temporal=$(mktemp)

# Extraer gene_name del archivo de entrada y guardar líneas únicas en archivo temporal
while IFS=$'\t' read -r -a campos; do
    # Buscar gene_name dentro del campo gene_id o gene_name
    for campo in "${campos[@]}"; do
        if [[ $campo =~ gene_id\ \"([^\"]+)\" ]]; then
            gene_name="${BASH_REMATCH[1]}"
            break
        elif [[ $campo =~ gene_name\ \"([^\"]+)\" ]]; then
            gene_name="${BASH_REMATCH[1]}"
            break
        fi
    done
    
    # Verificar si ya se ha visto este gene_name
    if ! grep -q "^$gene_name$" "$archivo_temporal"; then
        echo "$gene_name" >> "$archivo_temporal"
        echo -e "${campos[@]}" >> tsses.gtf
    fi
done < "$archivo_entrada"

# Limpiar archivo temporal
rm "$archivo_temporal"

echo "Archivo procesado. Las líneas únicas basadas en gene_name se han guardado en tsses_sin_duplicados.gtf"
