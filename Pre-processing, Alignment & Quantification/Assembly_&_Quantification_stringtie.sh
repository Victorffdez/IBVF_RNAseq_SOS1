#!/bin/bash

# Ruta de la anotación GTF
annotation="./Annotation/Oryza_sativa.IRGSP-1.0.57.gtf"

# Creación del archivo mergelist.txt
mergelist_file="./Results/mergelist.txt"

# Recorre cada directorio de condiciones
for condition_dir in ./Results/*/; do
    condition_name=$(basename "$condition_dir")
    echo "Procesando condición: $condition_name"

    # Ensamblado de transcritos
    stringtie -G "$annotation" -o "$condition_dir$condition_name.gtf" -l "$condition_name" "$condition_dir$condition_name.bam"

    # Agrega la ruta al archivo mergelist.txt
    echo "$condition_dir$condition_name.gtf" >> "$mergelist_file"
done

# Ensamblado del transcriptoma completo
stringtie --merge -G "$annotation" -o ./Results/stringtie_merged.gtf "$mergelist_file"

# Cuantificación de niveles de expresión
for condition_dir in ./Results/*/; do
    condition_name=$(basename "$condition_dir")
    echo "Cuantificando expresión para condición: $condition_name"

    stringtie -e -B -G ./Results/stringtie_merged.gtf -o "$condition_dir""quantified.gtf" "$condition_dir$condition_name.bam"
done

echo "Proceso finalizado."
