#!/bin/bash

# Ruta al directorio donde están almacenados los archivos de referencia indexados
genome_dir="./Genome"

# Ruta al directorio donde se encuentran los directorios de muestras
samples_dir="./Samples"

# Ruta al directorio donde se almacenarán los resultados
results_dir="./Results"

# Lista de directorios de muestra específicos
specific_sample_dirs=("H25_RNa" "SOS1_HC" "SOS1_HNa" "SOS1_RC" "SOS1_RNa")

# Recorremos cada directorio de muestra específico
for sample_name in "${specific_sample_dirs[@]}"; do
    # Ruta al directorio de muestra
    sample_dir="$samples_dir/Trimmed_Data/$sample_name"
    
    if [ -d "$sample_dir" ]; then
        # Creamos el directorio correspondiente en "Results"
        mkdir -p "$results_dir/$sample_name"

        # Ejecutamos HISAT2 para mapear las lecturas
        hisat2 --dta -x "$genome_dir/index" -1 "$sample_dir"/*R1.fastq.gz -2 "$sample_dir"/*R2.fastq.gz -S "$results_dir/$sample_name/$sample_name.sam"

        # Convertimos el archivo SAM a BAM y ordenamos
        samtools view -bS "$results_dir/$sample_name/$sample_name.sam" | samtools sort -o "$results_dir/$sample_name/$sample_name.bam"

        # Creamos el índice para el archivo BAM
        samtools index "$results_dir/$sample_name/$sample_name.bam"

        # Eliminamos el archivo SAM para liberar espacio en disco
        rm "$results_dir/$sample_name/$sample_name.sam"
        
        echo "Análisis completado para $sample_name"
    else
        echo "Directorio de muestra $sample_name no encontrado."
    fi
done

echo "¡Todos los análisis completados!"
