#!/bin/bash

# Definimos la ruta del directorio base
base_dir="Samples"

# Nos aseguramos de estar en el directorio correcto. Este script se ejecuta en el directorio principal, en el que se encontrará Samples.
cd "$base_dir" || exit

# Obtenemos la lista de subdirectorios dentro de "Raw_Data"
subdirs=$(find Raw_Data -type d)

# Creamos los mismos subdirectorios vacíos dentro de "Fastqc_Files_Raw"
for subdir in $subdirs; do
  mkdir -p "Fastqc_Files_Raw/${subdir#Raw_Data/}"
done

# Realizamos el análisis FastQC para cada archivo en "Raw_Data" y lo almacenamos en su directorio correspondiente en "Fastqc_Files_Raw"
for subdir in $subdirs; do
  files=$(find "$subdir" -maxdepth 1 -type f -printf '%P\n')
  for file in $files; do
    fastqc "$subdir/$file" -o "Fastqc_Files_Raw/${subdir#Raw_Data/}"
  done
done

# Ejecutar MultiQC para generar un informe combinado de todos los archivos FastQC
multiqc_output_dir="MultiQC_RAW"
multiqc Fastqc_Files_Raw -o "$multiqc_output_dir"