#!/bin/bash

# =======================================================
# CONFIGURACIÓN DE PARÁMETROS GLOBALES
# =======================================================

# 1. FORZAR UTF-8 para evitar errores de codificación de emojis en Windows/Git Bash
export PYTHONIOENCODING=utf8

# 2. Determinar el intérprete de Python (PYTHON_EXEC)
# Intenta usar 'python3' (más explícito para Python 3), si no existe, usa 'python'.
if command -v python3 &> /dev/null; then
    PYTHON_EXEC="python3"
elif command -v python &> /dev/null; then
    PYTHON_EXEC="python"
else
    echo "ERROR: No se encontró el comando 'python' ni 'python3'. Asegúrate de que Python esté instalado y en tu PATH."
    exit 1
fi


# --- PARÁMETROS DE ENTRADA Y SALIDA ---
# PASO 1: Descarga de Genes Semilla
GSEAPY_LIBRARY="KEGG_2021_Human"
PATHWAY_SEARCH_TERM="autophagy"
# El archivo de entrada para DIAMOnD. Se crea en el Paso 1.
SEED_GENE_FILE="data/genes_ruta.txt" 

# PASO 2: Propagación DIAMOnD
# Archivo de red de interacciones (DEBE EXISTIR)
NETWORK_FILE="data/string_network_filtered_hugo-400.tsv" 
# Archivos de salida de DIAMOnD
DIAMOND_OUTPUT_FILE="results/diamond_results.tsv"
CONNECTED_SEEDS_FILE="results/connected_seed_genes.tsv"
DIAMOND_PLOT_FILE="diamond_network.png"

# =======================================================
# CONFIGURACIÓN DEL ENTORNO Y EJECUCIÓN
# =======================================================

echo "====================================================="
echo "INICIO DEL FLUJO DE TRABAJO AUTOMATIZADO"
echo "Librería: $GSEAPY_LIBRARY, Ruta: $PATHWAY_SEARCH_TERM"
echo "====================================================="

# Verificación de archivos esenciales
if [ ! -f "$NETWORK_FILE" ]; then
    echo "ERROR: El archivo de red '$NETWORK_FILE' no se encontró."
    echo "Asegúrate de que la carpeta 'network/' y el archivo de red existan."
    exit 1
fi

# -----------------------------------------------------------
# PASO 1: DESCARGA DE GENES SEMILLA (descargar_ruta.py)
# Output: data/genes_ruta.txt
# -----------------------------------------------------------
echo ""
echo "--- 1. EJECUTANDO DESCARGA DE RUTA ---"
# basename $SEED_GENE_FILE extrae 'genes_ruta.txt' para pasarlo como argumento
$PYTHON_EXEC scripts/descargar_ruta.py \
    --library "$GSEAPY_LIBRARY" \
    --pathway "$PATHWAY_SEARCH_TERM" \
    --output-file "$(basename $SEED_GENE_FILE)"

# $? es el código de salida del último comando ejecutado
if [ $? -ne 0 ]; then
    echo "ERROR (Paso 1): La descarga de la ruta falló. Abortando el flujo."
    exit 1
fi
echo "Paso 1 completado. Archivo de semillas: $SEED_GENE_FILE"


# -----------------------------------------------------------
# PASO 2: PROPAGACIÓN DIAMOnD (propagacion_diamond.py)
# Input: data/genes_ruta.txt y network/network.tsv
# Output: results/diamond_results.tsv y results/connected_seed_genes.tsv
# -----------------------------------------------------------
echo ""
echo "--- 2. EJECUTANDO PROPAGACIÓN DIAMOnD ---"

# Creamos la carpeta de resultados si no existe (aunque los scripts Python ya lo hacen, es buena práctica)
mkdir -p results

$PYTHON_EXEC scripts/propagacion_diamond.py \
    --seed-file "$SEED_GENE_FILE" \
    --input "$NETWORK_FILE" \
    --output "$(basename $DIAMOND_OUTPUT_FILE)" \
    --plot "$DIAMOND_PLOT_FILE"

if [ $? -ne 0 ]; then
    echo "ERROR (Paso 2): La propagación DIAMOnD falló. Abortando el flujo."
    exit 1
fi
echo "Paso 2 completado. Resultados DIAMOnD: $DIAMOND_OUTPUT_FILE"


# -----------------------------------------------------------
# PASO 3: ANÁLISIS DE RESULTADOS (analisis.py)
# Input: Archivos de resultados de DIAMOnD y el archivo de red original
# -----------------------------------------------------------
echo ""
echo "--- 3. EJECUTANDO SCRIPT DE ANÁLISIS ---"
$PYTHON_EXEC scripts/enriquecimiento_funcional.py \
    --connected-seeds "$CONNECTED_SEEDS_FILE" \
    --diamond-results "$DIAMOND_OUTPUT_FILE" \
    --network-file "$NETWORK_FILE"

if [ $? -ne 0 ]; then
    echo "ERROR (Paso 3): El script de análisis falló."
    exit 1
fi

echo "====================================================="
echo "FLUJO DE TRABAJO COMPLETADO CON ÉXITO"
echo "Resultados finales en la carpeta 'results/'"
echo "====================================================="