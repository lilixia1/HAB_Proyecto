# ----------------------------------------------------------------------
#          SCRIPT DE DESCARGA, FILTRADO Y MAPEO DE RED PPI (STRING DB)
# ----------------------------------------------------------------------
import requests
import os
import argparse
import sys
import pandas as pd
from tqdm import tqdm
import re

# --- Parámetros Globales de STRING DB ---
# ID de taxón para humano
ORGANISM_ID = 9606
# Versión de la base de datos (v12.0)
STRING_VERSION = "v12.0"
# URL base para la descarga de interacciones
STRING_LINKS_URL = f"https://stringdb-downloads.org/download/protein.links.{STRING_VERSION}/"
# URL base para la descarga de alias (mapeo a HUGO)
STRING_ALIAS_URL = f"https://stringdb-downloads.org/download/protein.aliases.{STRING_VERSION}/"
# Umbral de score combinado por defecto para el FILTRADO
SCORE_THRESHOLD = 700 

def download_and_process_aliases(organism: int) -> pd.Series:
    """
    Descarga el archivo de alias de STRING y crea un mapa de ID de STRING a HUGO Symbol.
    """
    print("\n--- Descargando Archivo de Alias (Mapeo a HUGO) ---")
    alias_filename = f"{organism}.protein.aliases.{STRING_VERSION}.txt.gz"
    download_url = STRING_ALIAS_URL + alias_filename
    
    try:
        # Leer el archivo de alias .gz directamente en un DataFrame
        alias_df = pd.read_csv(
            download_url, 
            compression='gzip', 
            sep='\t', 
            comment='#',
            header=None,
            usecols=[0, 1, 2], # STRING ID, Alias, Source
            names=['string_id', 'alias', 'source']
        )
    except Exception as e:
        print(f"[ERROR] Falló la descarga o lectura de alias: {e}")
        return pd.Series(dtype=str)

    # Filtrar para obtener solo los símbolos oficiales (HUGO)
    # Normalmente 'source' contiene 'Ensembl_gene', 'HUGO', 'UniProt', etc.
    # Buscamos específicamente los nombres de genes.
    hugo_map_df = alias_df[alias_df['source'].str.contains('HUGO|Ensembl|Gene_Name', case=False, na=False)]
    
    # Crear el mapa: {string_id: alias (HUGO Symbol)}
    # Usamos .drop_duplicates(subset='string_id', keep='first') para evitar conflictos, 
    # aunque el mapeo 1:1 no es perfecto en STRING, este es el mejor intento.
    hugo_map = hugo_map_df.drop_duplicates(subset='string_id', keep='first').set_index('string_id')['alias']
    
    print(f"Mapeo HUGO/Alias cargado: {len(hugo_map)} IDs únicos de STRING.")
    return hugo_map


def download_and_filter_string_network(organism: int, score_threshold: int, output_file: str, hugo_map: pd.Series):
    """
    Descarga las interacciones PPI, las filtra y mapea los IDs de STRING a HUGO.
    """
    links_filename = f"{organism}.protein.links.{STRING_VERSION}.txt.gz"
    download_url = STRING_LINKS_URL + links_filename
    
    print(f"\n--- Iniciando Descarga y Procesamiento de Interacciones PPI ---")
    print(f"URL de interacciones: {download_url}")
    print(f"Umbral de score combinado: >= {score_threshold}")
    
    try:
        # 1. Descargar y leer interacciones
        interactions_df = pd.read_csv(
            download_url, 
            compression='gzip', 
            sep=' ', 
            comment='#', 
            usecols=[0, 1, 2],
            names=['protein1', 'protein2', 'combined_score']
        )
    except Exception as e:
        print(f"[ERROR] Falló la lectura del archivo .gz de interacciones: {e}")
        sys.exit(1)

    original_count = len(interactions_df)

    # 2. Filtrar por el umbral de score
    interactions_df['combined_score'] = pd.to_numeric(interactions_df['combined_score'], errors='coerce')
    interactions_df = interactions_df.dropna(subset=['combined_score'])
    
    interactions_df_filtered = interactions_df[interactions_df['combined_score'] >= score_threshold]
    filtered_count_score = len(interactions_df_filtered)
    print(f"   Líneas originales: {original_count}")
    print(f"   Líneas después de filtro de score (>= {score_threshold}): {filtered_count_score}")


    # 3. Mapear IDs de STRING a HUGO Symbols
    # Los IDs de STRING tienen el prefijo '9606.' que hay que eliminar para el mapeo
    interactions_df_filtered['protein1'] = interactions_df_filtered['protein1'].apply(lambda x: re.sub(r'^\d+\.', '', str(x)))
    interactions_df_filtered['protein2'] = interactions_df_filtered['protein2'].apply(lambda x: re.sub(r'^\d+\.', '', str(x)))

    # Aplicar el mapeo
    print("Mapeando IDs de STRING a símbolos HUGO...")
    interactions_df_filtered['protein1'] = interactions_df_filtered['protein1'].map(hugo_map)
    interactions_df_filtered['protein2'] = interactions_df_filtered['protein2'].map(hugo_map)
    
    # 4. Eliminar interacciones que no pudieron mapearse a HUGO
    mapped_count = len(interactions_df_filtered)
    interactions_df_filtered = interactions_df_filtered.dropna(subset=['protein1', 'protein2'])
    final_count = len(interactions_df_filtered)
    
    print(f"   Líneas mapeadas a HUGO (finales): {final_count} (se eliminaron {mapped_count - final_count} interacciones sin mapeo).")

    # 5. Guardar el resultado filtrado y mapeado
    try:
        # Crear la carpeta de salida si no existe
        output_dir = os.path.dirname(output_file)
        os.makedirs(output_dir, exist_ok=True)
        
        # Guardar el DataFrame
        interactions_df_filtered.to_csv(
            output_file, 
            sep="\t", 
            index=False, 
            header=True # Mantenemos la cabecera (protein1, protein2, combined_score)
        )
        print(f"\nRed de interacciones (HUGO, Score >= {score_threshold}) guardada en: {output_file}")
    except Exception as e:
        print(f"[ERROR] No se pudo guardar el archivo en {output_file}. Error: {e}")
        sys.exit(1)


def main():
    
    parser = argparse.ArgumentParser(
        description="Descarga, filtra y mapea a HUGO la red PPI de STRING DB."
    )
    parser.add_argument(
        '--organism', 
        type=int, 
        default=ORGANISM_ID, 
        help=f"ID de taxón (e.g., 9606 para Humano) (default: {ORGANISM_ID})."
    )
    parser.add_argument(
        '--score', 
        type=int, 
        default=SCORE_THRESHOLD, 
        help=f"Umbral mínimo de Combined Score [0-1000] (default: {SCORE_THRESHOLD})."
    )
    parser.add_argument(
        '--output-file', 
        type=str, 
        required=True, # Ahora es obligatorio
        help="Ruta y nombre del archivo de salida para la red filtrada (ej: data/network.tsv)."
    )
    args = parser.parse_args()

    # 1. Obtener el mapa de IDs de STRING a HUGO
    hugo_map = download_and_process_aliases(args.organism)
    if hugo_map.empty:
        print("Error: No se pudo cargar el mapeo de alias. Abortando.")
        sys.exit(1)

    # 2. Descargar, filtrar y mapear la red PPI
    download_and_filter_string_network(args.organism, args.score, args.output_file, hugo_map)


if __name__ == '__main__':
    main()