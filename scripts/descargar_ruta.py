import gseapy as gp
import os
import argparse
import sys

def main():
    # 1. Configuración de Argumentos
    parser = argparse.ArgumentParser(
        description="Descarga genes de una ruta específica de una librería de GSEAPY.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '--library', 
        type=str, 
        required=True, 
        help="Librería de rutas a usar (ej: KEGG_2021_Human, Reactome_2022)."
    )
    parser.add_argument(
        '--pathway', 
        type=str, 
        required=True, 
        help="Nombre o parte del nombre de la ruta a buscar (ej: autophagy)."
    )
    parser.add_argument(
        '--output-file', 
        type=str, 
        default="genes_ruta.txt", 
        help="Nombre del archivo de salida en la carpeta 'data/'."
    )
    args = parser.parse_args()

    LIBRARY = args.library
    BUSCADA = args.pathway
    OUTPUT_FILENAME = args.output_file

    print(f"--- 1. INICIANDO DESCARGA DE RUTA ---")
    print(f"Librería: {LIBRARY}, Ruta Buscada: {BUSCADA}")
    
    # 2. Descargar todas las rutas de esa librería
    try:
        genesets = gp.get_library(name=LIBRARY, organism="Human")
    except Exception as e:
        print(f"ERROR: No se pudo descargar la librería {LIBRARY}. {e}")
        sys.exit(1)


    ruta_encontrada = None
    genes_ruta = None

    for nombre, genes in genesets.items():
        if BUSCADA.lower() in nombre.lower():
            ruta_encontrada = nombre
            genes_ruta = genes
            break

    if ruta_encontrada is None:
        print(f"ERROR: No se encontró ninguna ruta en {LIBRARY} que contenga: {BUSCADA}")
        sys.exit(1)

    print(f"Ruta encontrada: {ruta_encontrada}")
    print(f"Número de genes: {len(genes_ruta)}")

    # 3. Guardar a un txt (un gen por línea)
    # Se utiliza la lógica de ruta relativa del usuario (asume que se ejecuta desde una carpeta 'scripts')
    try:
        script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
    except IndexError:
        script_dir = os.path.dirname(os.path.abspath(__file__))

    basedir = os.path.dirname(script_dir)
    datadir = os.path.join(basedir, "data")

    # Crear la carpeta si no existe
    os.makedirs(datadir, exist_ok=True)

    output_path = os.path.join(datadir, OUTPUT_FILENAME)

    with open(output_path, "w") as f:
        for g in genes_ruta:
            f.write(g + "\n")

    print(f"Genes guardados en: {output_path}")

if __name__ == '__main__':
    main()