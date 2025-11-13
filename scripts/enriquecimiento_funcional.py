import os
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt

# Configuración

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")

os.makedirs(RESULTS_DIR, exist_ok=True)

GENES_SEMILLA = os.path.join(RESULTS_DIR, "connected_seed_genes.tsv")
GENES_DIAMOND = os.path.join(RESULTS_DIR, "diamond_results.tsv")

# Funciones auxiliares

# Carga la lista de genes y devuelve los nombres únicos de genes en una lista
def cargar_genes(file_path, columna=None):
    """Carga una lista de genes desde archivo .txt o .tsv"""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"No se encontró el archivo: {file_path}")
    
    if file_path.endswith(".txt"):
        with open(file_path) as f:
            genes = [line.strip() for line in f if line.strip()]
    elif file_path.endswith(".tsv"):
        df = pd.read_csv(file_path, sep="\t")
        if columna is None:
            columna = df.columns[0]
        genes = df[columna].dropna().unique().tolist()
    else:
        raise ValueError("Formato de archivo no reconocido (usa .txt o .tsv)")

    print(f"Genes cargados ({len(genes)}): {file_path}")
    return genes

# Enriquecimiento funcional con GSEAPY (GOKEGG) para una lista de genes
def realizar_enriquecimiento(genes, conjunto_nombre, libreria="KEGG_2021_Human", outdir=RESULTS_DIR):
    """Ejecuta análisis de enriquecimiento con GSEAPY"""
    enr = gp.enrichr(
        gene_list=genes,
        gene_sets=[libreria, "GO_Biological_Process_2021"],
        organism="Human",
        outdir=None,  # No guardar automáticamente
        cutoff=0.05
    )

    # Guardar resultados
    resultado_path = os.path.join(outdir, f"enriquecimiento_{conjunto_nombre}.tsv")
    enr.results.to_csv(resultado_path, sep="\t", index=False)
    print(f"Resultados guardados en: {resultado_path}")

    return enr.results

# Miramos qué terminos se repiten entre dos listas y mostramos los comunes y los únicos
def comparar_enriquecimientos(df1, df2):
    """Compara los términos enriquecidos entre dos conjuntos"""
    set1 = set(df1["Term"].head(20))
    set2 = set(df2["Term"].head(20))
    
    comun = set1 & set2
    unicos1 = set1 - set2
    unicos2 = set2 - set1

    print("\n--- Comparación de Términos Enriquecidos ---")
    print(f"Términos comunes: {len(comun)}")
    print(f"Términos únicos en Semillas: {len(unicos1)}")
    print(f"Términos únicos en Candidatos: {len(unicos2)}")

    return comun, unicos1, unicos2

# Gráfico con los términos más significativos 
import numpy as np

def graficar_top_terms(df, conjunto_nombre, top_n=10):
    """Crea gráfico de barras de términos más significativos (robusto a nombres de columnas)."""
    if "Adjusted P-value" in df.columns:
        df = df.copy()
        df["__score__"] = -np.log10(np.clip(df["Adjusted P-value"].astype(float), 1e-300, None))
        score_label = "-log10(FDR)"
        sort_key = "Adjusted P-value"   # menor es mejor
        ascending = True
    elif "P-value" in df.columns:
        df = df.copy()
        df["__score__"] = -np.log10(np.clip(df["P-value"].astype(float), 1e-300, None))
        score_label = "-log10(p-valor)"
        sort_key = "P-value"
        ascending = True
    elif "Combined Score" in df.columns:
        df = df.copy()
        df["__score__"] = df["Combined Score"].astype(float)
        score_label = "Combined Score"
        sort_key = "Combined Score"
        ascending = False
    else:
        raise KeyError("No se encuentran columnas esperadas: 'Adjusted P-value', 'P-value' o 'Combined Score'.")

    df_top = df.sort_values(sort_key, ascending=ascending).head(top_n)

    plt.figure(figsize=(8, 6))
    plt.barh(df_top["Term"], df_top["__score__"])
    plt.xlabel(score_label)
    plt.title(f"Top {top_n} términos enriquecidos: {conjunto_nombre}")
    plt.gca().invert_yaxis()

    output_file = os.path.join(RESULTS_DIR, f"top_terms_{conjunto_nombre}.png")
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    print(f"Gráfico guardado en: {output_file}")



# Main

def main():
    print("\n=== Enriquecimiento Funcional ===")
    
    # 1. Cargar genes
    genes_semilla = cargar_genes(GENES_SEMILLA)
    genes_diamond = cargar_genes(GENES_DIAMOND, columna="HUGO_Symbol")
    
    # 2. Enriquecimiento para cada grupo
    enr_semilla = realizar_enriquecimiento(genes_semilla, "semillas")
    enr_diamond = realizar_enriquecimiento(genes_diamond, "candidatos")

    # 3. Comparar resultados
    comun, unicos_semilla, unicos_diamond = comparar_enriquecimientos(enr_semilla, enr_diamond)

    # 4. Visualizar
    graficar_top_terms(enr_semilla, "Semillas")
    graficar_top_terms(enr_diamond, "Candidatos")

    # 5. Guardar resumen comparativo
    resumen_path = os.path.join(RESULTS_DIR, "comparacion_enriquecimiento.txt")
    with open(resumen_path, "w") as f:
        f.write("=== Comparación de Enriquecimiento ===\n\n")
        f.write(f"Términos comunes ({len(comun)}):\n" + "\n".join(comun) + "\n\n")
        f.write(f"Términos únicos en Semillas ({len(unicos_semilla)}):\n" + "\n".join(unicos_semilla) + "\n\n")
        f.write(f"Términos únicos en Candidatos ({len(unicos_diamond)}):\n" + "\n".join(unicos_diamond) + "\n")
    
    print(f"\nResumen comparativo guardado en: {resumen_path}")
    print("=== Enriquecimiento completado ===")

if __name__ == "__main__":
    main()
