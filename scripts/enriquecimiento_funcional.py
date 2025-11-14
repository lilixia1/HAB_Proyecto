import os
import re
import numpy as np
import pandas as pd
import gseapy as gp
import networkx as nx
import matplotlib.pyplot as plt
import argparse 
import sys


try:
    import community as community_louvain
    HAS_LOUVAIN = True
except Exception:
    community_louvain = None
    HAS_LOUVAIN = False

# Configuración base (Se mantiene el cálculo de directorios)
try:
    script_path = os.path.abspath(sys.argv[0])
except IndexError:
    script_path = os.path.abspath(__file__)

BASE_DIR = os.path.dirname(os.path.dirname(script_path))
DATA_DIR = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

PPI_SCORE_UMBRAL = 700

# Funciones auxiliares

# Carga la lista de genes y devuelve los nombres únicos de genes en una lista
def cargar_genes(file_path, columna=None):
    """Carga una lista de genes desde archivo .txt o .tsv"""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"No se encontró el archivo: {file_path}")
    
    if file_path.endswith(".txt"):
        with open(file_path, encoding="utf-8-sig") as f:
            genes = [l.strip().upper() for l in f if l.strip()]
    elif file_path.endswith(".tsv"):
        df = pd.read_csv(file_path, sep="\t")
        col = columna or df.columns[0]
        genes = df[col].dropna().astype(str).str.upper().tolist()
    else:
        raise ValueError("Formato de archivo no reconocido (usa .txt o .tsv)")

    # Normalización mínima 
    clean = []
    for g in genes:
        g1 = (g.replace("Ø", "0")
                .replace("—", "-").replace("–", "-").replace("−", "-")
                .strip())
        clean.append(g1)
    genes = list(dict.fromkeys(clean))  # únicos preservando orden

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
def graficar_top_terms(df, conjunto_nombre, top_n=10):
    """Crea gráfico de barras de términos más significativos (robusto a nombres de columnas)."""
    
    df = df.copy()

    if "Adjusted P-value" in df.columns:
        df["__score__"] = -np.log10(np.clip(df["Adjusted P-value"].astype(float), 1e-300, None))
        score_label = "-log10(FDR)"
        sort_key = "Adjusted P-value"   # menor es mejor
        ascending = True
    elif "P-value" in df.columns:
        df["__score__"] = -np.log10(np.clip(df["P-value"].astype(float), 1e-300, None))
        score_label = "-log10(p-valor)"
        sort_key = "P-value"
        ascending = True
    elif "Combined Score" in df.columns:
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

def construir_grafo(ppi_file, umbral=PPI_SCORE_UMBRAL):
    """Crea un grafo NetworkX desde un TSV PPI (filtra por score si existe)."""
    df = pd.read_csv(ppi_file, sep="\t")
    
    if umbral is not None:
        df = df[df["combined_score"] >= umbral]

    G = nx.Graph()
    for _, row in df.iterrows():
        G.add_edge(str(row["protein1_hugo"]).upper(),
                   str(row["protein2_hugo"]).upper(),
                   weight=row["combined_score"])
    return G

def calcular_propiedades(G, semillas, candidatos):
    """Calcula grado, centralidades y modularidad en la subred semillas∪candidatos."""
    sub_nodes = (set(semillas) | set(candidatos)) & set(G.nodes())
    subgraph = G.subgraph(sub_nodes).copy()

    # Métricas en la subred
    grado = dict(subgraph.degree())
    deg_cent = nx.degree_centrality(subgraph)
    betw = nx.betweenness_centrality(subgraph, normalized=True)

    df = pd.DataFrame({
        "Gen": list(sub_nodes),
        "Grado": [grado[n] for n in sub_nodes],
        "Centralidad": [deg_cent[n] for n in sub_nodes],
        "Betweenness": [betw[n] for n in sub_nodes],
        "Tipo": ["Semilla" if n in semillas else "Candidato" for n in sub_nodes]
    }).sort_values(["Grado", "Centralidad"], ascending=False)

    modularidad = None
    if HAS_LOUVAIN and subgraph.number_of_edges() > 0:
        part = community_louvain.best_partition(subgraph, random_state=42)
        modularidad = community_louvain.modularity(part, subgraph)

    return df, modularidad, subgraph

def graficar_metrica(df, metrica, title, top_n=20):
    df_plot = df.sort_values(metrica, ascending=False).head(top_n)
    colores = df_plot["Tipo"].map({"Semilla": "tab:blue", "Candidato": "tab:orange"})
    plt.figure(figsize=(10, 5))
    plt.bar(df_plot["Gen"], df_plot[metrica], color=colores)
    plt.xticks(rotation=90, fontsize=8)
    plt.ylabel(metrica); plt.title(title)
    plt.tight_layout()
    out = os.path.join(RESULTS_DIR, f"{re.sub(r'\\W+', '_', title.lower())}.png")
    plt.savefig(out); plt.close()
    print(f"Gráfico guardado en: {out}")

# Main


def main():
    # --- 0. CONFIGURACIÓN DE ARGUMENTOS ---
    parser = argparse.ArgumentParser(
        description="Realiza el análisis estructural y de enriquecimiento funcional de los resultados DIAMOnD."
    )
    parser.add_argument(
        '--connected-seeds', 
        required=True, 
        help="Ruta al archivo de genes semilla conectados (connected_seed_genes.tsv)."
    )
    parser.add_argument(
        '--diamond-results', 
        required=True, 
        help="Ruta al archivo de resultados de DIAMOnD (diamond_results.tsv)."
    )
    parser.add_argument(
        '--network-file', 
        required=True, 
        help="Ruta al archivo de interacciones PPI original (p. ej., string_network_filtered.tsv)."
    )
    args = parser.parse_args()
    
    genes_semilla_path = args.connected_seeds
    genes_diamond_path = args.diamond_results
    interacciones_file_path = args.network_file
    
    print("\n=== Enriquecimiento Funcional + Análisis Estructural ===")
    
    # 1. Cargar genes
    try:
        genes_semilla = cargar_genes(genes_semilla_path)
        genes_diamond = cargar_genes(genes_diamond_path)
    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        sys.exit(1)
    
    if not genes_semilla and not genes_diamond:
        print("No hay genes semilla conectados ni genes DIAMOnD. Finalizando análisis.")
        return

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

    # 6. Análisis estructural
    print("\n--- Construyendo subred PPI ---")
    G = construir_grafo(interacciones_file_path, umbral=PPI_SCORE_UMBRAL)
    print(f"Nodos totales en la red: {G.number_of_nodes()}")
        
    df_struct, modularidad, subgraph = calcular_propiedades(G, genes_semilla, genes_diamond)
    
    if df_struct.empty:
        print("No se pudo realizar el análisis estructural (subred vacía).")
        return

    # Guardar métricas y subred
    estructural_path = os.path.join(RESULTS_DIR, "analisis_estructural.tsv")
    df_struct.to_csv(estructural_path, sep="\t", index=False)
    print(f"Resultados estructurales guardados en: {estructural_path}")
    if modularidad is not None:
        print(f"Modularidad de la subred (Louvain): {modularidad:.3f}")
    else:
        print("Modularidad no calculada (instala 'python-louvain' para obtenerla).")

    # Gráficos de métricas
    graficar_metrica(df_struct, "Grado", "Top 20 genes por grado en subred")
    graficar_metrica(df_struct, "Centralidad", "Top 20 genes por centralidad de grado en subred")
    graficar_metrica(df_struct, "Betweenness", "Top 20 genes por Betweenness en subred")

    # Exportar subred a GraphML (Cytoscape/Gephi)
    sub_path = os.path.join(RESULTS_DIR, "subred_enriquecida.graphml")
    nx.write_graphml(subgraph, sub_path)
    print(f"Subred exportada a: {sub_path}")

    print("=== Análisis completado ===")

if __name__ == "__main__":
    main()