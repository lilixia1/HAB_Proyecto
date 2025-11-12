# ----------------------------------------------------------------------
#                          LIBRERIAS Y SETUP
# ----------------------------------------------------------------------
import os
import argparse
import pandas as pd
from tqdm import tqdm
import networkx as nx
import scipy.stats
import matplotlib.pyplot as plt
import numpy as np
import scipy.special
import operator

# --- Parámetros Globales ---
nodos_añadidos = 200
especie = 'human' # Se mantiene solo por consistencia, pero no se usa para mapeo
UMBRAL_SCORE = 700

# ----------------------------------------------------------------------
#                   FUNCIONES DE LECTURA Y CONVERSIÓN
# ----------------------------------------------------------------------
def importar_genes(file_path: str) -> list[str]:
    """
    Carga una lista de genes desde un archivo donde hay un gen por linea
    """
    print(f"\nCargando genes desde: {file_path}")
    
    if not os.path.exists(file_path):
        print(f"[ERROR] Archivo no encontrado: {file_path}")
        return []
        
    try:
        with open(file_path, 'r') as f:
            lista_genes = []
            for linea in f:
                lista_genes.append(linea.strip())
                
    except Exception as e:
        print(f"[ERROR] Error al leer el archivo {file_path}: {e}")
        return []
    
    print(f"Genes cargados ({len(lista_genes)}): {lista_genes}")

    return lista_genes


def cargar_red_conocida(fichero, umbral=700):
    """Carga las interacciones de un archivo TSV y nombra las columnas como HUGO."""
    try:
        interactions = pd.read_csv(fichero, sep="\t")
        interactions.columns = [
            "protein1_hugo", "protein2_hugo", "combined_score"
        ]

        # FILTRADO CRÍTICO: Eliminar interacciones por debajo del umbral, sino la propagacion sera demasiado extensa
        lineas_originales = len(interactions)
        interactions = interactions[interactions['combined_score'] >= umbral]
        return interactions
    
    except Exception as e:
        print(f"Error al cargar los datos de {fichero}: {e}")
        return None

def construir_red(interactions_df, SEED_GENES_ENTREZ):
    """Construye un grafo NetworkX usando los HUGO ID y el score combinado como peso."""
    G = nx.Graph()
    for _, row in interactions_df.iterrows():
        # Asumiendo que 'combined_score' ya está entre 0 y 1000, lo dividimos por 1000.
        G.add_edge(row["protein1_hugo"], row["protein2_hugo"], weight=row["combined_score"] / 1000.0)

    seed_interactions = interactions_df[
        interactions_df["protein1_hugo"].isin(SEED_GENES_ENTREZ) | 
        interactions_df["protein2_hugo"].isin(SEED_GENES_ENTREZ)
    ]
    print(f"Interacciones encontradas para genes semilla: {seed_interactions.shape[0]}")
    if seed_interactions.empty:
        print("Advertencia: No se encontraron interacciones para los genes semilla.")
    return G

# ----------------------------------------------------------------------
#                         FUNCIONES DIAMOnD
# ----------------------------------------------------------------------

def compute_all_gamma_ln(N):
    gamma_ln = {i: scipy.special.gammaln(i) for i in range(1, N + 1)}
    return gamma_ln

def gauss_hypergeom(x, r, b, n, gamma_ln):
    max_index = len(gamma_ln) - 1
    if r + b > max_index or x < 0 or r < x or b < (n - x) or n < x:
        return 0 

    log_p = (gamma_ln[r] - (gamma_ln[x] + gamma_ln[r - x]) +
             gamma_ln[b] - (gamma_ln[n - x] + gamma_ln[b - n]) -
             gamma_ln[r + b])
             
    return np.exp(log_p)

def pvalue(kb, k, N, s, gamma_ln):
    """Calcula el p-valor hipergeométrico."""
    return sum(gauss_hypergeom(n, s, N - s, k, gamma_ln) for n in range(kb, k + 1))

def diamond_iteration_of_first_X_nodes(G, S, X, alpha=1):
    added_nodes = []
    neighbors = {node: set(G.neighbors(node)) for node in G.nodes}
    degrees = dict(G.degree())
    cluster_nodes = set(S)
    
    if len(G.nodes) == 0:
        return []

    gamma_ln = compute_all_gamma_ln(len(G.nodes))

    while len(added_nodes) < X:
        min_p = float('inf')
        next_node = None
        
        candidate_nodes = set()
        for seed in cluster_nodes:
            candidate_nodes.update(neighbors.get(seed, set()))
        candidate_nodes -= cluster_nodes

        if not candidate_nodes:
            break

        for node in tqdm(candidate_nodes, desc=f"DIAMOnD Iteración {len(added_nodes) + 1}/{X}"):
            k = degrees.get(node, 0)
            kb = sum((1 for neighbor in neighbors[node] if neighbor in cluster_nodes))

            if k == 0 or k < kb: 
                continue

            try:
                p = pvalue(kb, k, len(G.nodes), len(cluster_nodes), gamma_ln)
            except Exception:
                continue

            if p < min_p:
                min_p = p
                next_node = node

        if next_node:
            added_nodes.append(next_node)
            cluster_nodes.add(next_node)
        else:
            break
            
    return added_nodes

# ----------------------------------------------------------------------
#                         FUNCIONES DE GUARDADO
# ----------------------------------------------------------------------
def guardar_resultados(seed_genes_hugo, diamond_genes_hugo, archivo_salida):
    """
    Guarda solo los genes añadidos por DIAMOnD en un archivo TSV. 
    """
    
    # 1. Crear el conjunto de genes añadidos, excluyendo cualquier duplicado con las semillas
    # (Aunque diamond_genes ya debería excluir las semillas, es una doble verificación)
    seed_set = set(seed_genes_hugo)
    diamond_only_genes = [gene for gene in diamond_genes_hugo if gene not in seed_set]

    data = []
    for hugo_symbol in diamond_only_genes:
        data.append({
            'HUGO_Symbol': hugo_symbol
        })

    results_df = pd.DataFrame(data)
    results_df.to_csv(archivo_salida, sep="\t", index=False)
    print(f"\nResultados de DIAMOnD guardados en: {archivo_salida} ({len(results_df)} genes)")

def graficar_red_enriquecida(G, seed_genes_hugo, diamond_genes_hugo, output_image_file):
    
    all_nodes = set(seed_genes_hugo) | set(diamond_genes_hugo)
    # 1. Crear el subgrafo
    subgraph = G.subgraph(all_nodes)
    
    if len(subgraph.nodes) < 2:
        print("Advertencia: Grafo enriquecido demasiado pequeño para dibujar.")
        return

    # 2. Calcular posiciones (solo para los nodos del subgrafo)
    pos = nx.spring_layout(subgraph, k=0.15, iterations=20) 
    plt.figure(figsize=(12, 12))

    # 3. FILTRADO CRÍTICO: Asegurarse de que solo se dibuje lo que está en 'pos'
    
    # Filtrar Semillas
    seed_genes_presentes = [n for n in seed_genes_hugo if n in pos] 
    
    # Filtrar DIAMOnD (los que no son semilla)
    diamond_only = list((set(diamond_genes_hugo) - set(seed_genes_hugo)) & set(pos.keys()))
    
    # --- DIBUJADO ---

    # Nodos de semillas (Azul)
    nx.draw_networkx_nodes(subgraph, pos, nodelist=seed_genes_presentes, node_color='blue', node_size=800, label="Seed Genes", alpha=0.8)

    # Nodos DIAMOnD (Naranja)
    nx.draw_networkx_nodes(subgraph, pos, nodelist=diamond_only, node_color='orange', node_size=500, label="DIAMOnD Added", alpha=0.7)
    
    # ... (el resto del código sigue igual)
    
    # Etiquetas:
    node_labels = {n: n for n in subgraph.nodes()} 
    nx.draw_networkx_labels(subgraph, pos, labels=node_labels, font_size=8, font_weight='bold')
    plt.legend(loc="upper left", markerscale=0.7)
    plt.title(f"Red Enriquecida con DIAMOnD ({len(subgraph.nodes)} Nodos)")
    plt.axis('off')
    
    try:
        plt.savefig(output_image_file, format='png', bbox_inches='tight')
        print(f"Grafo de DIAMOnD guardado como imagen en: {output_image_file}")
    except Exception as e:
        print(f"Error al guardar la imagen: {e}")
        
    plt.close()


# ----------------------------------------------------------------------
#                             MAIN CLI
# ----------------------------------------------------------------------

def main():
    
    # 1. Configuración de argparse para CLI
    parser = argparse.ArgumentParser(
        description="Propagación de genes en redes de interacción usando DIAMOnD."
    )
    parser.add_argument(
        '--seed-file', 
        required=True, 
        help="Ruta al archivo de texto que contiene los genes semilla (HUGO)."
    )
    parser.add_argument(
        '--input', 
        required=True, 
        help="Ruta al archivo de entrada de la red (HUGO/TSV)."
    )
    parser.add_argument(
        '--output', 
        default='diamond_results.tsv', 
        help="Nombre del archivo de resultados TSV."
    )
    parser.add_argument(
        '--plot',
        default='diamond_network.png', 
        help="Nombre del archivo para la imagen de la red DIAMOnD."
    )
    args = parser.parse_args()
    
    print(f"---Iniciando Propagación DIAMOnD ---")
    
    # 2. Creación de la Carpeta de Resultados
    RESULTS_DIR = "results"
    if not os.path.exists(RESULTS_DIR):
        os.makedirs(RESULTS_DIR)
        print(f"📁 Carpeta de resultados creada: {RESULTS_DIR}")

    output_tsv_path = os.path.join(RESULTS_DIR, args.output)
    output_plot_path = os.path.join(RESULTS_DIR, args.plot)
    print(f"Archivo de entrada de la red: {args.input}")
    
    ## 3. CARGA DE DATOS
    genes_semilla_hugo = importar_genes(args.seed_file)
    
    if not genes_semilla_hugo:
        return
        
    # Carga y filtrado (usando el alto umbral definido)
    interacciones = cargar_red_conocida(args.input)
    if interacciones is None or interacciones.empty:
        return

    ## 4. CONSTRUIR GRAFO
    # Se usa el dataframe filtrado y los símbolos HUGO directamente
    red = construir_red(interacciones, genes_semilla_hugo)
    
    if len(red.nodes) == 0:
        print("El grafo está vacío. No se puede ejecutar la propagación.")
        return

    # Determinamos el número real de nodos a añadir
    n = min(nodos_añadidos, len(red.nodes) - len(genes_semilla_hugo))
    if n <= 0:
        print("No hay nodos para añadir o la red es muy pequeña.")
        return

    ## 5. EJECUTAR DIAMOnD
    print(f"\n---Ejecutando DIAMOnD para añadir {n} nodos ---")
    # S (semillas) y el resultado son HUGO IDs
    diamond_genes = diamond_iteration_of_first_X_nodes(red, genes_semilla_hugo, n)

    ## 6. GUARDAR Y GRAFICAR RESULTADOS
    guardar_resultados(genes_semilla_hugo, diamond_genes, output_tsv_path)
    graficar_red_enriquecida(red, genes_semilla_hugo, diamond_genes, output_plot_path)

if __name__ == '__main__':
    main()
