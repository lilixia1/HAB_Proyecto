# ----------------------------------------------------------------------
#                                                  LIBRERIAS Y SETUP
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
import sys # Importado para manejo de rutas

# --- Parámetros Globales ---
nodos_añadidos = 200
especie = 'human' 
UMBRAL_SCORE = 700

# ----------------------------------------------------------------------
#                                   FUNCIONES DE LECTURA Y CONVERSIÓN
# ----------------------------------------------------------------------
def importar_genes(file_path: str) -> list[str]:
        """
        Carga una lista de genes desde un archivo donde hay un gen por linea,
    eliminando duplicados y líneas vacías.
        """
        print(f"\nCargando genes desde: {file_path}")
        
        if not os.path.exists(file_path):
                print(f"[ERROR] Archivo no encontrado: {file_path}")
                return []
                
        try:
                with open(file_path, 'r') as f:
                        lista_genes = []
                        for linea in f:
                                cleaned_gene = linea.strip()
                                if cleaned_gene:
                                        lista_genes.append(cleaned_gene)
                        
                        # Eliminar duplicados manteniendo el orden de aparición
                        lista_genes = list(dict.fromkeys(lista_genes))

        except Exception as e:
                print(f"[ERROR] Error al leer el archivo {file_path}: {e}")
                return []
        
        print(f"Genes cargados ({len(lista_genes)}): {lista_genes}")

        return lista_genes


def cargar_red_conocida(fichero, umbral=UMBRAL_SCORE):
        """Carga las interacciones de un archivo TSV, filtra por score, y nombra las columnas."""
        print(f"\nCargando datos de {fichero} y aplicando umbral de score >= {umbral}...")
        try:
                interactions = pd.read_csv(fichero, sep="\t")
                
                # Se asegura de manejar solo las tres columnas relevantes
                if interactions.shape[1] >= 3:
                        interactions = interactions.iloc[:, [0, 1, 2]]
                
                interactions.columns = [
                        "protein1_hugo", "protein2_hugo", "combined_score"
                ]

                # Asegurar que el score es numérico antes de filtrar
                interactions['combined_score'] = pd.to_numeric(interactions['combined_score'], errors='coerce')
                interactions.dropna(subset=['combined_score'], inplace=True)
                
                # FILTRADO CRÍTICO
                lineas_originales = len(interactions)
                interactions = interactions[interactions['combined_score'] >= umbral]
                
                print(f"   Líneas originales: {lineas_originales}")
                print(f"   Líneas retenidas (score >= {umbral}): {len(interactions)}")
                
                return interactions
        
        except Exception as e:
                print(f"Error al cargar los datos de {fichero}: {e}")
                return None

def construir_red(interactions_df, SEED_GENES_HUGO):
        """
    Construye un grafo NetworkX y devuelve el grafo (G) y los genes semilla
    que están presentes en él (genes_semilla_valid).
    """
        G = nx.Graph()
        for _, row in interactions_df.iterrows():
                # Asumiendo que 'combined_score' ya está entre 0 y 1000, lo dividimos por 1000.
                G.add_edge(row["protein1_hugo"], row["protein2_hugo"], weight=row["combined_score"] / 1000.0)

        # Identificar qué genes semilla están realmente en el grafo (conectados a algo)
        genes_semilla_valid = [gene for gene in SEED_GENES_HUGO if gene in G]
        
        print(f"Genes semilla encontrados en la red: {len(genes_semilla_valid)}/{len(SEED_GENES_HUGO)}")
        
        return G, genes_semilla_valid

# ----------------------------------------------------------------------
#                                                 FUNCIONES DIAMOnD
# ----------------------------------------------------------------------

def compute_all_gamma_ln(N):
        gamma_ln = {i: scipy.special.gammaln(i) for i in range(1, N + 1)}
        return gamma_ln

def gauss_hypergeom(x, r, b, n, gamma_ln):
        max_index = len(gamma_ln) - 1
        # Comprobaciones de límites
        if r + b > max_index or r < 0 or b < 0 or n < 0 or x < 0 or r < x or b < (n - x) or n < x:
                return 0 

        try:
                # Cálculo del logaritmo de la probabilidad hipergeométrica
                log_p = (gamma_ln[r] - (gamma_ln[x] + gamma_ln[r - x]) + gamma_ln[b] - (gamma_ln[n - x] + gamma_ln[b - n]) - gamma_ln[r + b])
                return np.exp(log_p)
        except KeyError:
                return 0 
        except Exception:
                return 0


def pvalue(kb, k, N, s, gamma_ln):
        """Calcula el p-valor hipergeométrico."""
        sum_p = 0.0
        for n in range(kb, k + 1):
                sum_p += gauss_hypergeom(n, s, N - s, k, gamma_ln)
        return sum_p


def diamond_iteration_of_first_X_nodes(G, S_valid, X, alpha=1):
        """
    Ejecuta el algoritmo DIAMOnD usando SOLO los genes semilla válidos (S_valid)
    como cluster inicial.
    """
        added_nodes = []
        neighbors = {node: set(G.neighbors(node)) for node in G.nodes}
        degrees = dict(G.degree())
        cluster_nodes = set(S_valid) # USAMOS SOLO LAS SEMILLAS VÁLIDAS
        
        if len(G.nodes) == 0 or not S_valid:
                print("Grafo vacío o no hay genes semilla válidos para iniciar DIAMOnD.")
                return []

        gamma_ln = compute_all_gamma_ln(len(G.nodes))
        N = len(G.nodes)

        while len(added_nodes) < X:
                min_p = float('inf')
                next_node = None
                
                candidate_nodes = set()
                for seed in cluster_nodes:
                        candidate_nodes.update(neighbors.get(seed, set()))
                candidate_nodes -= cluster_nodes

                if not candidate_nodes:
                        print("Todos los nodos vecinos han sido añadidos. Deteniendo la propagación.")
                        break

                for node in tqdm(candidate_nodes, desc=f"DIAMOnD Iteración {len(added_nodes) + 1}/{X} (Cluster: {len(cluster_nodes)})"):
                        k = degrees.get(node, 0)
                        # Contar las conexiones al clúster (kb)
                        kb = sum((1 for neighbor in neighbors[node] if neighbor in cluster_nodes))

                        s = len(cluster_nodes) # Tamaño actual del cluster

                        if k == 0 or kb == 0: 
                                continue

                        try:
                                p = pvalue(kb, k, N, s, gamma_ln)
                        except Exception:
                                continue

                        if p < min_p:
                                min_p = p
                                next_node = node

                if next_node:
                        added_nodes.append(next_node)
                        cluster_nodes.add(next_node)
                else:
                        print(f"No se encontró ningún candidato significativo en la iteración {len(added_nodes) + 1}. Deteniendo.")
                        break
                        
        return added_nodes

# ----------------------------------------------------------------------
#                                                 FUNCIONES DE GUARDADO
# ----------------------------------------------------------------------
def guardar_genes_semilla_conectados(seed_genes_valid, output_file):
        """
        Guarda la lista de genes semilla que realmente tienen conexiones
        en la red filtrada (genes válidos).
        """
        data = []
        for hugo_symbol in seed_genes_valid:
                data.append({'HUGO_Symbol': hugo_symbol})

        results_df = pd.DataFrame(data)
        results_df.to_csv(output_file, sep="\t", index=False)
        print(f"Genes semilla CONECTADOS guardados en: {output_file} ({len(seed_genes_valid)} genes)")


def analizar_y_guardar_genes_aislados(seed_genes_hugo, seed_genes_valid, umbral, output_file):
        """
        Identifica, imprime y guarda los genes semilla de la lista inicial 
        que NO están en la lista de genes semilla válidos (es decir, no están
    en el grafo filtrado por umbral).
        """
        
        valid_seeds_set = set(seed_genes_valid)
        # Semillas Iniciales - Semillas Válidas = Semillas Aisladas
        isolated_seeds = set(seed_genes_hugo) - valid_seeds_set
        
        if not isolated_seeds:
                print("\nTodos los genes semilla pudieron ser validados (ninguno aislado de la red).")
                return

        # --- CONTROL DE IMPRESIÓN ---
        print(f"\nCONTROL: {len(isolated_seeds)} Genes Semilla AISLADOS de la red principal:")
        for gene in sorted(list(isolated_seeds)):
                print(f"   -> {gene}")
        print("-----------------------------------------------------")

        # --- GUARDAR EN ARCHIVO ---
        data = []
        for hugo_symbol in isolated_seeds:
                data.append({
                        'HUGO_Symbol': hugo_symbol, 
                        'Tipo': 'Seed_Gene_Aislado',
                        'Comentario': f'No conectado a la red con score >= {umbral}'
                })
                
        results_df = pd.DataFrame(data)
        results_df.to_csv(output_file, sep="\t", index=False)
        print(f"{len(isolated_seeds)} genes semilla aislados guardados en: {output_file}")


def guardar_resultados(seed_genes_hugo, diamond_genes_hugo, archivo_salida):
        """
        Guarda solo los genes añadidos por DIAMOnD en un archivo TSV. 
        """
        
        seed_set = set(seed_genes_hugo)
        diamond_only_genes = [gene for gene in diamond_genes_hugo if gene not in seed_set]

        data = []
        for hugo_symbol in diamond_only_genes:
                data.append({
                        'HUGO_Symbol': hugo_symbol
                })

        results_df = pd.DataFrame(data)
        results_df.to_csv(archivo_salida, sep="\t", index=False, header=False)
        print(f"\nResultados de DIAMOnD guardados en: {archivo_salida} ({len(results_df)} genes)")

def graficar_red_enriquecida(G, seed_genes_valid, diamond_genes_hugo, output_image_file):
        """
    Genera la visualización de la red enriquecida. Solo incluye los genes 
    semilla válidos (conectados).
    """
        # Usamos SOLO los genes semilla válidos y los genes añadidos por DIAMOnD
        all_nodes = set(seed_genes_valid) | set(diamond_genes_hugo)
        
        # 1. Crear el subgrafo
        subgraph = G.subgraph(all_nodes)
        
        if len(subgraph.nodes) < 2:
                print("Advertencia: Grafo enriquecido demasiado pequeño para dibujar.")
                return

        # 2. Calcular posiciones (solo para los nodos del subgrafo)
        pos = None
        try:
                # Aumentar iteraciones para estabilidad
                pos = nx.spring_layout(subgraph, k=0.15, iterations=50) 
        except Exception as e:
                print(f"Error al calcular el layout: {e}. No se puede dibujar el grafo.")
                return None 
        
        plt.figure(figsize=(18, 18))

        # 3. FILTRADO: Nos aseguramos de que solo dibujamos lo que está en 'pos'
        pos_keys = set(pos.keys())
        
        # Nodos de semillas (Azul) - presentes en 'pos'
        seed_genes_presentes = list(set(seed_genes_valid) & pos_keys)
        
        # Nodos DIAMOnD (Naranja) - añadidos, que no son semilla Y presentes en 'pos'
        diamond_only = list((set(diamond_genes_hugo) - set(seed_genes_valid)) & pos_keys)
        
        # --- DIBUJADO ---

        nx.draw_networkx_nodes(subgraph, pos, 
                           nodelist=seed_genes_presentes, 
                           node_color='blue', 
                           node_size=800, 
                           label="Seed Genes Conectados", 
                           alpha=0.8)

        nx.draw_networkx_nodes(subgraph, pos, 
                           nodelist=diamond_only, 
                           node_color='orange', 
                           node_size=500, 
                           label="DIAMOnD Added", 
                           alpha=0.7)
        
        # Enlaces (Edges)
        nx.draw_networkx_edges(subgraph, pos, alpha=0.4, edge_color='gray')
        
        # Etiquetas (Solo para nodos con posición)
        node_labels = {n: n for n in pos.keys()} 
                
        nx.draw_networkx_labels(subgraph, pos, labels=node_labels, font_size=8, font_weight='bold')

        plt.legend(loc="upper left", markerscale=0.7)
        total_nodes = len(seed_genes_presentes) + len(diamond_only)
        plt.title(f"Red Enriquecida con DIAMOnD ({total_nodes} Nodos, Semillas Conectadas: {len(seed_genes_presentes)})", fontsize=16)
        plt.axis('off')
        
        try:
                plt.savefig(output_image_file, format='png', dpi=300, bbox_inches='tight')
                print(f"Grafo de DIAMOnD guardado como imagen en: {output_image_file}")
        except Exception as e:
                print(f"Error al guardar la imagen: {e}")
                
        plt.close()


# ----------------------------------------------------------------------
#                                                         MAIN CLI
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
        
        print(f"--- Iniciando Propagación DIAMOnD para {nodos_añadidos} Nodos ---")
        
        # 2. Creación de la Carpeta de Resultados
        # Se usa sys.argv[0] para obtener la ruta del script
        try:
                script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
        except IndexError:
                # Fallback si se ejecuta interactivamente
                script_dir = os.path.dirname(os.path.abspath(__file__))
        
        project_root = os.path.dirname(script_dir) 
        RESULTS_DIR = os.path.join(project_root, "results")

        if not os.path.exists(RESULTS_DIR):
                os.makedirs(RESULTS_DIR)
                print(f"Carpeta de resultados creada: {RESULTS_DIR}")

        output_tsv_path = os.path.join(RESULTS_DIR, args.output)
        output_plot_path = os.path.join(RESULTS_DIR, args.plot)
        # RUTA DEL ARCHIVO DE GENES AISLADOS
        output_isolated_path = os.path.join(RESULTS_DIR, 'isolated_seed_genes.tsv')
        # NUEVA RUTA DEL ARCHIVO DE GENES CONECTADOS
        output_connected_path = os.path.join(RESULTS_DIR, 'connected_seed_genes.tsv')

        ## 3. CARGA DE DATOS
        genes_semilla_hugo = importar_genes(args.seed_file)
        
        if not genes_semilla_hugo:
                return
                
        # Carga y filtrado (usando el alto umbral definido)
        interacciones = cargar_red_conocida(args.input)
        if interacciones is None or interacciones.empty:
                return

        ## 4. CONSTRUIR GRAFO y obtener semillas VÁLIDAS
        # La función construir_red devuelve el grafo y la lista de semillas válidas (conectadas)
        red, genes_semilla_valid = construir_red(interacciones, genes_semilla_hugo)
        
        # Generar el archivo de genes conectados
        guardar_genes_semilla_conectados(genes_semilla_valid, output_connected_path)

        # Análisis temprano si no hay semillas válidas
        if not genes_semilla_valid:
                print("Ningún gen semilla está conectado a la red con el umbral especificado. Abortando DIAMOnD.")
                analizar_y_guardar_genes_aislados(genes_semilla_hugo, genes_semilla_valid, UMBRAL_SCORE, output_isolated_path)
                return
        
        # Determinamos el número real de nodos a añadir
        n = min(nodos_añadidos, len(red.nodes) - len(genes_semilla_valid))
        if n <= 0:
                print("No hay nodos para añadir o la red es muy pequeña respecto al set de semillas válidas.")
                # Aún analizamos y guardamos los aislados para dejar constancia
                analizar_y_guardar_genes_aislados(genes_semilla_hugo, genes_semilla_valid, UMBRAL_SCORE, output_isolated_path)
                return

        ## 5. EJECUTAR DIAMOnD (Usando solo las semillas VÁLIDAS)
        print(f"\n--- Ejecutando DIAMOnD para añadir {n} nodos (Cluster inicial: {len(genes_semilla_valid)} genes conectados) ---")
        # Pasamos solo los genes válidos a la función DIAMOnD
        diamond_genes = diamond_iteration_of_first_X_nodes(red, genes_semilla_valid, n)

        ## 6. GUARDAR Y GRAFICAR RESULTADOS
        guardar_resultados(genes_semilla_hugo, diamond_genes, output_tsv_path)
        
        # Graficar, usando solo los genes VÁLIDOS y los añadidos por DIAMOnD
        graficar_red_enriquecida(red, genes_semilla_valid, diamond_genes, output_plot_path)

        ## 7. GUARDAR GENES AISLADOS Y MOSTRAR CONTROL
        analizar_y_guardar_genes_aislados(genes_semilla_hugo, genes_semilla_valid, UMBRAL_SCORE, output_isolated_path)


if __name__ == '__main__':
        main()