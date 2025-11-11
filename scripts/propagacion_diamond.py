# ----------------------------------------------------------------------
#                          LIBRERIAS Y SETUP
# ----------------------------------------------------------------------
import os
import argparse
import pandas as pd
from mygene import MyGeneInfo
from tqdm import tqdm
import networkx as nx
import scipy.stats
import matplotlib.pyplot as plt
import numpy as np
import scipy.special
import requests
import operator
import importlib
import subprocess
import sys


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