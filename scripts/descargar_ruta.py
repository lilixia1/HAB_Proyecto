import gseapy as gp
import os

# 1. Elige la librería de rutas
# Puedes probar con: "KEGG_2021_Human", "Reactome_2022", "KEGG_2019_Human", etc.
LIBRARY = "KEGG_2021_Human"

# 2. Nombre (o parte del nombre) de la ruta que quieres
BUSCADA = "apoptosis"  

# 3. Descargar todas las rutas de esa librería
genesets = gp.get_library(name=LIBRARY, organism="Human")

ruta_encontrada = None
genes_ruta = None

for nombre, genes in genesets.items():
    if BUSCADA.lower() in nombre.lower():
        ruta_encontrada = nombre
        genes_ruta = genes
        break

if ruta_encontrada is None:
    raise ValueError(f"No se encontró ninguna ruta que contenga: {BUSCADA}")

print(f"Ruta encontrada: {ruta_encontrada}")
print(f"Número de genes: {len(genes_ruta)}")

# 4. Guardar a un txt (un gen por línea)
# Ruta base del proyecto (un nivel por encima del script actual)
basedir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Carpeta donde quieres guardar los datos
datadir = os.path.join(basedir, "data")

# Crear la carpeta si no existe
os.makedirs(datadir, exist_ok=True)

output_file = "genes_ruta.txt"
output_path = os.path.join(datadir, output_file)

with open(output_path, "w") as f:
    for g in genes_ruta:
        f.write(g + "\n")

print(f"Genes guardados en: data/{output_path}")

