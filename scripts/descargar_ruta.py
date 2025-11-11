import gseapy as gp

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
output_file = "genes_ruta.txt"
with open(output_file, "w") as f:
    for g in genes_ruta:
        f.write(g + "\n")

print(f"Genes guardados en: {output_file}")

